#!/usr/bin/env python3
"""Chain-partner analysis for Tau K18 slab trajectories.

This standalone analysis is intentionally separate from ``ana_v4``: it tests one
specific hypothesis (many contacts concentrated on few partners) without
changing the established production summaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


THRESHOLDS = (1, 5, 10)


def _chain_layout(topology: md.Topology) -> tuple[np.ndarray, np.ndarray]:
    """Return atom-to-chain ids and rectangular chain atom indices."""
    indices = [np.asarray([a.index for a in chain.atoms], dtype=int)
               for chain in topology.chains]
    lengths = {len(x) for x in indices}
    if len(lengths) != 1:
        raise ValueError("This analysis requires equal-length chains")
    chain_atoms = np.stack(indices)
    atom_chain = np.empty(topology.n_atoms, dtype=int)
    for ci, idx in enumerate(chain_atoms):
        atom_chain[idx] = ci
    return atom_chain, chain_atoms


def _minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def frame_chain_contacts(
    xyz: np.ndarray,
    box: np.ndarray,
    atom_chain: np.ndarray,
    n_chains: int,
    cutoff_nm: float,
) -> np.ndarray:
    """Count residue contacts for every unordered chain pair in one frame."""
    wrapped = np.mod(xyz, box)
    atom_pairs = cKDTree(wrapped, boxsize=box).query_pairs(
        cutoff_nm, output_type="ndarray"
    )
    # query_pairs includes the radius boundary; enforce the strict d < cutoff
    # definition explicitly, using minimum-image distances in nanometres.
    if atom_pairs.size:
        delta = (wrapped[atom_pairs[:, 0]].astype(np.float64)
                 - wrapped[atom_pairs[:, 1]].astype(np.float64))
        delta = _minimum_image(delta, np.asarray(box, dtype=np.float64))
        squared_distance = np.einsum("ij,ij->i", delta, delta)
        atom_pairs = atom_pairs[squared_distance < cutoff_nm ** 2]
    counts = np.zeros((n_chains, n_chains), dtype=np.int32)
    if atom_pairs.size == 0:
        return counts
    ci = atom_chain[atom_pairs[:, 0]]
    cj = atom_chain[atom_pairs[:, 1]]
    keep = ci != cj
    lo = np.minimum(ci[keep], cj[keep])
    hi = np.maximum(ci[keep], cj[keep])
    pair_id = lo * n_chains + hi
    bincount = np.bincount(pair_id, minlength=n_chains * n_chains)
    upper = bincount.reshape(n_chains, n_chains)
    counts += upper + upper.T
    return counts


def analyze_trajectory(
    variant: str,
    pdb: Path,
    dcd: Path,
    z_den_nm: float,
    z_dil_nm: float,
    start: int = 750,
    step: int = 3,
    cutoff_nm: float = 1.0,
    local_radius_nm: float = 5.0,
    block_frames: int = 50,
    max_frames: int | None = None,
) -> pd.DataFrame:
    """Return one row per dense chain, frame, and partner threshold."""
    traj = md.load_dcd(str(dcd), top=str(pdb), stride=step)[start // step:]
    if max_frames is not None:
        traj = traj[:max_frames]
    atom_chain, chain_atoms = _chain_layout(traj.topology)
    n_chains = len(chain_atoms)
    records: list[dict] = []

    for fi, (xyz, box) in enumerate(zip(traj.xyz, traj.unitcell_lengths)):
        chain_xyz = xyz[chain_atoms]
        com = chain_xyz.mean(axis=1)
        centered_z = _minimum_image(com[:, 2] - box[2] / 2, box[2])
        dense = np.abs(centered_z) < z_den_nm
        dilute = np.abs(centered_z) > z_dil_nm

        pair_counts = frame_chain_contacts(
            xyz, box, atom_chain, n_chains, cutoff_nm
        )
        # The hypothesis concerns the dense network: exclude contacts to chains
        # outside the dense core rather than mixing phase and interface effects.
        pair_counts[~dense, :] = 0
        pair_counts[:, ~dense] = 0

        dcom = _minimum_image(com[:, None, :] - com[None, :, :], box)
        com_dist = np.linalg.norm(dcom, axis=-1)
        local_neighbors = ((com_dist < local_radius_nm) & dense[:, None]
                           & dense[None, :])
        np.fill_diagonal(local_neighbors, False)

        for ci in np.flatnonzero(dense):
            weights = pair_counts[ci]
            total = int(weights.sum())
            for threshold in THRESHOLDS:
                partners = int(np.count_nonzero(weights >= threshold))
                records.append({
                    "variant": variant,
                    "frame": fi,
                    "block": fi // block_frames,
                    "chain": int(ci),
                    "phase": "dense",
                    "is_dilute": bool(dilute[ci]),
                    "threshold": threshold,
                    "total_contacts": total,
                    "partner_count": partners,
                    "contacts_per_partner": total / partners if partners else np.nan,
                    "local_neighbors_5nm": int(local_neighbors[ci].sum()),
                })
    return pd.DataFrame.from_records(records)


METRICS = ["total_contacts", "partner_count", "contacts_per_partner",
           "local_neighbors_5nm"]


def frame_partner_means(rows: pd.DataFrame) -> pd.DataFrame:
    """One equally weighted observation per variant, threshold and frame.

    contacts_per_partner remains the mean of valid per-chain ratios; zero-
    partner chains have NaN ratios and are excluded only from that metric.
    """
    keys = ["variant", "threshold", "frame", "block"]
    if rows.empty:
        raise ValueError("Cannot summarize empty chain records")
    if rows.groupby(keys[:-1])["block"].nunique().ne(1).any():
        raise ValueError("Each frame must belong to exactly one time block")
    grouped = rows.groupby(keys, sort=True)
    frame = grouped[METRICS].mean()
    frame["n_dense"] = grouped.size()
    frame["n_valid_contacts_per_partner"] = grouped["contacts_per_partner"].count()
    return frame.reset_index()


def summarize_blocks(rows: pd.DataFrame) -> pd.DataFrame:
    """Frame-equal means and frame SD, with block-based uncertainty estimates.

    Legacy *_std and *_count retain their meanings (block SD and valid block
    count). *_mean is always computed directly from frame means, including
    when blocks have unequal lengths. Block SEM assumes sufficiently weak
    correlation between block means; block length must be checked separately.
    """
    keys = ["variant", "threshold"]
    frame = frame_partner_means(rows)
    per_block = frame.groupby(keys + ["block"], as_index=False)[METRICS].mean()
    summary = per_block.groupby(keys)[METRICS].agg(["mean", "std", "count"])
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    frames = frame.groupby(keys)
    for metric in METRICS:
        summary[f"{metric}_mean"] = frames[metric].mean()
        summary[f"{metric}_std_frames"] = frames[metric].std()
        summary[f"{metric}_count_frames"] = frames[metric].count()
        summary[f"{metric}_std_blocks"] = summary[f"{metric}_std"]
        summary[f"{metric}_sem_blocks"] = (
            summary[f"{metric}_std_blocks"] /
            np.sqrt(summary[f"{metric}_count"].clip(lower=1))
        )
    summary["n_frames"] = frames.size()
    summary["n_blocks"] = frames["block"].nunique()
    block_sizes = frame.groupby(keys + ["block"]).size()
    summary["block_frames_min"] = block_sizes.groupby(level=keys).min()
    summary["block_frames_max"] = block_sizes.groupby(level=keys).max()
    return summary.reset_index()


def plot_partner_test(rows: pd.DataFrame, output_stem: Path) -> None:
    """Quantitative-grid figure: primary test plus threshold robustness."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
        "axes.spines.right": False, "axes.spines.top": False,
        "axes.linewidth": 0.8, "legend.frameon": False,
    })
    colors = {"A18_T260": "#8a8f98", "B01_T260": "#3b73b9",
              "B10_T260": "#d07a3f"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)

    primary = rows[rows["threshold"] == 1].copy()
    # Frame averages prevent frames with more dense chains receiving more weight.
    frame = (primary.groupby(["variant", "frame"], as_index=False)
             [["total_contacts", "partner_count"]].mean())
    for variant, data in frame.groupby("variant"):
        axes[0].scatter(data["total_contacts"], data["partner_count"], s=8,
                        alpha=0.22, linewidths=0,
                        color=colors.get(variant, "#555555"), label=variant)
    axes[0].set(xlabel="Residue contacts per chain",
                ylabel="Distinct partners per chain")
    axes[0].legend(title="Variant", markerscale=1.8)

    summary = summarize_blocks(rows)
    variants = list(dict.fromkeys(rows["variant"]))
    for variant in variants:
        stat = (summary[summary["variant"] == variant]
                .sort_values("threshold"))
        axes[1].errorbar(stat["threshold"].to_numpy(), stat["partner_count_mean"],
                         yerr=stat["partner_count_sem_blocks"], marker="o", ms=4, capsize=2,
                         color=colors.get(variant, "#555555"), label=variant)
    axes[1].set(xlabel="Minimum residue contacts defining a partner",
                ylabel="Distinct partners per chain", xticks=THRESHOLDS)

    for label, ax in zip(("a", "b"), axes):
        ax.text(-0.14, 1.04, label, transform=ax.transAxes,
                fontweight="bold", fontsize=8)

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--dcd", type=Path, required=True)
    parser.add_argument("--z-den", type=float, required=True)
    parser.add_argument("--z-dil", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start", type=int, default=750)
    parser.add_argument("--step", type=int, default=3)
    parser.add_argument("--max-frames", type=int, default=None,
                        help="Local smoke test only; omit for the full server run")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    output = args.output_dir / f"{args.variant}_chain_partner_rows.csv.gz"
    summary = args.output_dir / f"{args.variant}_chain_partner_summary.csv"
    if output.exists() and summary.exists() and not args.force:
        print(f"[skip] {args.variant}: outputs already exist")
        return
    if not args.pdb.is_file():
        raise FileNotFoundError(f"Missing topology: {args.pdb}")
    if not args.dcd.is_file():
        raise FileNotFoundError(f"Missing trajectory: {args.dcd}")

    rows = analyze_trajectory(args.variant, args.pdb, args.dcd, args.z_den,
                              args.z_dil, start=args.start, step=args.step,
                              max_frames=args.max_frames)
    if rows.empty:
        raise RuntimeError(f"{args.variant}: no dense-chain samples were produced")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp.gz")
    rows.to_csv(temporary, index=False, compression="gzip")
    temporary.replace(output)
    summarize_blocks(rows).to_csv(summary, index=False)
    print(f"[ok] {args.variant}: {rows.frame.nunique()} frames -> {output}")


if __name__ == "__main__":
    main()
