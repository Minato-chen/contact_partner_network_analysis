#!/usr/bin/env python3
"""Exploratory persistence of dense-phase chain-pair contacts (isolated test)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd

from dependencies.partner_analysis import _chain_layout, _minimum_image, frame_chain_contacts


VARIANTS = {"A18": 11.25, "B01": 8.91, "B10": 9.81}
COLORS = {"A18": "#888888", "B01": "#3975b8", "B10": "#d06c32"}


def load_edges(pdb: Path, dcd: Path, z_den: float, max_frames: int,
               cutoff: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    traj = md.load_dcd(str(dcd), top=str(pdb))[-max_frames:]
    atom_chain, chain_atoms = _chain_layout(traj.topology)
    n = len(chain_atoms)
    pair_i, pair_j = np.triu_indices(n, 1)
    edges = np.zeros((traj.n_frames, len(pair_i)), dtype=np.uint16)
    dense_pairs = np.zeros_like(edges, dtype=bool)
    for frame, (xyz, box) in enumerate(zip(traj.xyz, traj.unitcell_lengths)):
        com = xyz[chain_atoms].mean(axis=1)
        z = _minimum_image(com[:, 2] - box[2] / 2, box[2])
        dense = np.abs(z) < z_den
        counts = frame_chain_contacts(xyz, box, atom_chain, n, cutoff)
        edges[frame] = counts[pair_i, pair_j]
        dense_pairs[frame] = dense[pair_i] & dense[pair_j]
        if (frame + 1) % 50 == 0:
            print(f"  processed {frame + 1}/{traj.n_frames} frames", flush=True)
    return edges, dense_pairs


def summarize(variant: str, edges: np.ndarray, dense_pairs: np.ndarray,
              thresholds: tuple[int, ...], max_lag: int,
              block_frames: int) -> pd.DataFrame:
    # Each time origin is assigned to a contiguous block. Pooled pair origins
    # are descriptive; blocks, not individual pairs, are variability units.
    n_frames = len(edges)
    rows = []
    for threshold in thresholds:
        active = (edges >= threshold) & dense_pairs
        for block, start in enumerate(range(0, n_frames, block_frames)):
            stop = min(start + block_frames, n_frames)
            for lag in range(max_lag + 1):
                last_start = min(stop, n_frames - lag)
                origins = active[start:last_start]
                if len(origins) == 0:
                    continue
                at_risk = origins.sum()
                if not at_risk:
                    continue
                continuous = origins.copy()
                dense_continuous = dense_pairs[start:last_start].copy()
                for offset in range(1, lag + 1):
                    continuous &= active[start + offset:last_start + offset]
                    dense_continuous &= dense_pairs[start + offset:last_start + offset]
                # Primary: leaving the dense core counts as loss of the edge.
                # Sensitivity: only origins where both chains stay dense.
                resident_origins = origins & dense_continuous
                resident_n = resident_origins.sum()
                rows.append({
                    "variant": variant, "threshold": threshold, "block": block,
                    "lag": lag, "at_risk": int(at_risk),
                    "survived": int(continuous.sum()),
                    "resident_at_risk": int(resident_n),
                    "resident_survived": int((continuous & resident_origins).sum()),
                    "dense_exit": int(at_risk - resident_n),
                })
    return pd.DataFrame(rows)


def plot(rows: pd.DataFrame, output: Path) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9,
        "axes.spines.right": False, "axes.spines.top": False,
        "legend.frameon": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), sharey=True,
                             constrained_layout=True)
    for ax, threshold in zip(axes, (1, 5)):
        for variant in VARIANTS:
            subset = rows[(rows.variant == variant) & (rows.threshold == threshold)]
            pooled = subset.groupby("lag")[["at_risk", "survived"]].sum()
            survival = pooled.survived / pooled.at_risk
            ax.plot(survival.index, survival, marker="o", markersize=3,
                    color=COLORS[variant], label=variant)
        ax.set(xlabel="Lag (sampled frames)", title=f"≥{threshold} residue contacts",
               ylim=(0, 1.03), xlim=(0, rows.lag.max()))
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Continuous chain-pair contact survival")
    axes[1].legend(title="Variant")
    output.parent.mkdir(parents=True, exist_ok=True)
    for extension, kwargs in ((".png", {"dpi": 250}), (".svg", {}), (".pdf", {})):
        fig.savefig(output.with_suffix(extension), bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=Path("."))
    parser.add_argument("--centered-root", type=Path,
                        default=Path("isolated_partner_analysis/work_local_test"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("isolated_partner_analysis/results_contact_persistence_T260_local"))
    parser.add_argument("--max-frames", type=int, default=750)
    parser.add_argument("--max-lag", type=int, default=20)
    parser.add_argument("--block-frames", type=int, default=50)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for variant, z_den in VARIANTS.items():
        print(f"{variant}:", flush=True)
        pdb = args.data_root / f"tau_k18_129aa_{variant}_T260_5/top.pdb"
        dcd = args.centered_root / f"{variant}_T260_centered.dcd"
        edges, dense = load_edges(pdb, dcd, z_den, args.max_frames)
        np.savez_compressed(args.output_dir / f"{variant}_edge_arrays.npz",
                            contact_counts=edges, both_dense=dense)
        rows = summarize(variant, edges, dense, (1, 5), args.max_lag,
                         args.block_frames)
        all_rows.append(rows)
    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(args.output_dir / "persistence_blocks.csv", index=False)
    pooled = combined.groupby(["variant", "threshold", "lag"], as_index=False)[
        ["at_risk", "survived", "resident_at_risk", "resident_survived", "dense_exit"]
    ].sum()
    pooled["survival"] = pooled.survived / pooled.at_risk
    pooled["resident_survival"] = pooled.resident_survived / pooled.resident_at_risk
    pooled["dense_exit_fraction"] = pooled.dense_exit / pooled.at_risk
    pooled.to_csv(args.output_dir / "persistence_summary.csv", index=False)
    plot(combined, args.output_dir / "contact_persistence_local")
    print(pooled[pooled.lag.isin((1, 5, 10, 20))].to_string(index=False))


if __name__ == "__main__":
    main()
