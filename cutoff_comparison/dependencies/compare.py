"""Finite distance sensitivity study; one candidate search per frame."""
from pathlib import Path
import numpy as np
import pandas as pd
from common import results as files

OBJECTS = ("central_chain", "dense_chains")


def cutoff_name(value):
    return f"{value:.1f}".replace(".", "p") + "nm"


def aggregate_distribution(sums, counts, parameters):
    """Return direct frame-equal means and block fractions with correct weights."""
    rows = []
    for key, total in sums.items():
        obj, threshold, block = key
        count = counts[key]
        for position, value in enumerate(total):
            rows.append(
                dict(
                    object=obj,
                    threshold=threshold,
                    block=block,
                    position=position,
                    fraction=float(value / count) if count else np.nan,
                    valid_frames=count,
                    total_fraction=float(value),
                )
            )
    blocks = pd.DataFrame(rows)
    keys = ["object", "threshold", "position"]
    summary = blocks.groupby(keys)[["total_fraction", "valid_frames"]].sum()
    summary["mean"] = summary.total_fraction / summary.valid_frames.replace(0, np.nan)
    stats = blocks.groupby(keys).fraction.agg(["std", "count"])
    summary["std_blocks"] = stats["std"]
    summary["count_blocks"] = stats["count"]
    summary["sem_blocks"] = stats["std"] / np.sqrt(stats["count"].replace(0, np.nan))
    summary = summary.drop(columns="total_fraction").rename(
        columns={"valid_frames": "count_frames"}
    )
    return summary.reset_index(), blocks.drop(columns="total_fraction")


def compare(simulation, centered_root, phase, targets, selection_output, parameters, io):
    cutoffs = parameters["cutoffs_nm"]
    width = float(files.read_json(phase / "boundary.json")["dense_half_width_nm"])
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Invalid measured dense boundary")
    frame_rows = {r: [] for r in cutoffs}
    profile_sums = {r: {} for r in cutoffs}
    profile_counts = {r: {} for r in cutoffs}
    partner_sums = {r: {} for r in cutoffs}
    partner_counts = {r: {} for r in cutoffs}
    selection_buffer = []
    for info, xyz, box, chains, atom_chain, residue, com in io.centered_frames(
        simulation, centered_root, parameters["chunk_frames"]
    ):
        n, n_res = chains.shape
        if width >= box[2] / 2:
            raise ValueError("Dense width exceeds half the box")
        z = io.minimum_image(com[:, 2] - box[2] / 2, box[2])
        dense = np.abs(z) < width
        eligible = np.flatnonzero(dense)
        central = int(eligible[np.argmin(np.abs(z[eligible]))]) if len(eligible) else -1
        selected = dict(
            dense_chains=eligible,
            central_chain=np.array([central], int) if central >= 0 else np.array([], int),
        )
        selection_buffer.extend(
            info | dict(chain_id=c, dense=bool(dense[c]), central=c == central) for c in range(n)
        )
        block = info["frame"] // parameters["block_frames"]
        pairs = io.close_pairs(xyz, box, max(cutoffs))
        wrapped = np.mod(xyz, box)
        delta = io.minimum_image(
            wrapped[pairs[:, 0]].astype(float) - wrapped[pairs[:, 1]].astype(float), box
        )
        squared = np.einsum("ij,ij->i", delta, delta)
        ci, cj = atom_chain[pairs[:, 0]], atom_chain[pairs[:, 1]]
        ri, rj = residue[pairs[:, 0]], residue[pairs[:, 1]]
        same = ci == cj
        separated = np.abs(ri - rj) > parameters["exclude_sequence_separation"]
        previous_intra = previous_inter = None
        for cutoff in cutoffs:
            close = squared < cutoff**2
            internal = close & same & separated
            intrachain = np.bincount(ci[internal], minlength=n)
            degree = np.zeros((n, n_res), dtype=np.int32)
            np.add.at(degree, (ci[internal], ri[internal]), 1)
            np.add.at(degree, (ci[internal], rj[internal]), 1)
            external = close & ~same
            a, b = np.minimum(ci[external], cj[external]), np.maximum(ci[external], cj[external])
            upper = np.bincount(a * n + b, minlength=n * n).reshape(n, n)
            inter = upper + upper.T
            if previous_intra is not None:
                if (intrachain < previous_intra).any() or (inter < previous_inter).any():
                    raise AssertionError("Contact counts decreased as distance increased")
            previous_intra, previous_inter = intrachain, inter
            inter[~dense, :] = 0
            inter[:, ~dense] = 0
            for obj, chosen in selected.items():
                count = len(chosen)
                profile = (degree[chosen] > 0).mean(axis=0) if count else np.full(n_res, np.nan)
                key = (obj, 0, block)
                profile_sums[cutoff].setdefault(key, np.zeros(n_res))
                profile_counts[cutoff].setdefault(key, 0)
                if count:
                    profile_sums[cutoff][key] += profile
                    profile_counts[cutoff][key] += 1
                totals = inter[chosen].sum(axis=1)
                for threshold in parameters["partner_thresholds"]:
                    partners = (inter[chosen] >= threshold).sum(axis=1)
                    frame_rows[cutoff].append(
                        info
                        | dict(
                            block=block,
                            object=obj,
                            threshold=threshold,
                            n_chains=count,
                            intra_contacts=float(intrachain[chosen].mean()) if count else np.nan,
                            inter_contacts=float(totals.mean()) if count else np.nan,
                            partner_count=float(partners.mean()) if count else np.nan,
                            participation=float(profile.mean()) if count else np.nan,
                        )
                    )
                    key = (obj, threshold, block)
                    partner_sums[cutoff].setdefault(key, np.zeros(n))
                    partner_counts[cutoff].setdefault(key, 0)
                    if count:
                        partner_sums[cutoff][key] += np.bincount(partners, minlength=n) / count
                        partner_counts[cutoff][key] += 1
        if (info["frame"] + 1) % parameters["chunk_frames"] == 0:
            io.append_csv(selection_buffer, selection_output / "chain_selection.csv.gz")
            selection_buffer.clear()
            print(f'{simulation["id"]}: {info["frame"]+1} frames × four distances', flush=True)
    if selection_buffer:
        io.append_csv(selection_buffer, selection_output / "chain_selection.csv.gz")
    if not frame_rows[cutoffs[0]]:
        raise ValueError("No frames")
    for cutoff, target in targets.items():
        frame = pd.DataFrame(frame_rows[cutoff])
        frame.to_csv(target / "frame_means.csv", index=False)
        io.summarize(
            frame,
            ["intra_contacts", "inter_contacts", "partner_count", "participation"],
            ["simulation_id", "object", "threshold"],
            target,
        )
        for filename, sums, counts in [
            ("residue_probability", profile_sums[cutoff], profile_counts[cutoff]),
            ("partner_distribution", partner_sums[cutoff], partner_counts[cutoff]),
        ]:
            summary, blocks = aggregate_distribution(sums, counts, parameters)
            name = "residue_index" if filename == "residue_probability" else "partner_count"
            if filename == "residue_probability":
                summary = summary.drop(columns="threshold")
                blocks = blocks.drop(columns="threshold")
            summary.rename(columns={"position": name}).to_csv(
                target / (filename + ".csv"), index=False
            )
            blocks.rename(columns={"position": name}).to_csv(
                target / (filename + "_blocks.csv"), index=False
            )
