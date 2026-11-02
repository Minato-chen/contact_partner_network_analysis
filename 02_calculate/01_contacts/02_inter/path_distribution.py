"""Shortest-path distributions derived from saved chain-pair contact counts."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


def pair_distribution(adjacency):
    """Fractions over unordered node pairs; index 0 means unreachable.

    Index k>0 means shortest path k. Isolates remain nodes. Graphs with
    fewer than two nodes have no pairs and return NaNs, not zero fractions.
    """
    n = len(adjacency)
    if n < 2:
        return np.full(max(n, 1), np.nan)
    distances = shortest_path(csr_matrix(adjacency), directed=False, unweighted=True)
    distances = distances[np.triu_indices(n, 1)]
    reachable = np.isfinite(distances)
    counts = np.bincount(distances[reachable].astype(int), minlength=n).astype(float)
    counts[0] = (~reachable).sum()
    return counts / len(distances)


def summarize_paths(source, thresholds, block_frames):
    """Stream saved edge chunks; return frame-equal means and block fractions.

    Never pool pairs across frames: frames with more dense chains do not
    receive extra weight. Infinity in the tables denotes unreachable pairs.
    """
    source = Path(source)
    sums, valid = {}, {}
    n_nodes = None
    for chunk in sorted((source / "edges").glob("*.npz")):
        with np.load(chunk) as saved:
            frames = saved["frame"]
            dense_masks = saved["dense_chains"]
            if n_nodes is None:
                n_nodes = dense_masks.shape[1]
            if dense_masks.shape[1] != n_nodes:
                raise ValueError("Chain count differs between saved edge chunks")
            pair_i, pair_j = saved["pair_i"], saved["pair_j"]
            for frame, dense, contacts in zip(frames, dense_masks, saved["contact_counts"]):
                matrix = np.zeros((n_nodes, n_nodes), dtype=np.int32)
                matrix[pair_i, pair_j] = contacts
                matrix[pair_j, pair_i] = contacts
                matrix = matrix[np.ix_(dense, dense)]
                block = int(frame) // block_frames
                for threshold in thresholds:
                    key = (threshold, block)
                    if key not in sums:
                        sums[key] = np.zeros(max(n_nodes, 1))
                        valid[key] = 0
                    if len(matrix) < 2:
                        continue
                    fraction = pair_distribution(matrix >= threshold)
                    sums[key][: len(fraction)] += fraction
                    valid[key] += 1
    if n_nodes is None:
        raise ValueError(f"No saved network edge chunks in {source}")
    rows = []
    for (threshold, block), total in sums.items():
        for step, value in enumerate(total):
            rows.append(
                dict(
                    threshold=threshold,
                    block=block,
                    path_length=step if step else np.inf,
                    fraction=value / valid[(threshold, block)]
                    if valid[(threshold, block)]
                    else np.nan,
                    valid_frames=valid[(threshold, block)],
                    total_fraction=value,
                )
            )
    blocks = pd.DataFrame(rows)
    keys = ["threshold", "path_length"]
    totals = blocks.groupby(keys, sort=True)[["total_fraction", "valid_frames"]].sum()
    summary = totals.assign(mean=totals.total_fraction / totals.valid_frames.replace(0, np.nan))
    summary = summary.rename(columns={"valid_frames": "count_frames"}).drop(
        columns="total_fraction"
    )
    statistics = blocks.groupby(keys).fraction.agg(["std", "count"])
    summary["std_blocks"] = statistics["std"]
    summary["count_blocks"] = statistics["count"]
    summary["sem_blocks"] = statistics["std"] / np.sqrt(statistics["count"].replace(0, np.nan))
    return summary.reset_index(), blocks.drop(columns="total_fraction")


def display_bins(blocks, max_steps=10):
    """Fold long paths for plotting; recompute tail SD from block fractions."""
    if isinstance(max_steps, bool) or not isinstance(max_steps, int) or max_steps < 1:
        raise ValueError("max_steps must be a positive integer")
    rows = []
    for (simulation, threshold), group in blocks.groupby(["simulation_id", "threshold"]):
        categories = [*range(1, max_steps + 1), max_steps + 1, 0]
        for category in categories:
            values, counts = [], []
            for _, block in group.groupby("block"):
                count = int(block.valid_frames.iloc[0])
                if category == 0:
                    selected = np.isinf(block.path_length)
                elif category == max_steps + 1:
                    selected = np.isfinite(block.path_length) & block.path_length.gt(max_steps)
                else:
                    selected = block.path_length.eq(category)
                values.append(float(block.loc[selected, "fraction"].sum()) if count else np.nan)
                counts.append(count)
            values, counts = np.asarray(values), np.asarray(counts)
            valid = counts > 0
            mean = (
                float(np.sum(values[valid] * counts[valid]) / counts[valid].sum())
                if valid.any()
                else np.nan
            )
            sd = float(values[valid].std(ddof=1)) if valid.sum() > 1 else np.nan
            rows.append(
                dict(
                    simulation_id=simulation,
                    threshold=threshold,
                    bin=category,
                    label="Unreachable"
                    if category == 0
                    else f">{max_steps}"
                    if category == max_steps + 1
                    else str(category),
                    mean=mean,
                    std_blocks=sd,
                )
            )
    return pd.DataFrame(rows)
