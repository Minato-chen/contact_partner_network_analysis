"""Verified stage-1 inputs, bounded trajectory reads and shared numeric utilities."""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import mdtraj as md
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from common import results as files


STAGE = Path(__file__).resolve().parents[1]


def validate_parameters(parameters):
    for key in ("contact_cutoff_nm",):
        if key in parameters and (not np.isfinite(parameters[key]) or parameters[key] <= 0):
            raise ValueError(f"{key} must be finite and positive")
    for key in ("chunk_frames", "block_frames"):
        value = parameters[key]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{key} must be a positive integer")
    for key in ("partner_thresholds", "network_thresholds"):
        if key in parameters:
            values = parameters[key]
            if (
                not values
                or len(values) != len(set(values))
                or any(isinstance(v, bool) or not isinstance(v, int) or v <= 0 for v in values)
            ):
                raise ValueError(f"{key} requires distinct positive integers")
    if "exclude_sequence_separation" in parameters:
        value = parameters["exclude_sequence_separation"]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError("exclude_sequence_separation must be a nonnegative integer")


def centered_source(simulation, centered_root):
    return files.completed(Path(centered_root) / simulation["id"])


def centered_reference(simulation, centered_root):
    return files.read_json(centered_source(simulation, centered_root) / "metadata.json")


def phase_source(simulation, centered_root, results_root=files.RESULTS):
    source = files.completed(Path(results_root) / simulation["id"] / "01_basic/phase_boundary")
    saved = files.read_json(source / "metadata.json")
    if saved["upstream"] != centered_reference(simulation, centered_root):
        raise files.InputUnavailable(
            "Centered output changed; rerun 01_basic/phase_boundary.py first"
        )
    return source


def minimum_image(delta, box):
    return delta - box * np.round(delta / box)


def close_pairs(xyz, box, cutoff):
    wrapped = np.mod(xyz, box)
    pairs = cKDTree(wrapped, boxsize=box).query_pairs(cutoff, output_type="ndarray")
    if len(pairs):
        delta = minimum_image(
            wrapped[pairs[:, 0]].astype(float) - wrapped[pairs[:, 1]].astype(float), box
        )
        pairs = pairs[np.einsum("ij,ij->i", delta, delta) < cutoff**2]
    return pairs


def centered_frames(simulation, centered_root, chunk_frames):
    source = centered_source(simulation, centered_root)
    chain_lists = files.read_json(source / "chains.json")
    lengths = {len(ids) for ids in chain_lists}
    if len(lengths) != 1:
        raise ValueError("Contact maps currently require equal-length linear chains")
    chains = np.asarray(chain_lists, dtype=int)
    n_atoms = chains.size
    if not np.array_equal(np.sort(chains.ravel()), np.arange(n_atoms)):
        raise ValueError("Chain layout must partition every atom exactly once")
    atom_chain = np.empty(n_atoms, dtype=int)
    atom_residue = np.empty(n_atoms, dtype=int)
    for index, ids in enumerate(chains):
        atom_chain[ids] = index
        atom_residue[ids] = np.arange(len(ids))
    frames = pd.read_csv(source / "frame_index.csv")
    expected = files.read_json(source / "metadata.json")["nframes_written"]
    if len(frames) != expected or not np.array_equal(frames["frame"], np.arange(expected)):
        raise ValueError("Stage-1 frame index is incomplete or unordered")
    fi = 0
    # ASCII aliases preserve compatibility with native DCD libraries.
    with TemporaryDirectory(prefix="contact_read_") as temporary:
        alias = Path(temporary) / "trajectory.dcd"
        topology = Path(temporary) / "topology.pdb"
        alias.symlink_to((source / "centered.dcd").resolve())
        topology.symlink_to((source / "topology.pdb").resolve())
        for chunk in md.iterload(str(alias), top=str(topology), chunk=chunk_frames):
            if chunk.n_atoms != n_atoms:
                raise ValueError("Topology and stage-1 chain layout disagree")
            if chunk.unitcell_angles is None or not np.allclose(
                chunk.unitcell_angles, 90, atol=1e-4
            ):
                raise ValueError("Contact analysis requires orthorhombic boxes")
            if chunk.unitcell_lengths is None or not np.all(
                np.isfinite(chunk.unitcell_lengths)
            ):
                raise ValueError("Missing or nonfinite periodic box")
            for xyz, box in zip(chunk.xyz, chunk.unitcell_lengths):
                if np.any(box <= 0) or not np.all(np.isfinite(xyz)):
                    raise ValueError("Invalid coordinates or periodic box")
                com = xyz[chains].mean(axis=1)
                if fi >= expected:
                    raise ValueError("More trajectory frames than upstream frame records")
                row = frames.iloc[fi]
                info = dict(
                    simulation_id=simulation["id"],
                    frame=fi,
                    source_frame=int(row.source_frame),
                    time_ps=float(row.time_ps),
                )
                yield info, xyz, box, chains, atom_chain, atom_residue, com
                fi += 1
    if fi != expected:
        raise ValueError(f"Expected {expected} centered frames, read {fi}")


def trajectory(simulation, centered_root, chunk_frames, results_root=files.RESULTS):
    phase = phase_source(simulation, centered_root, results_root)
    width = float(files.read_json(phase / "boundary.json")["dense_half_width_nm"])
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Invalid upstream dense half width")
    for info, xyz, box, chains, atom_chain, residue, com in centered_frames(
        simulation, centered_root, chunk_frames
    ):
        if width >= box[2] / 2:
            raise ValueError("Upstream dense half width exceeds the current box")
        z = minimum_image(com[:, 2] - box[2] / 2, box[2])
        yield info, xyz, box, chains, atom_chain, residue, com, np.abs(z) < width


def append_csv(rows, path):
    frame = pd.DataFrame.from_records(rows)
    frame.to_csv(path, mode="a", header=not Path(path).exists(), index=False)


def summarize(frame, metrics, keys, output):
    """Direct frame-equal means; SD of frame means and SEM of block means."""
    per_block = frame.groupby(keys + ["block"], dropna=False)[metrics].mean().reset_index()
    summary = frame.groupby(keys, dropna=False)[metrics].agg(["mean", "std", "count"])
    summary.columns = [f"{metric}_{stat}_frames" for metric, stat in summary.columns]
    block = per_block.groupby(keys, dropna=False)[metrics].agg(["std", "count"])
    for metric in metrics:
        summary[f"{metric}_std_blocks"] = block[(metric, "std")]
        summary[f"{metric}_count_blocks"] = block[(metric, "count")]
        summary[f"{metric}_sem_blocks"] = block[(metric, "std")] / np.sqrt(
            block[(metric, "count")]
        )
    summary.reset_index().to_csv(output / "summary.csv", index=False)
    per_block.to_csv(output / "block_means.csv", index=False)
