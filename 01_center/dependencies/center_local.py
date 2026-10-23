#!/usr/bin/env python3
"""Two-pass slab centering and whole linear-chain reconstruction.

MDAnalysis coordinates are in Angstrom. Only fixed orthorhombic boxes and
one bead per residue, ordered linear chains (one per segment) are supported.
"""
from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import MDAnalysis as mda
import numpy as np


def _density_runs(counts: np.ndarray) -> list[np.ndarray]:
    occupied = np.flatnonzero(counts > 0)
    if not occupied.size:
        return []
    runs = list(np.split(occupied, np.flatnonzero(np.diff(occupied) > 1) + 1))
    if len(runs) > 1 and occupied[0] == 0 and occupied[-1] == len(counts) - 1:
        runs = [np.concatenate((runs[-1], runs[0]))] + runs[1:-1]
    return runs


def _largest_density_patch(z: np.ndarray, counts: np.ndarray,
                           box_length: float) -> tuple[np.ndarray, np.ndarray]:
    """Merge periodic end runs and unwrap their coordinates before averaging."""
    runs = _density_runs(counts)
    if not runs:
        return np.array([]), np.array([])
    best = max(runs, key=lambda ids: counts[ids].sum())
    positions = z[best].copy()
    wraps = np.flatnonzero(np.diff(best) < 0)
    if wraps.size:
        positions[wraps[0] + 1:] += box_length
    return positions, counts[best]


def _box(dimensions: np.ndarray) -> np.ndarray:
    if dimensions is None or not np.all(np.isfinite(dimensions)):
        raise ValueError('A finite periodic box is required')
    if np.any(dimensions[:3] <= 0) or not np.allclose(dimensions[3:], 90, atol=1e-4):
        raise ValueError('Only positive orthorhombic boxes are supported')
    return np.asarray(dimensions[:3], dtype=float)


def _chains(universe: mda.Universe) -> list[np.ndarray]:
    chains = []
    for segment in universe.segments:
        ids = segment.atoms.indices
        if len(ids) < 2 or any(len(res.atoms) != 1 for res in segment.residues):
            raise ValueError('Each segment must be an ordered linear chain, one bead per residue')
        chains.append(ids)
    if not chains or sum(map(len, chains)) != len(universe.atoms):
        raise ValueError('Segments must cover all beads')
    return chains


def _make_chains_whole(xyz: np.ndarray, chains: list[np.ndarray],
                       box: np.ndarray, max_bond_A: float) -> tuple[np.ndarray, float]:
    """Reconstruct along sequential backbone bonds, then translate whole chains."""
    result = np.asarray(xyz, dtype=float).copy()
    largest = 0.0
    for ids in chains:
        chain = result[ids].copy()
        delta = np.diff(chain, axis=0)
        delta -= box * np.round(delta / box)
        lengths = np.linalg.norm(delta, axis=1)
        largest = max(largest, float(lengths.max()))
        if not np.all(np.isfinite(lengths)) or lengths.max() > max_bond_A:
            raise ValueError('Backbone exceeds max_bond_nm; check topology order and box')
        chain[1:] = chain[0] + np.cumsum(delta, axis=0)
        # A bead may lie outside the box; the chain itself must stay whole.
        chain -= box * np.floor(chain.mean(axis=0) / box)
        result[ids] = chain
    displacement = result - xyz
    displacement -= box * np.round(displacement / box)
    if np.max(np.abs(displacement)) > 1e-3:
        raise AssertionError('Whole-chain reconstruction changed periodic bead positions')
    return result, largest


def center(pdb: Path, dcd: Path, output: Path, start: int, step: int, *,
           expected_frames: int | None = None,
           max_bond_nm: float = 1.0, stop: int | None = None) -> dict:
    """Center selected frames and restore whole chains; return sampling metadata."""
    pdb, dcd, output = map(Path, (pdb, dcd, output))
    if output.resolve() in (pdb.resolve(), dcd.resolve()):
        raise ValueError('Output must not overwrite an input')
    if start < 0 or step <= 0 or not np.isfinite(max_bond_nm) or max_bond_nm <= 0:
        raise ValueError('Require start >= 0, step > 0 and finite max_bond_nm > 0')
    if expected_frames is not None and expected_frames <= 0:
        raise ValueError('expected_frames must be positive')
    if stop is not None and stop <= start:
        raise ValueError('stop must be greater than start (exclusive stop frame)')
    universe = mda.Universe(str(pdb), str(dcd))
    try:
        n_total = len(universe.trajectory)
        if stop is not None and stop > n_total:
            raise ValueError(f'stop={stop} exceeds trajectory length {n_total}')
        frames = range(start, n_total if stop is None else stop, step)
        if not frames or (expected_frames is not None and len(frames) != expected_frames):
            raise ValueError(f'Selected {len(frames)} frames; expected {expected_frames}; total={n_total}')
        universe.trajectory[frames[0]]
        box = _box(universe.dimensions)
        chains = _chains(universe)
        # Uniform bins, at most 1 Angstrom wide, cover exactly one period.
        edges = np.linspace(0.0, box[2], int(np.ceil(box[2])) + 1)
        centers = (edges[:-1] + edges[1:]) / 2
        output.parent.mkdir(parents=True, exist_ok=True)
        records = []
        with TemporaryDirectory(prefix='contact_center_') as tmp:
            tmp = Path(tmp)
            work = tmp / 'centered.dcd'
            with mda.Writer(str(work), universe.atoms.n_atoms) as writer:
                for ordinal, frame in enumerate(frames):
                    ts = universe.trajectory[frame]
                    if not np.allclose(_box(ts.dimensions), box, rtol=0, atol=1e-3):
                        raise ValueError('Box varies with time; fixed-box centering is required')
                    xyz = universe.atoms.positions.astype(float)
                    if not np.all(np.isfinite(xyz)):
                        raise ValueError(f'Nonfinite coordinates in frame {frame}')
                    # Histogram periodic positions, so out-of-box beads cannot be omitted.
                    xyz = np.mod(xyz, box)
                    counts, _ = np.histogram(xyz[:, 2], bins=edges)
                    shift1 = float(box[2] / 2 - centers[np.argmax(counts)])
                    xyz[:, 2] += shift1
                    xyz = np.mod(xyz, box)
                    counts, _ = np.histogram(xyz[:, 2], bins=edges)
                    patch_z, patch_h = _largest_density_patch(centers, counts, box[2])
                    if not patch_h.size:
                        raise ValueError(f'No occupied density region in frame {frame}')
                    midpoint = float(np.average(patch_z, weights=patch_h))
                    shift2 = float(box[2] / 2 - midpoint)
                    xyz[:, 2] += shift2
                    xyz = np.mod(xyz, box)
                    whole, _ = _make_chains_whole(xyz, chains, box, max_bond_nm * 10)
                    universe.atoms.positions = whole
                    writer.write(universe.atoms)
                    records.append(dict(output_frame=ordinal, source_frame=frame,
                        shift_peak_nm=shift1 / 10, shift_patch_nm=shift2 / 10))
            report = tmp / 'diagnostics.csv'
            with report.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=records[0].keys())
                writer.writeheader()
                writer.writerows(records)
            shutil.copyfile(work, output)
            shutil.copyfile(report, output.with_suffix('.diagnostics.csv'))
        return dict(nframes_total=n_total, nframes_written=len(frames),
                    box_A=box.tolist(), n_chains=len(chains))
    finally:
        universe.trajectory.close()
