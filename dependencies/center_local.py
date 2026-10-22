#!/usr/bin/env python3
"""Two-pass slab alignment with verified whole linear-chain coordinates.

MDAnalysis coordinates are in Angstrom. Only fixed orthorhombic boxes and
one bead per residue, ordered linear chains (one per segment) are supported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import MDAnalysis as mda
import numpy as np

ALGORITHM_VERSION = 'slab-two-pass-whole-linear-v1'


def _largest_density_patch(z: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Keep the existing linear nonempty-bin rule, including the final run."""
    runs = np.split(np.flatnonzero(counts > 0),
                    np.flatnonzero(np.diff(np.flatnonzero(counts > 0)) > 1) + 1)
    if not runs or not runs[0].size:
        return np.array([]), np.array([])
    best = max(runs, key=lambda ids: counts[ids].sum())
    return z[best], counts[best]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _input_identity(path: Path, *, content_hash: bool = False) -> dict:
    path = path.resolve()
    info = path.stat()
    result = dict(path=str(path), size=info.st_size, mtime_ns=info.st_mtime_ns)
    if content_hash:
        result['sha256'] = _sha256(path)
    return result


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


def _request(pdb: Path, dcd: Path, start: int, step: int,
             expected_frames: int | None, max_bond_nm: float) -> dict:
    return dict(algorithm=ALGORITHM_VERSION, implementation_sha256=_sha256(Path(__file__)),
                topology=_input_identity(pdb, content_hash=True),
                trajectory=_input_identity(dcd), start=start, step=step,
                expected_frames=expected_frames, max_bond_nm=max_bond_nm)


def _validate_output(pdb: Path, dcd: Path, chains: list[np.ndarray],
                     box: np.ndarray, n_expected: int, max_bond_A: float) -> None:
    """Reread every output frame, without loading the full DCD into memory."""
    check = mda.Universe(str(pdb), str(dcd))
    try:
        if len(check.trajectory) != n_expected:
            raise ValueError('Written trajectory has an unexpected number of frames')
        for ts in check.trajectory:
            if not np.allclose(_box(ts.dimensions), box, rtol=0, atol=1e-3):
                raise ValueError('Written trajectory box changed')
            for ids in chains:
                xyz = check.atoms.positions[ids].astype(float)
                if not np.all(np.isfinite(xyz)):
                    raise ValueError('Written coordinates are not finite')
                if np.linalg.norm(np.diff(xyz, axis=0), axis=1).max() > max_bond_A:
                    raise ValueError('Written chain is not whole')
                center = xyz.mean(axis=0)
                if np.any(center < -1e-3) or np.any(center >= box + 1e-3):
                    raise ValueError('Written chain center is outside the primary box')
    finally:
        check.trajectory.close()


def center(pdb: Path, dcd: Path, output: Path, start: int, step: int, *,
           expected_frames: int | None = None, force: bool = False,
           max_bond_nm: float = 1.0) -> None:
    """Write a verified DCD, frame diagnostics and a reusable provenance record.

    max_bond_nm is a generous topology/continuity check, NOT a contact cutoff.
    Existing files without matching provenance are refused unless force=True.
    """
    pdb, dcd, output = map(Path, (pdb, dcd, output))
    if output.resolve() in (pdb.resolve(), dcd.resolve()):
        raise ValueError('Output must not overwrite an input')
    if start < 0 or step <= 0 or not np.isfinite(max_bond_nm) or max_bond_nm <= 0:
        raise ValueError('Require start >= 0, step > 0 and finite max_bond_nm > 0')
    if expected_frames is not None and expected_frames <= 0:
        raise ValueError('expected_frames must be positive')
    before = _request(pdb, dcd, start, step, expected_frames, max_bond_nm)
    metadata = output.with_suffix('.metadata.json')
    diagnostics = output.with_suffix('.diagnostics.csv')
    if output.exists() and not force:
        try:
            saved = json.loads(metadata.read_text())
            valid = (saved['request'] == before
                     and saved['output_sha256'] == _sha256(output)
                     and saved['diagnostics_sha256'] == _sha256(diagnostics))
        except (OSError, ValueError, KeyError):
            valid = False
        if not valid:
            raise ValueError(f'{output}: missing or mismatched provenance; use a new output '
                             'or --force to explicitly rebuild')
        print(f'[skip] verified trajectory and settings: {output}')
        return
    universe = mda.Universe(str(pdb), str(dcd))
    try:
        n_total = len(universe.trajectory)
        frames = range(start, n_total, step)
        if not frames or (expected_frames is not None and len(frames) != expected_frames):
            raise ValueError(f'Selected {len(frames)} frames; expected {expected_frames}; total={n_total}')
        box = _box(universe.dimensions)
        chains = _chains(universe)
        # Preserve the existing 1-Angstrom histogram on the current fixed box.
        edges = np.arange(0.0, box[2] + 1.0, 1.0)
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
                    patch_z, patch_h = _largest_density_patch(centers, counts)
                    if not patch_h.size:
                        raise ValueError(f'No occupied density region in frame {frame}')
                    midpoint = float(np.average(patch_z, weights=patch_h))
                    shift2 = float(box[2] / 2 - midpoint)
                    # Diagnose ambiguity without silently changing the existing patch rule.
                    occupied = np.flatnonzero(counts > 0)
                    runs = np.split(occupied, np.flatnonzero(np.diff(occupied) > 1) + 1)
                    masses = sorted((int(counts[ids].sum()) for ids in runs), reverse=True)
                    xyz[:, 2] += shift2
                    xyz = np.mod(xyz, box)
                    whole, max_bond = _make_chains_whole(xyz, chains, box, max_bond_nm * 10)
                    universe.atoms.positions = whole
                    writer.write(universe.atoms)
                    final_counts, _ = np.histogram(np.mod(whole[:, 2], box[2]), bins=edges)
                    final_z, final_h = _largest_density_patch(centers, final_counts)
                    final_offset = float(np.average(final_z, weights=final_h) - box[2] / 2)
                    records.append(dict(output_frame=ordinal, source_frame=frame,
                        shift_peak_nm=shift1 / 10, shift_patch_nm=shift2 / 10,
                        patch_left_before_second_shift_nm=float(patch_z[0] - .5) / 10,
                        patch_right_before_second_shift_nm=float(patch_z[-1] + .5) / 10,
                        patch_width_nm=float(patch_z[-1] - patch_z[0] + 1) / 10,
                        patch_fraction=float(patch_h.sum() / len(xyz)),
                        second_patch_fraction=(masses[1] / len(xyz) if len(masses) > 1 else 0),
                        patch_touches_box_edge=bool(patch_z[0] == centers[0] or patch_z[-1] == centers[-1]),
                        final_patch_center_offset_nm=final_offset / 10,
                        max_backbone_bond_nm=max_bond / 10))
            _validate_output(pdb, work, chains, box, len(frames), max_bond_nm * 10)
            if _request(pdb, dcd, start, step, expected_frames, max_bond_nm) != before:
                raise ValueError('Inputs changed while centering')
            report = tmp / 'diagnostics.csv'
            with report.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=records[0].keys())
                writer.writeheader()
                writer.writerows(records)
            manifest = dict(request=before, nframes_total=n_total, nframes_written=len(frames),
                            box_A=box.tolist(), n_chains=len(chains),
                            output_sha256=_sha256(work), diagnostics_sha256=_sha256(report))
            # Metadata is the completion marker: publish it last. Failed/interrupted
            # writes cannot subsequently pass reuse verification.
            for source, target in [(work, output), (report, diagnostics)]:
                stage = target.with_name(target.name + '.pending')
                try:
                    with source.open('rb') as src, stage.open('wb') as dest:
                        shutil.copyfileobj(src, dest)
                    os.replace(stage, target)
                finally:
                    stage.unlink(missing_ok=True)
            stage = metadata.with_name(metadata.name + '.pending')
            try:
                stage.write_text(json.dumps(manifest, indent=2) + '\n')
                os.replace(stage, metadata)
            finally:
                stage.unlink(missing_ok=True)
        print(f'[ok] centered and verified {len(frames)} frames -> {output}')
    finally:
        universe.trajectory.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdb', type=Path, required=True)
    parser.add_argument('--dcd', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--start', type=int, default=750)
    parser.add_argument('--step', type=int, default=3)
    parser.add_argument('--expected-frames', type=int)
    parser.add_argument('--max-bond-nm', type=float, default=1.0)
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    center(args.pdb, args.dcd, args.output, args.start, args.step,
           expected_frames=args.expected_frames, force=args.force, max_bond_nm=args.max_bond_nm)


if __name__ == '__main__':
    main()
