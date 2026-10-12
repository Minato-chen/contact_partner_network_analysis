#!/usr/bin/env python3
"""Standalone slab centering for local tests; does not import ana_v4."""

from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np


def _largest_density_patch(z: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    best_z: list[float] = []
    best_h: list[float] = []
    current_z: list[float] = []
    current_h: list[float] = []
    for zi, value in zip(z, counts):
        if value > 0:
            current_z.append(float(zi))
            current_h.append(float(value))
        else:
            if sum(current_h) > sum(best_h):
                best_z, best_h = current_z, current_h
            current_z, current_h = [], []
    if sum(current_h) > sum(best_h):
        best_z, best_h = current_z, current_h
    return np.asarray(best_z), np.asarray(best_h)


def center(pdb: Path, dcd: Path, output: Path, start: int, step: int) -> None:
    if output.exists():
        print(f"[skip] centered trajectory exists: {output}")
        return
    universe = mda.Universe(str(pdb), str(dcd), in_memory=True)
    atoms = universe.atoms
    lz = float(universe.dimensions[2])
    edges = np.arange(0.0, lz + 1.0, 1.0)
    centers = edges[:-1] + 0.5

    bonds = []
    for segment in universe.segments:
        indices = segment.atoms.indices
        bonds.extend((int(a), int(b)) for a, b in zip(indices[:-1], indices[1:]))
    universe.add_TopologyAttr("bonds", bonds)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.stem + ".tmp.dcd")
    with mda.Writer(str(temporary), atoms.n_atoms) as writer:
        for ts in universe.trajectory[start::step]:
            counts, _ = np.histogram(atoms.positions[:, 2], bins=edges)
            z_peak = centers[np.argmax(counts)]
            atoms.translate([0.0, 0.0, -z_peak + 0.5 * lz])
            transformations.wrap(atoms)(ts)

            counts, _ = np.histogram(atoms.positions[:, 2], bins=edges)
            patch_z, patch_h = _largest_density_patch(centers, counts)
            if patch_h.size and patch_h.sum() > 0:
                z_mid = np.average(patch_z, weights=patch_h)
                atoms.translate([0.0, 0.0, -z_mid + 0.5 * lz])
                transformations.wrap(atoms)(ts)
            transformations.unwrap(atoms)(ts)
            writer.write(atoms)
    temporary.replace(output)
    print(f"[ok] centered {len(universe.trajectory[start::step])} frames -> {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--dcd", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, default=750)
    parser.add_argument("--step", type=int, default=3)
    args = parser.parse_args()
    center(args.pdb, args.dcd, args.output, args.start, args.step)


if __name__ == "__main__":
    main()
