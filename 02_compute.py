#!/usr/bin/env python3
"""Step 2: create chain-level CSVs and chain-pair NPZs from centred T260 DCDs."""
from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from dependencies.contact_persistence_test import load_edges
from dependencies.partner_analysis import analyze_trajectory, summarize_blocks

DATA_ROOT = Path('/Users/minato/Desktop/slab_idr')
HERE = Path(__file__).resolve().parent
CENTERED = HERE / 'centered'
RESULTS = HERE / 'results'
VARIANTS = {'A18': (11.25, 37.5), 'B01': (8.91, 29.7), 'B10': (9.81, 32.7)}
EXPECTED_FRAMES = 750
BLOCK_FRAMES = 50


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    frame.to_csv(temporary, index=False, compression='gzip' if path.suffix == '.gz' else None)
    temporary.replace(path)


def write_npz(contacts: np.ndarray, both_dense: np.ndarray, path: Path) -> None:
    if contacts.shape != both_dense.shape or len(contacts) != EXPECTED_FRAMES:
        raise ValueError(f'{path.name}: unexpected array shape')
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('wb') as handle:
        np.savez_compressed(handle, contact_counts=contacts, both_dense=both_dense)
    temporary.replace(path)


def partner_distributions(rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Frame-normalized fractions for thresholds 1 and 5, matching notebook 02."""
    selected = rows[rows.threshold.isin((1, 5))][
        ['variant', 'threshold', 'frame', 'block', 'chain', 'partner_count']]
    keys = ['variant', 'threshold', 'frame', 'block']
    counts = selected.groupby(keys + ['partner_count']).size().rename('n_chains').reset_index()
    populations = selected.groupby(keys).size().rename('n_dense').reset_index()
    counts = counts.merge(populations, on=keys, validate='many_to_one')
    counts['fraction'] = counts.n_chains / counts.n_dense
    index = pd.MultiIndex.from_product(
        [[f'{v}_T260' for v in VARIANTS], (1, 5), range(EXPECTED_FRAMES),
         range(int(selected.partner_count.max()) + 1)],
        names=['variant', 'threshold', 'frame', 'partner_count'])
    fractions = (counts.set_index(['variant', 'threshold', 'frame', 'partner_count'])
                 ['fraction'].reindex(index, fill_value=0).rename('fraction').reset_index())
    fractions = fractions.merge(selected[keys].drop_duplicates(),
                                on=['variant', 'threshold', 'frame'], validate='many_to_one')
    check = fractions.groupby(['variant', 'threshold', 'frame']).fraction.sum()
    if not np.allclose(check, 1):
        raise AssertionError('Partner-count fractions must sum to one per frame')
    mean = (fractions.groupby(['variant', 'threshold', 'partner_count'], as_index=False)
            .fraction.mean().rename(columns={'fraction': 'mean_frame_fraction'}))
    block = (fractions.groupby(['variant', 'threshold', 'block', 'partner_count'], as_index=False)
             .fraction.mean().rename(columns={'fraction': 'mean_block_fraction'}))
    by_frame = (selected.groupby(keys, as_index=False)
                .agg(n_dense=('chain', 'size'), mean_partners=('partner_count', 'mean')))
    diagnostic = (by_frame.groupby(['variant', 'threshold', 'block'], as_index=False)
                  .agg(mean_n_dense=('n_dense', 'mean'), mean_partners=('mean_partners', 'mean')))
    return mean, block, diagnostic


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true', help='Recompute existing outputs')
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for variant, (z_den, z_dil) in VARIANTS.items():
        name = f'tau_k18_129aa_{variant}_T260_5'
        pdb = DATA_ROOT / name / 'top.pdb'
        dcd = CENTERED / f'{variant}_T260_centered.dcd'
        for path in (pdb, dcd):
            if not path.is_file():
                raise FileNotFoundError(path)
        rows_path = RESULTS / f'{variant}_T260_chain_partner_rows.csv.gz'
        edges_path = RESULTS / f'{variant}_edge_arrays.npz'
        with TemporaryDirectory(prefix='tau_k18_read_') as temp_dir:
            # mdtraj's DCD reader requires an ASCII pathname on this machine.
            alias = Path(temp_dir) / f'{variant}_T260_centered.dcd'
            alias.symlink_to(dcd)
            if rows_path.exists() and not args.force:
                print(f'[reuse] {rows_path.name}', flush=True)
                rows = pd.read_csv(rows_path)
            else:
                print(f'[compute] {variant}: chain contacts', flush=True)
                rows = analyze_trajectory(f'{variant}_T260', pdb, alias, z_den, z_dil,
                                          start=0, step=1, block_frames=BLOCK_FRAMES)
                if rows.frame.nunique() != EXPECTED_FRAMES:
                    raise ValueError(f'{variant}: expected {EXPECTED_FRAMES} frames')
                write_csv(rows, rows_path)
            if rows.frame.nunique() != EXPECTED_FRAMES:
                raise ValueError(f'{variant}: saved rows have unexpected frame count')
            all_rows.append(rows)
            if edges_path.exists() and not args.force:
                print(f'[reuse] {edges_path.name}', flush=True)
            else:
                print(f'[compute] {variant}: chain-pair arrays', flush=True)
                contacts, both_dense = load_edges(pdb, alias, z_den, EXPECTED_FRAMES)
                write_npz(contacts, both_dense, edges_path)
    summary = summarize_blocks(pd.concat(all_rows, ignore_index=True))
    write_csv(summary, RESULTS / 'chain_partner_summary.csv')
    mean, blocks, diagnostic = partner_distributions(pd.concat(all_rows, ignore_index=True))
    write_csv(mean, RESULTS / 'partner_count_distribution.csv')
    write_csv(blocks, RESULTS / 'partner_count_distribution_by_block.csv')
    write_csv(diagnostic, RESULTS / 'partner_count_time_blocks.csv')
    print(f'[ok] {EXPECTED_FRAMES} frames per variant; CSV and NPZ outputs: {RESULTS}')


if __name__ == '__main__':
    main()
