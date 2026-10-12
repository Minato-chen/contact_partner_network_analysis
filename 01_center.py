#!/usr/bin/env python3
"""Step 1: sample and centre the three original T260 slab trajectories."""
from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from dependencies.center_local import center

DATA_ROOT = Path('/Users/minato/Desktop/slab_idr')
HERE = Path(__file__).resolve().parent
VARIANTS = ('A18', 'B01', 'B10')
START = 750
STEP = 3
EXPECTED_FRAMES = 750


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=VARIANTS, action='append',
                        help='Run selected variants; default: all three')
    parser.add_argument('--force', action='store_true', help='Rebuild existing centred DCDs')
    args = parser.parse_args()
    output_dir = HERE / 'centered'
    output_dir.mkdir(parents=True, exist_ok=True)
    for variant in args.variant or VARIANTS:
        name = f'tau_k18_129aa_{variant}_T260_5'
        pdb = DATA_ROOT / name / 'top.pdb'
        raw = DATA_ROOT / name / f'{name}.dcd'
        output = output_dir / f'{variant}_T260_centered.dcd'
        for path in (pdb, raw):
            if not path.is_file():
                raise FileNotFoundError(path)
        if output.exists() and not args.force:
            print(f'[skip] {variant}: {output}')
            continue
        # MDAnalysis/DCD libraries can reject non-ASCII output paths. Write to an
        # ASCII temporary path, then move the completed DCD into this folder.
        with TemporaryDirectory(prefix='tau_k18_center_') as temp_dir:
            temporary = Path(temp_dir) / f'{variant}_T260_centered.dcd'
            center(pdb, raw, temporary, START, STEP)
            temporary.replace(output)
        print(f'[ok] {variant}: original frames {START},{START+STEP},...; output {output}')
    print('Step 1 complete. Step 2: python 02_compute.py')


if __name__ == '__main__':
    main()
