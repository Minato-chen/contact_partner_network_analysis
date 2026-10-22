#!/usr/bin/env python3
"""Step 1: sample, centre and verify the three original T260 slab trajectories."""
from __future__ import annotations
import argparse
from pathlib import Path
from dependencies.center_local import center

DATA_ROOT = Path('/Users/minato/Desktop/slab_idr')
HERE = Path(__file__).resolve().parent
VARIANTS = ('A18', 'B01', 'B10')
START = 750
STEP = 3
EXPECTED_FRAMES = 750


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=VARIANTS, action='append')
    parser.add_argument('--force', action='store_true', help='Explicitly replace existing outputs')
    parser.add_argument('--output-dir', type=Path, default=HERE / 'centered',
                        help='Use a separate directory to preserve existing trajectories')
    args = parser.parse_args()
    for variant in args.variant or VARIANTS:
        name = f'tau_k18_129aa_{variant}_T260_5'
        pdb = DATA_ROOT / name / 'top.pdb'
        raw = DATA_ROOT / name / f'{name}.dcd'
        for path in (pdb, raw):
            if not path.is_file():
                raise FileNotFoundError(path)
        center(pdb, raw, args.output_dir / f'{variant}_T260_centered.dcd', START, STEP,
               expected_frames=EXPECTED_FRAMES, force=args.force)
    print('Step 1 complete. Rebuilt trajectories require matching downstream results.')


if __name__ == '__main__':
    main()
