#!/usr/bin/env python3
"""Stage 2 upstream: density profile and measured dense-core boundary."""
from pathlib import Path
import math
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import project
import result_files as files
from datetime import datetime, timezone

# Same half-height convention as ana_v4; no fixed-width fallback.
PARAMETERS = dict(
    bin_width_nm=0.1,
    peak_halfwidth_bins=5,
    dense_coefficient=0.6,
    min_peak_concentration_mM=1.0,
    chunk_frames=25,
)
OUTPUT = HERE / 'output'
CENTERED = project.ROOT / '01_center/output'
io = project.load_module('boundary_input_io', HERE.parent / 'dependencies/analysis_io.py')
profile = project.load_module('density_profile', HERE / 'dependencies/profile.py')


def validate():
    for key in ('bin_width_nm', 'min_peak_concentration_mM'):
        if not math.isfinite(PARAMETERS[key]) or PARAMETERS[key] <= 0:
            raise ValueError(f'{key} must be positive and finite')
    if not 0 < PARAMETERS['dense_coefficient'] <= 1:
        raise ValueError('dense_coefficient must be in (0, 1]')
    for key in ('peak_halfwidth_bins', 'chunk_frames'):
        value = PARAMETERS[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < (1 if key == 'chunk_frames' else 0):
            raise ValueError(f'Invalid integer {key}')


def calculate(simulation, centered_root=CENTERED, output=OUTPUT):
    validate()
    stage = files.reset_output(output, simulation['id'])
    upstream = io.centered_reference(simulation, centered_root)
    profile.compute(simulation, centered_root, stage, PARAMETERS, io)
    files.write_json(stage / 'metadata.json', dict(simulation_id=simulation['id'],
        parameters=PARAMETERS, upstream=upstream,
        created_at=datetime.now(timezone.utc).isoformat()))
    (stage / 'COMPLETE').touch()
    print(f'[ok] {stage}', flush=True)


def main():
    parser = project.parser(__doc__)
    parser.add_argument('--centered-dir', type=Path, default=CENTERED)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    validate()
    simulations = project.select_simulations(args.config, args.simulation)
    failed = False
    outcomes = []
    for simulation in simulations:
        try:
            calculate(simulation, args.centered_dir, args.output_dir)
            outcomes.append((simulation['id'], ['SUCCESS']))
        except Exception as exc:
            failed = True
            outcomes.append((simulation['id'], ['FAILED']))
            print(f'[failed] {simulation["id"]}: {exc}', file=sys.stderr)
    project.print_summary(outcomes, [HERE.name])
    return failed


if __name__ == '__main__':
    raise SystemExit(main())
