#!/usr/bin/env python3
"""Stage 2: intermolecular contacts. Edit scientific defaults below."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import project
import result_files as files
from datetime import datetime, timezone

PARAMETERS = dict(
    contact_cutoff_nm=0.7,  # Strict d < cutoff; partner counts use contacts >= threshold.
    partner_thresholds=[1, 5, 10],
    network_thresholds=[1, 5, 10, 20, 30, 35],
    block_frames=50,
    chunk_frames=25,
)
OUTPUT = HERE / 'output'
CENTERED = project.ROOT / '01_center/output'
io = project.load_module('contact_analysis_io', HERE.parent / 'dependencies/analysis_io.py')
contacts = project.load_module('intermolecular_contacts', HERE / 'dependencies/contacts.py')


def validate():
    io.validate_parameters(PARAMETERS)


def calculate(simulation, centered_root=CENTERED, output=OUTPUT):
    validate()
    stage = files.reset_output(output, simulation['id'])
    upstream = files.read_json(io.phase_source(simulation, centered_root) / 'metadata.json')
    contacts.compute(simulation, centered_root, stage, PARAMETERS, io)
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
