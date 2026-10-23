#!/usr/bin/env python3
"""Stage 1: batch centering of simulations registered in simulations.yaml."""
from importlib import import_module
from pathlib import Path
import sys
import math
import csv
import shutil

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import inputs
from common import results as files
from datetime import datetime, timezone

# Sampling and centering parameters belong here, never in the raw simulation registry.
PARAMETERS = dict(start=750, stop=None, step=3, expected_frames=750, max_bond_nm=1.0)
# None selects all registered simulations. CLI --simulation overrides this selection.
SIMULATIONS = None
# Optional per-simulation overrides, for example {'A18_T260': {'start': 0}}.
SIMULATION_PARAMETERS = {}
OUTPUT = HERE / 'output'
center_module = import_module('01_center.dependencies.center_local')


def settings(simulation):
    overrides = SIMULATION_PARAMETERS.get(simulation['id'], {})
    if set(overrides) - set(PARAMETERS):
        raise ValueError(f'Unknown centering parameters: {set(overrides) - set(PARAMETERS)}')
    values = PARAMETERS | overrides
    for key in ('start', 'step'):
        value = values[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < (0 if key == 'start' else 1):
            raise ValueError(f'{key} must be an integer with start >= 0 and step > 0')
    stop = values['stop']
    if stop is not None and (isinstance(stop, bool) or not isinstance(stop, int) or stop <= values['start']):
        raise ValueError('stop must be an integer greater than start, or null')
    count = values['expected_frames']
    if count is not None and (isinstance(count, bool) or not isinstance(count, int) or count <= 0):
        raise ValueError('expected_frames must be positive or null')
    if not math.isfinite(values['max_bond_nm']) or values['max_bond_nm'] <= 0:
        raise ValueError('max_bond_nm must be finite and positive')
    return values


def calculate(simulation, output=OUTPUT):
    parameters = settings(simulation)
    time_step = simulation.get('time_step_ps')
    origin = simulation.get('time_origin_ps', 0.0)
    if time_step is not None and (not math.isfinite(time_step) or time_step <= 0):
        raise ValueError('time_step_ps must be positive and finite')
    if not math.isfinite(origin):
        raise ValueError('time_origin_ps must be finite')
    stage = files.reset_output(output, simulation['id'])
    sampling = center_module.center(Path(simulation['topology']), Path(simulation['trajectory']),
                                   stage / 'centered.dcd', **parameters)
    shutil.copyfile(simulation['topology'], stage / 'topology.pdb')
    with (stage / 'centered.diagnostics.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    with (stage / 'frame_index.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['frame', 'source_frame', 'time_ps'])
        writer.writeheader()
        for row in rows:
            source_frame = int(row['source_frame'])
            writer.writerow(dict(frame=int(row['output_frame']), source_frame=source_frame,
                time_ps=origin + source_frame * time_step if time_step is not None else ''))
    import MDAnalysis as mda
    universe = mda.Universe(str(stage / 'topology.pdb'))
    try:
        chains = center_module._chains(universe)
        files.write_json(stage / 'chains.json', [ids.tolist() for ids in chains])
    finally:
        universe.trajectory.close()
    files.write_json(stage / 'metadata.json', dict(simulation_id=simulation['id'],
        inputs={key: simulation[key] for key in ('topology', 'trajectory')},
        parameters=parameters, time_step_ps=time_step, time_origin_ps=origin,
        coordinate_unit='angstrom', created_at=datetime.now(timezone.utc).isoformat(), **sampling))
    (stage / 'COMPLETE').touch()
    print(f'[ok] {stage}', flush=True)


def main():
    parser = inputs.parser(__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    simulations = inputs.select_simulations(args.config, args.simulation or SIMULATIONS)
    for simulation in simulations:
        settings(simulation)
    failed = False
    outcomes = []
    for simulation in simulations:
        try:
            calculate(simulation, args.output_dir)
            outcomes.append((simulation['id'], ['SUCCESS']))
        except Exception as exc:
            failed = True
            outcomes.append((simulation['id'], ['FAILED']))
            print(f'[failed] {simulation["id"]}: {exc}', file=sys.stderr)
    files.print_summary(outcomes, [HERE.name])
    return failed


if __name__ == '__main__':
    raise SystemExit(main())
