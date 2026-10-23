#!/usr/bin/env python3
"""Stage 2: compute enabled contact analyses for all configured simulations."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import project

# Enable only explicitly supported analyses. Scientific parameters live in their run.py.
ANALYSES = ('phase_boundary', 'intramolecular', 'intermolecular')
# None selects all registered simulations. CLI --simulation overrides this selection.
SIMULATIONS = None


def main():
    parser = project.parser(__doc__)
    parser.add_argument('--analysis', choices=ANALYSES, action='append')
    parser.add_argument('--centered-dir', type=Path, default=HERE.parent / '01_center/output')
    args = parser.parse_args()
    simulations = project.select_simulations(args.config, args.simulation or SIMULATIONS)
    selected = list(args.analysis or ANALYSES)
    # Boundary calculation is a required upstream task, even for a single contact analysis.
    selected = ['phase_boundary'] + [name for name in selected if name != 'phase_boundary']
    modules = [project.load_module(name + '_entry', HERE / name / 'run.py') for name in selected]
    for module in modules:
        module.validate()
    failed = False
    outcomes = []
    for simulation in simulations:
        statuses = []
        boundary_failed = False
        for name, module in zip(selected, modules):
            if boundary_failed:
                statuses.append('SKIPPED')
                continue
            try:
                module.calculate(simulation, args.centered_dir)
                statuses.append('SUCCESS')
            except Exception as exc:
                failed = True
                statuses.append('FAILED')
                print(f'[failed] {simulation["id"]}/{name}: {exc}', file=sys.stderr, flush=True)
                # The two contact analyses are independent once the boundary succeeds.
                boundary_failed = name == 'phase_boundary'
        outcomes.append((simulation['id'], statuses))
    project.print_summary(outcomes, selected)
    return failed


if __name__ == '__main__':
    raise SystemExit(main())
