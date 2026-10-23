"""Run intra analysis for selected simulations and objects."""
from importlib import import_module
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import inputs
from common import results as files

analysis = import_module("02_calculate.02_contacts.01_intra.intra")
calculate = analysis.calculate


def main():
    cli = inputs.parser(__doc__)
    cli.add_argument("--centered-dir", type=Path, default=ROOT / "01_center/output")
    cli.add_argument("--results-dir", type=Path, default=files.RESULTS)
    cli.add_argument(
        "--object", dest="objects", action="append", choices=["central_chain", "dense_chains"]
    )
    args = cli.parse_args()
    rows = []
    for simulation in inputs.select_simulations(args.config, args.simulation):
        try:
            calculate(simulation, args.centered_dir, args.results_dir, objects=args.objects)
            status = "SUCCESS"
        except files.InputUnavailable as exc:
            status = "SKIPPED"
            print(f'[skipped] {simulation["id"]}: {exc}', file=sys.stderr, flush=True)
        except Exception as exc:
            status = "FAILED"
            print(f'[failed] {simulation["id"]}: {exc}', file=sys.stderr, flush=True)
        rows.append((simulation["id"], [status]))
    files.print_summary(rows, [HERE.name])
    return int(any(statuses[0] != "SUCCESS" for _, statuses in rows))


if __name__ == "__main__":
    raise SystemExit(main())
