"""Run calculations in this folder and report each outcome."""
from importlib import import_module
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import inputs
from common import results as files


def entries():
    return [
        (p.name, p / "run.py")
        for p in sorted(HERE.iterdir())
        if p.is_dir() and (p / "run.py").is_file()
    ]


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    tasks=None,
    objects=None,
):
    children = entries()
    if tasks is not None:
        unknown = set(tasks) - {name for name, _ in children}
        if unknown:
            raise ValueError(f"Unknown calculations: {sorted(unknown)}")
        children = [(name, path) for name, path in children if name in tasks]
    if not children:
        raise ValueError(f"No calculations in {HERE}")
    statuses = []
    for name, path in children:
        try:
            module = import_module(".".join(path.relative_to(ROOT).with_suffix("").parts))
            module.calculate(simulation, centered_root, results_root, objects=objects)
            statuses.append("SUCCESS")
        except files.InputUnavailable as exc:
            statuses.append("SKIPPED")
            print(f'[skipped] {simulation["id"]}/{name}: {exc}', file=sys.stderr, flush=True)
        except Exception as exc:
            statuses.append("FAILED")
            print(f'[failed] {simulation["id"]}/{name}: {exc}', file=sys.stderr, flush=True)
    files.print_summary([(simulation["id"], statuses)], [name for name, _ in children])
    if "FAILED" in statuses:
        raise RuntimeError(f"One or more calculations failed in {HERE}")
    if "SKIPPED" in statuses:
        raise files.InputUnavailable(f"One or more calculations skipped in {HERE}")


def main():
    cli = inputs.parser(__doc__)
    cli.add_argument("--centered-dir", type=Path, default=ROOT / "01_center/output")
    cli.add_argument("--results-dir", type=Path, default=files.RESULTS)
    cli.add_argument(
        "--object", dest="objects", action="append", choices=["central_chain", "dense_chains"]
    )
    cli.add_argument("--task", action="append", choices=[name for name, _ in entries()])
    args = cli.parse_args()
    rows = []
    for simulation in inputs.select_simulations(args.config, args.simulation):
        try:
            calculate(
                simulation,
                args.centered_dir,
                args.results_dir,
                objects=args.objects,
                tasks=args.task,
            )
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
