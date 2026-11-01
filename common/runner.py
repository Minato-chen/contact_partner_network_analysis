"""Shared discovery and execution helpers for calculations."""
from importlib import import_module
from pathlib import Path
import sys

from common import inputs
from common import results as files


def find_module_entries(folder):
    """Find calculation scripts in a flat folder (e.g. 02_basic_calculate)."""
    folder = Path(folder).resolve()
    entries = [
        (p.stem, p)
        for p in folder.glob("*.py")
        if p.name != "run.py" and not p.name.startswith("_")
    ]
    return sorted(entries, key=lambda item: (item[0] != "phase_boundary", item[0]))


def find_dir_entries(folder):
    """Find calculation subdirectories containing a run.py (e.g. 02_calculate)."""
    folder = Path(folder).resolve()
    return [
        (p.name, p / "run.py")
        for p in sorted(folder.iterdir())
        if p.is_dir() and (p / "run.py").is_file()
    ]


def run_batch(folder, simulation, centered_root, results_root, tasks=None, objects=None, entry_type="modules"):
    folder = Path(folder).resolve()
    root = next(p for p in (folder, *folder.parents) if (p / "simulations.yaml").is_file())
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))

    if entry_type == "modules":
        children = find_module_entries(folder)
    else:
        children = find_dir_entries(folder)

    if tasks is not None:
        unknown = set(tasks) - {name for name, _ in children}
        if unknown:
            raise ValueError(f"Unknown calculations: {sorted(unknown)}")
        children = [(name, path) for name, path in children if name in tasks]
    if not children:
        raise ValueError(f"No calculations in {folder}")

    statuses = []
    for name, path in children:
        try:
            mod_parts = path.relative_to(root).with_suffix("").parts
            module = import_module(".".join(mod_parts))
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
        raise RuntimeError(f"One or more calculations failed in {folder}")
    if "SKIPPED" in statuses:
        raise files.InputUnavailable(f"One or more calculations skipped in {folder}")


def make_cli_main(folder, doc, entry_type="modules", single_task=None):
    folder = Path(folder).resolve()
    root = next(p for p in (folder, *folder.parents) if (p / "simulations.yaml").is_file())

    def main(task=single_task):
        cli = inputs.parser(doc)
        cli.add_argument("--centered-dir", type=Path, default=root / "01_center/output")
        cli.add_argument("--results-dir", type=Path, default=files.RESULTS)
        cli.add_argument(
            "--object", dest="objects", action="append", choices=["central_chain", "dense_chains"]
        )
        cli.add_argument(
            "--task", action="append", help="Task name; repeat to select several"
        )
        args = cli.parse_args()
        selected_tasks = [task] if task is not None else args.task
        rows = []
        for simulation in inputs.select_simulations(args.config, args.simulation):
            try:
                run_batch(
                    folder,
                    simulation,
                    args.centered_dir,
                    args.results_dir,
                    tasks=selected_tasks,
                    objects=args.objects,
                    entry_type=entry_type,
                )
                status = "SUCCESS"
            except files.InputUnavailable as exc:
                status = "SKIPPED"
                print(f'[skipped] {simulation["id"]}: {exc}', file=sys.stderr, flush=True)
            except Exception as exc:
                status = "FAILED"
                print(f'[failed] {simulation["id"]}: {exc}', file=sys.stderr, flush=True)
            rows.append((simulation["id"], [status]))
        files.print_summary(rows, [folder.name])
        return int(any(statuses[0] != "SUCCESS" for _, statuses in rows))

    return main
