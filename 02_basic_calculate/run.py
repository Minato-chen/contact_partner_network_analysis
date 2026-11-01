"""Run all basic property calculations (phase boundary, Rg, csat, SCD, kappa)."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files
from common.runner import make_cli_main, run_batch


def entries():
    from common.runner import find_module_entries
    return find_module_entries(HERE)


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    tasks=None,
    objects=None,
):
    return run_batch(
        HERE,
        simulation,
        centered_root,
        results_root,
        tasks=tasks,
        objects=objects,
        entry_type="modules",
    )


main = make_cli_main(HERE, __doc__, entry_type="modules")

if __name__ == "__main__":
    raise SystemExit(main())
