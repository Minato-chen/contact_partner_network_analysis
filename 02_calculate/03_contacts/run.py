"""Run all direct child tasks in numbered order."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "project.py").is_file())
sys.path.insert(0, str(ROOT))
import project
import result_files as files

SIMULATIONS = None


def calculate(simulation, centered_root=ROOT / "01_center/output", tasks=None):
    project.run_children(HERE, simulation, centered_root, selected=tasks)


if __name__ == "__main__":
    raise SystemExit(project.entry_main(HERE, calculate, parent=True, simulations=SIMULATIONS))
