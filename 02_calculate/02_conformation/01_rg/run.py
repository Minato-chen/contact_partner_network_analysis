"""Independent calculation entry; scientific defaults are defined below."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "project.py").is_file())
sys.path.insert(0, str(ROOT))
import project
import result_files as files

PARAMETERS = dict(chunk_frames=25, n_terminal_mass_Da=2.0, c_terminal_mass_Da=16.0)
OUTPUT = HERE / "output"
io = project.load_module("rg_input", ROOT / "02_calculate/dependencies/analysis_io.py")
rg = project.load_module("rg_calculation", HERE / "dependencies/rg.py")
MASS_FILE = HERE / "dependencies/residue_masses.csv"


def calculate(simulation, centered_root=ROOT / "01_center/output", output=OUTPUT):
    stage = files.reset_output(output, simulation["id"])
    if not isinstance(PARAMETERS["chunk_frames"], int) or PARAMETERS["chunk_frames"] <= 0:
        raise ValueError("chunk_frames must be positive")
    inputs = rg.compute(simulation, centered_root, stage, PARAMETERS, io, MASS_FILE)
    files.finish(stage, simulation, PARAMETERS, inputs)


if __name__ == "__main__":
    raise SystemExit(project.entry_main(HERE, calculate))
