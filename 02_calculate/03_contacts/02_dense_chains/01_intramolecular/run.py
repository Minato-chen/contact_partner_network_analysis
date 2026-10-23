"""Independent calculation entry; scientific defaults are defined below."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "project.py").is_file())
sys.path.insert(0, str(ROOT))
import project
import result_files as files

OBJECT = "dense_chains"
PARAMETERS = dict(block_frames=50)
OUTPUT = HERE / "output"
io = project.load_module("dense_chains_intra_io", ROOT / "02_calculate/dependencies/analysis_io.py")
shared = HERE.parents[1] / "dependencies"
contact_io = project.load_module("dense_chains_intra_files", shared / "contact_io.py")
analysis = project.load_module("dense_chains_intra_analysis", shared / "intra.py")


def calculate(simulation, centered_root=ROOT / "01_center/output", output=OUTPUT):
    stage = files.reset_output(output, simulation["id"])
    if PARAMETERS["block_frames"] <= 0:
        raise ValueError("block_frames must be positive")
    inputs = analysis.compute(simulation, centered_root, stage, PARAMETERS, OBJECT, io, contact_io)
    files.finish(stage, simulation, PARAMETERS, inputs)


if __name__ == "__main__":
    raise SystemExit(project.entry_main(HERE, calculate))
