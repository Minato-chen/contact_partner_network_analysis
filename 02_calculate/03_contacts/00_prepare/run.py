"""Independent calculation entry; scientific defaults are defined below."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "project.py").is_file())
sys.path.insert(0, str(ROOT))
import project
import result_files as files

PARAMETERS = dict(
    contact_cutoff_nm=1.0, exclude_sequence_separation=2, chunk_frames=25, block_frames=50
)
OUTPUT = HERE / "output"
io = project.load_module("prepare_io", ROOT / "02_calculate/dependencies/analysis_io.py")
prepare = project.load_module("prepare_contacts", HERE / "dependencies/prepare.py")
network = project.load_module("intra_network", HERE.parent / "dependencies/network.py")


def calculate(simulation, centered_root=ROOT / "01_center/output", output=OUTPUT):
    stage = files.reset_output(output, simulation["id"])
    io.validate_parameters(PARAMETERS)
    inputs = prepare.compute(simulation, centered_root, stage, PARAMETERS, io, network)
    files.finish(stage, simulation, PARAMETERS, inputs)


if __name__ == "__main__":
    raise SystemExit(project.entry_main(HERE, calculate))
