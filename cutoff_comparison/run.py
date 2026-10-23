from importlib import import_module

"""Independent four-distance study. Does not modify formal analysis outputs."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import inputs
from common import results as files

PARAMETERS = dict(
    cutoffs_nm=[0.7, 0.8, 0.9, 1.0],
    exclude_sequence_separation=2,
    partner_thresholds=[1, 5, 10],
    block_frames=50,
    chunk_frames=25,
)
OUTPUT = HERE / "output"
io = import_module("02_calculate.common.analysis_io")
phase_entry = import_module("02_calculate.01_basic.phase_boundary")
study = import_module("cutoff_comparison.dependencies.compare")


def calculate(simulation, centered_root=ROOT / "01_center/output", output=OUTPUT):
    import numpy as np

    cutoffs = PARAMETERS["cutoffs_nm"]
    if (
        cutoffs != sorted(set(cutoffs))
        or not cutoffs
        or not np.isfinite(cutoffs).all()
        or min(cutoffs) <= 0
    ):
        raise ValueError("Distances must be distinct, sorted, finite and positive")
    io.validate_parameters(PARAMETERS)
    # Measure one boundary using the existing method, entirely under this study.
    phase_root = Path(output) / "upstream/phase_boundary"
    phase_entry.calculate(simulation, centered_root, output=phase_root)
    phase = files.completed(phase_root / simulation["id"])
    selection = files.reset_output(Path(output) / "upstream/selection", simulation["id"])
    targets = {
        r: files.reset_output(Path(output) / study.cutoff_name(r), simulation["id"])
        for r in cutoffs
    }
    study.compare(simulation, centered_root, phase, targets, selection, PARAMETERS, io)
    inputs = dict(
        centered=io.centered_reference(simulation, centered_root),
        phase=files.read_json(phase / "metadata.json"),
        selection="same chains for all distances",
    )
    files.finish(selection, simulation, dict(block_frames=PARAMETERS["block_frames"]), inputs)
    for cutoff, target in targets.items():
        files.finish(target, simulation, PARAMETERS | dict(contact_cutoff_nm=cutoff), inputs)


def main():
    cli = inputs.parser(__doc__)
    cli.add_argument("--centered-dir", type=Path, default=ROOT / "01_center/output")
    cli.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = cli.parse_args()
    rows = []
    for simulation in inputs.select_simulations(args.config, args.simulation):
        try:
            calculate(simulation, args.centered_dir, args.output_dir)
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
