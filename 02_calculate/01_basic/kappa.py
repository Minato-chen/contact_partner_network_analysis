"""Sequence charge patterning (kappa), calculated by localCIDER."""
from importlib import import_module
from importlib.metadata import version
from pathlib import Path
import math
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

read_sequence = import_module("02_calculate.01_basic._sequence").read_sequence


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    output = Path(results_root) / simulation["id"] / "01_basic/kappa.scalars.json"
    files.start_scalar(output, simulation)
    source = files.completed(Path(centered_root) / simulation["id"])
    sequence, n_chains = read_sequence(source)
    from localcider.sequenceParameters import SequenceParameters

    definition = dict(
        method="localcider.SequenceParameters.get_kappa", version=version("localcider")
    )
    files.start_scalar(
        output,
        simulation,
        parameters=definition,
        inputs=dict(
            sequence=sequence,
            n_chains=n_chains,
            centered=files.read_json(source / "metadata.json"),
        ),
    )
    value = SequenceParameters(sequence).get_kappa()
    if not math.isfinite(value) or value < 0:
        raise ValueError("localCIDER returned undefined kappa for this sequence")
    files.finish_scalar(
        output,
        [files.scalar_record("kappa", value, "dimensionless", definition=definition)],
    )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.01_basic.run").main(task="kappa"))
