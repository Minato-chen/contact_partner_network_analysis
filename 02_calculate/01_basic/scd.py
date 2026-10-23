"""Sequence charge decoration using the supplied pairwise formula."""
from importlib import import_module
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

read_sequence = import_module("02_calculate.01_basic._sequence").read_sequence
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}


def compute_scd(sequence):
    if not sequence:
        raise ValueError("SCD requires a nonempty sequence")
    q = [CHARGE.get(aa, 0) for aa in sequence]
    N = len(q)
    total = 0
    for m in range(1, N):
        for n in range(m):
            total += q[m] * q[n] * (m - n) ** 0.5
    return total / N


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    output = Path(results_root) / simulation["id"] / "01_basic/scd.scalars.json"
    files.start_scalar(output, simulation)
    source = files.completed(Path(centered_root) / simulation["id"])
    sequence, n_chains = read_sequence(source)
    definition = dict(
        formula="sum(n<m, q[m]*q[n]*sqrt(m-n))/N", charge=CHARGE, other_residues=0
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
    value = compute_scd(sequence)
    files.finish_scalar(
        output,
        [files.scalar_record("scd", value, "dimensionless", definition=definition)],
    )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.01_basic.run").main(task="scd"))
