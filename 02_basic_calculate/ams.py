"""Multiscale charge matrix spectral asymmetry (A_MS)."""
from importlib import import_module
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

read_sequence = import_module("02_basic_calculate._sequence").read_sequence
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
SCALES = (5.0, 10.0, 20.0)


def _lambda_fraction(mask, length_scale):
    positions = np.arange(len(mask), dtype=float)
    kernel = np.exp(-np.abs(positions[:, None] - positions[None, :]) / float(length_scale))
    matrix = mask[:, None] * kernel * mask[None, :]
    values = np.clip(np.linalg.eigvalsh(matrix), 0.0, None)
    trace = float(values.sum())
    if trace <= 0.0:
        raise ValueError("A_MS requires both positive and negative residues in sequence")
    return float(values[-1] / trace)


def compute_ams(sequence, scales=SCALES):
    if not sequence:
        raise ValueError("A_MS requires a nonempty sequence")
    q = np.array([CHARGE.get(aa, 0) for aa in sequence], dtype=int)
    pos = (q > 0).astype(float)
    neg = (q < 0).astype(float)
    if pos.sum() == 0 or neg.sum() == 0:
        raise ValueError("A_MS requires both positive and negative residues in sequence")
    scales_arr = np.asarray(scales, dtype=float)
    if scales_arr.ndim != 1 or len(scales_arr) < 2 or np.any(scales_arr <= 0):
        raise ValueError("scales must be a 1D sequence of at least two positive numbers")
    deltas = np.array(
        [_lambda_fraction(pos, ell) - _lambda_fraction(neg, ell) for ell in scales_arr],
        dtype=float,
    )
    log_scales = np.log(scales_arr)
    trapz_fn = getattr(np, "trapezoid", np.trapz)
    ams = float(trapz_fn(deltas, x=log_scales) / (log_scales[-1] - log_scales[0]))
    if not math.isfinite(ams):
        raise ValueError("Computed A_MS is not finite")
    return ams


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    output = Path(results_root) / simulation["id"] / "01_basic/ams.scalars.json"
    files.start_scalar(output, simulation)
    source = files.completed(Path(centered_root) / simulation["id"])
    sequence, n_chains = read_sequence(source)
    definition = dict(
        formula="(delta(5)+2*delta(10)+delta(20))/4 where delta(ell)=f_pos(ell)-f_neg(ell)",
        kernel="exp(-|i-j|/ell)",
        scales=list(SCALES),
        charge=CHARGE,
        other_residues=0,
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
    value = compute_ams(sequence)
    files.finish_scalar(
        output,
        [files.scalar_record("ams", value, "dimensionless", definition=definition)],
    )


if __name__ == "__main__":
    raise SystemExit(import_module("02_basic_calculate.run").main(task="ams"))
