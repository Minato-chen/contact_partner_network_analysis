"""Sequence charge spectrum, SCD sign-decomposition, and residue principal modes."""
from importlib import import_module
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

io = import_module("02_calculate.common.analysis_io")
read_sequence = import_module("02_basic_calculate._sequence").read_sequence

CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
PARAMETERS = dict(
    length_scales=[2.0, 5.0, 10.0, 20.0],
    integration_scales=[5.0, 10.0, 20.0],
    reference_scale=10.0,
)


def compute_scd_decomposition(q):
    n_res = len(q)
    idx = np.arange(n_res)
    m, n = np.meshgrid(idx, idx, indexing="ij")
    mask_upper = m > n
    weights = np.where(mask_upper, np.sqrt(np.maximum(m - n, 0)), 0.0)
    prod = q[:, None] * q[None, :]

    pos_pos = (q[:, None] > 0) & (q[None, :] > 0) & mask_upper
    neg_neg = (q[:, None] < 0) & (q[None, :] < 0) & mask_upper
    pos_neg = (prod < 0) & mask_upper

    scd_pp = float(np.sum(prod[pos_pos] * weights[pos_pos]) / n_res)
    scd_nn = float(np.sum(prod[neg_neg] * weights[neg_neg]) / n_res)
    scd_pn = float(np.sum(prod[pos_neg] * weights[pos_neg]) / n_res)
    return dict(
        scd_pos_pos=scd_pp,
        scd_neg_neg=scd_nn,
        scd_same_sign=scd_pp + scd_nn,
        scd_pos_neg=scd_pn,
        scd_total=scd_pp + scd_nn + scd_pn,
    )


def _charge_matrix(mask, length_scale):
    positions = np.arange(len(mask), dtype=float)
    kernel = np.exp(-np.abs(positions[:, None] - positions[None, :]) / float(length_scale))
    return mask[:, None] * kernel * mask[None, :]


def _spectral_features(matrix, prefix):
    values, vectors = np.linalg.eigh(matrix)
    order = np.argsort(values)[::-1]
    values = np.clip(values[order], 0.0, None)
    vectors = vectors[:, order]
    trace = float(values.sum())
    if trace <= 0.0:
        raise ValueError(f"{prefix} charge matrix has non-positive trace")
    probs = values / trace
    nonzero = probs[probs > 0]
    entropy = float(-np.sum(nonzero * np.log(nonzero)))
    eff_rank = float(np.exp(entropy))
    principal = vectors[:, 0]
    weights = principal**2
    weights = weights / weights.sum()
    axis = np.linspace(0.0, 1.0, matrix.shape[0]) if matrix.shape[0] > 1 else np.array([0.0])
    center = float(np.sum(axis * weights))
    spread = float(np.sqrt(np.sum((axis - center) ** 2 * weights)))
    ipr = float(np.sum(weights**2))
    return (
        dict(
            trace=trace,
            lambda_max=float(values[0]),
            lambda_fraction=float(values[0] / trace),
            effective_rank=eff_rank,
            spectral_entropy=entropy,
            principal_ipr=ipr,
            principal_center=center,
            principal_spread=spread,
        ),
        weights,
    )


def _cross_svd_features(pos_positions, neg_positions, seq_len, length_scale):
    distances = np.abs(pos_positions[:, None] - neg_positions[None, :])
    matrix = np.exp(-distances / float(length_scale))
    left, svals, right_t = np.linalg.svd(matrix, full_matrices=False)
    energy = svals**2
    total_energy = float(energy.sum())
    probs = energy / total_energy
    nonzero = probs[probs > 0]
    entropy = float(-np.sum(nonzero * np.log(nonzero)))
    left_w = left[:, 0] ** 2
    right_w = right_t[0, :] ** 2
    left_w /= left_w.sum()
    right_w /= right_w.sum()
    denom = max(seq_len - 1, 1)
    pos_center = float(np.sum((pos_positions / denom) * left_w))
    neg_center = float(np.sum((neg_positions / denom) * right_w))
    return dict(
        cross_mean=float(matrix.mean()),
        cross_sigma_max=float(svals[0]),
        cross_sigma_energy_fraction=float(energy[0] / total_energy),
        cross_effective_rank=float(np.exp(entropy)),
        cross_positive_principal_center=pos_center,
        cross_negative_principal_center=neg_center,
        cross_principal_center_separation=abs(pos_center - neg_center),
    )


def compute(simulation, centered_root, output, parameters):
    source = io.centered_source(simulation, centered_root)
    sequence, n_chains = read_sequence(source)
    q = np.array([CHARGE.get(aa, 0) for aa in sequence], dtype=int)
    pos = (q > 0).astype(float)
    neg = (q < 0).astype(float)
    pos_idx = np.flatnonzero(q > 0)
    neg_idx = np.flatnonzero(q < 0)
    if len(pos_idx) == 0 or len(neg_idx) == 0:
        raise ValueError("Sequence requires both positive and negative residues")

    scd_parts = compute_scd_decomposition(q)
    pd.DataFrame([dict(simulation_id=simulation["id"], **scd_parts)]).to_csv(
        output / "scd_decomposition.csv", index=False
    )

    scale_rows = []
    residue_table = pd.DataFrame(
        dict(
            residue_index=np.arange(len(sequence), dtype=int),
            residue_1based=np.arange(1, len(sequence) + 1, dtype=int),
            amino_acid=list(sequence),
            charge=q,
        )
    )

    for ell in parameters["length_scales"]:
        pos_feat, pos_w = _spectral_features(_charge_matrix(pos, ell), "positive")
        neg_feat, neg_w = _spectral_features(_charge_matrix(neg, ell), "negative")
        cross_feat = _cross_svd_features(pos_idx, neg_idx, len(sequence), ell)
        f_pos = pos_feat["lambda_fraction"]
        f_neg = neg_feat["lambda_fraction"]
        scale_rows.append(
            dict(
                simulation_id=simulation["id"],
                length_scale=float(ell),
                f_pos=f_pos,
                f_neg=f_neg,
                delta_lambda_fraction=f_pos - f_neg,
                mean_lambda_fraction=0.5 * (f_pos + f_neg),
                pos_effective_rank=pos_feat["effective_rank"],
                neg_effective_rank=neg_feat["effective_rank"],
                pos_principal_ipr=pos_feat["principal_ipr"],
                neg_principal_ipr=neg_feat["principal_ipr"],
                pos_principal_center=pos_feat["principal_center"],
                neg_principal_center=neg_feat["principal_center"],
                pos_principal_spread=pos_feat["principal_spread"],
                neg_principal_spread=neg_feat["principal_spread"],
                **cross_feat,
            )
        )
        tag = f"ell{int(ell) if float(ell).is_integer() else ell}"
        residue_table[f"pos_mode_weight_{tag}"] = pos_w
        residue_table[f"neg_mode_weight_{tag}"] = neg_w

    scale_df = pd.DataFrame(scale_rows)
    scale_df.to_csv(output / "spectrum_by_scale.csv", index=False)
    residue_table.to_csv(output / "residue_modes.csv", index=False)

    return (
        dict(
            sequence=sequence,
            n_chains=n_chains,
            centered=io.centered_reference(simulation, centered_root),
        ),
        scd_parts,
        scale_df,
    )


def summary_records(scd_parts, scale_df, parameters):
    trapz_fn = getattr(np, "trapezoid", np.trapz)
    int_scales = np.asarray(parameters["integration_scales"], dtype=float)
    sub = scale_df.set_index("length_scale").loc[int_scales]
    log_s = np.log(int_scales)
    span = float(log_s[-1] - log_s[0])

    f_pos_ms = float(trapz_fn(sub["f_pos"].to_numpy(), x=log_s) / span)
    f_neg_ms = float(trapz_fn(sub["f_neg"].to_numpy(), x=log_s) / span)
    ams_val = float(trapz_fn(sub["delta_lambda_fraction"].to_numpy(), x=log_s) / span)
    cms_val = float(trapz_fn(sub["mean_lambda_fraction"].to_numpy(), x=log_s) / span)

    base_def = dict(charge=CHARGE, kernel="exp(-|i-j|/ell)", integration_scales=list(int_scales))
    records = [
        files.scalar_record("spectrum.scd_pos_pos", scd_parts["scd_pos_pos"], "dimensionless", definition=base_def),
        files.scalar_record("spectrum.scd_neg_neg", scd_parts["scd_neg_neg"], "dimensionless", definition=base_def),
        files.scalar_record("spectrum.scd_same_sign", scd_parts["scd_same_sign"], "dimensionless", definition=base_def),
        files.scalar_record("spectrum.scd_pos_neg", scd_parts["scd_pos_neg"], "dimensionless", definition=base_def),
        files.scalar_record("spectrum.f_pos_ms", f_pos_ms, "dimensionless", definition=base_def),
        files.scalar_record("spectrum.f_neg_ms", f_neg_ms, "dimensionless", definition=base_def),
        files.scalar_record("spectrum.ams", ams_val, "dimensionless", definition=base_def),
        files.scalar_record("spectrum.cms", cms_val, "dimensionless", definition=base_def),
    ]

    for _, row in scale_df.iterrows():
        ell = int(row["length_scale"]) if float(row["length_scale"]).is_integer() else float(row["length_scale"])
        cond = dict(ell=ell)
        records.extend(
            [
                files.scalar_record("spectrum.f_pos", row["f_pos"], "dimensionless", conditions=cond, definition=base_def),
                files.scalar_record("spectrum.f_neg", row["f_neg"], "dimensionless", conditions=cond, definition=base_def),
                files.scalar_record(
                    "spectrum.delta_lambda_fraction",
                    row["delta_lambda_fraction"],
                    "dimensionless",
                    conditions=cond,
                    definition=base_def,
                ),
            ]
        )

    ref_row = scale_df.set_index("length_scale").loc[float(parameters["reference_scale"])]
    ref_def = base_def | dict(reference_scale=float(parameters["reference_scale"]))
    for key in (
        "pos_principal_center",
        "neg_principal_center",
        "pos_principal_spread",
        "neg_principal_spread",
        "pos_principal_ipr",
        "neg_principal_ipr",
        "cross_sigma_energy_fraction",
        "cross_principal_center_separation",
    ):
        records.append(
            files.scalar_record(f"spectrum.{key}", ref_row[key], "dimensionless", definition=ref_def)
        )
    return records


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    parent = Path(results_root) / simulation["id"] / "03_ams"
    target = files.reset_output(parent, "01_spectrum")
    inputs, scd_parts, scale_df = compute(simulation, centered_root, target, PARAMETERS)
    files.finish(
        target,
        simulation,
        PARAMETERS,
        inputs,
        scalars=summary_records(scd_parts, scale_df, PARAMETERS),
    )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.02_ams.01_spectrum.run").main())
