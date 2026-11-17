"""Couple sequence charge modes (A_MS) with intramolecular and intermolecular contacts."""
from importlib import import_module
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

PARAMETERS = dict(block_frames=50, long_range_separation=12, reference_scale=10.0)
OBJECTS = ["central_chain", "dense_chains"]
io = import_module("02_calculate.common.analysis_io")
store = import_module("02_calculate.01_contacts.common.contact_io")


def compute(
    simulation, centered_root, output, parameters, object_name, io, store, results_root
):
    sid = simulation["id"]
    spectrum_dir = files.completed(
        Path(results_root) / sid / "03_ams/01_spectrum",
        hint=f"请先运行: python 02_calculate/02_ams/01_spectrum/run.py --simulation {sid}",
    )
    intra = files.completed(store.result_path(simulation, "01_intra", results_root) / object_name)
    inter = files.completed(store.result_path(simulation, "02_inter", results_root) / object_name)
    joint = files.completed(store.result_path(simulation, "03_joint", results_root) / object_name)

    current_centered = io.centered_reference(simulation, centered_root)
    spec_meta = store.reference(spectrum_dir)
    if spec_meta["inputs"]["centered"] != current_centered:
        raise files.InputUnavailable("Spectrum centered input changed; rerun 01_spectrum")

    joint_meta = store.reference(joint)
    if (
        joint_meta["inputs"]["intramolecular"] != store.reference(intra)
        or joint_meta["inputs"]["intermolecular"] != store.reference(inter)
    ):
        raise files.InputUnavailable("Joint contact upstream changed; rerun 01_contacts")

    modes = pd.read_csv(spectrum_dir / "residue_modes.csv")
    q = modes["charge"].to_numpy(dtype=int)
    n_res = len(q)
    n_neg = int(np.sum(q < 0))
    n_pos = int(np.sum(q > 0))
    if n_neg == 0 or n_pos == 0:
        raise ValueError("Coupling requires both positive and negative residues")

    ref_tag = f"ell{int(parameters['reference_scale']) if float(parameters['reference_scale']).is_integer() else parameters['reference_scale']}"
    w_pos = modes[f"pos_mode_weight_{ref_tag}"].to_numpy(dtype=float)
    w_neg = modes[f"neg_mode_weight_{ref_tag}"].to_numpy(dtype=float)

    prob = np.load(intra / "contact_probability.npy")
    if prob.shape != (n_res, n_res):
        raise ValueError("Intramolecular contact matrix dimension mismatch with sequence")

    exclude_sep = int(store.reference(intra)["parameters"]["exclude_sequence_separation"])
    lr_sep = int(parameters["long_range_separation"])
    idx = np.arange(n_res)
    sep = np.abs(idx[:, None] - idx[None, :])
    valid_mask = sep > exclude_sep
    lr_mask = sep >= lr_sep
    sr_mask = valid_mask & (sep < lr_sep)

    prod = q[:, None] * q[None, :]
    opp_mask = prod < 0
    same_mask = prod > 0
    neu_mask = prod == 0

    intra_total = float(0.5 * prob[valid_mask].sum())
    intra_attract = float(0.5 * prob[valid_mask & opp_mask].sum())
    intra_attract_sr = float(0.5 * prob[sr_mask & opp_mask].sum())
    intra_attract_lr = float(0.5 * prob[lr_mask & opp_mask].sum())
    intra_repulse = float(0.5 * prob[valid_mask & same_mask].sum())
    intra_neutral = float(0.5 * prob[valid_mask & neu_mask].sum())
    intra_attract_frac = float(intra_attract / intra_total) if intra_total > 0 else np.nan

    res_total_deg = (prob * valid_mask).sum(axis=1)
    res_attract_deg = (prob * (valid_mask & opp_mask)).sum(axis=1)
    res_repulse_deg = (prob * (valid_mask & same_mask)).sum(axis=1)
    res_neutral_deg = (prob * (valid_mask & neu_mask)).sum(axis=1)

    res_prof_df = pd.read_csv(intra / "residue_profile.csv.gz")
    res_part_mean = (
        res_prof_df.groupby("residue_index")["fraction"]
        .mean()
        .reindex(np.arange(n_res))
        .to_numpy(dtype=float)
    )

    residue_coupling = pd.DataFrame(
        dict(
            simulation_id=sid,
            residue_index=modes["residue_index"],
            residue_1based=modes["residue_1based"],
            amino_acid=modes["amino_acid"],
            charge=q,
            pos_mode_weight=w_pos,
            neg_mode_weight=w_neg,
            intra_degree=res_total_deg,
            intra_attract_degree=res_attract_deg,
            intra_repulse_degree=res_repulse_deg,
            intra_neutral_degree=res_neutral_deg,
            intra_participation=res_part_mean,
        )
    )
    residue_coupling.to_csv(output / "residue_coupling_profile.csv", index=False)

    paired = pd.read_csv(joint / "paired_chain_frames.csv.gz")
    valid_zipper = (paired["partner_count"] > 0) & (paired["attract_contacts"] > 0)
    paired["max_partner_attract_share"] = np.where(
        valid_zipper,
        paired["max_attract_contacts_per_partner"] / paired["attract_contacts"],
        np.nan,
    )
    paired["partner_per_attract_contact"] = np.where(
        paired["attract_contacts"] > 0,
        paired["partner_count"] / paired["attract_contacts"],
        np.nan,
    )
    paired["partner_per_neg_residue"] = paired["partner_count"] / float(n_neg)
    paired["inter_attract_to_intra_ratio"] = np.where(
        paired["contact_count"] > 0,
        paired["attract_contacts"] / (2.0 * paired["contact_count"]),
        np.nan,
    )

    metrics = [
        "max_partner_attract_share",
        "partner_per_attract_contact",
        "partner_per_neg_residue",
        "inter_attract_to_intra_ratio",
    ]
    grid = pd.read_csv(inter / "frame_means.csv")[
        ["simulation_id", "frame", "source_frame", "time_ps", "threshold", "n_chains"]
    ]
    means = (
        paired.groupby(["simulation_id", "source_frame", "threshold"])[metrics]
        .mean()
        .reset_index()
    )
    frames = grid.merge(
        means,
        on=["simulation_id", "source_frame", "threshold"],
        how="left",
        validate="one_to_one",
    )
    frames["block"] = frames["frame"] // parameters["block_frames"]
    frames.to_csv(output / "frame_means.csv", index=False)
    io.summarize(frames, metrics, ["simulation_id", "threshold"], output)

    inter_summary = pd.read_csv(inter / "summary.csv").iloc[0]
    inter_attract_mean = float(inter_summary["attract_contacts_mean_frames"])
    denom_att = 2.0 * intra_attract + inter_attract_mean
    inter_attract_partition = float(inter_attract_mean / denom_att) if denom_att > 0 else np.nan

    decomp_row = dict(
        simulation_id=sid,
        object=object_name,
        intra_total_contacts=intra_total,
        intra_attract_contacts=intra_attract,
        intra_short_range_attract_contacts=intra_attract_sr,
        intra_long_range_attract_contacts=intra_attract_lr,
        intra_repulse_contacts=intra_repulse,
        intra_neutral_contacts=intra_neutral,
        intra_attract_fraction=intra_attract_frac,
        inter_attract_contacts=inter_attract_mean,
        inter_attract_partition=inter_attract_partition,
        pos_intra_degree=float(res_total_deg[q > 0].mean()),
        neg_intra_degree=float(res_total_deg[q < 0].mean()),
        neutral_intra_degree=float(res_total_deg[q == 0].mean()) if np.any(q == 0) else np.nan,
        pos_intra_attract_degree=float(res_attract_deg[q > 0].mean()),
        neg_intra_attract_degree=float(res_attract_deg[q < 0].mean()),
        pos_mode_intra_degree=float(np.sum(w_pos * res_total_deg)),
        neg_mode_intra_degree=float(np.sum(w_neg * res_total_deg)),
        pos_mode_intra_attract_degree=float(np.sum(w_pos * res_attract_deg)),
        neg_mode_intra_attract_degree=float(np.sum(w_neg * res_attract_deg)),
    )
    pd.DataFrame([decomp_row]).to_csv(output / "intra_charge_decomposition.csv", index=False)

    return (
        dict(
            object=object_name,
            spectrum=spec_meta,
            intramolecular=store.reference(intra),
            intermolecular=store.reference(inter),
            joint=joint_meta,
        ),
        decomp_row,
    )


def summary_records(stage, object_name, inputs, decomp_row, parameters):
    stage = Path(stage)
    intra_params = inputs["intramolecular"]["parameters"]
    definition = dict(
        contact_cutoff_nm=intra_params["contact_cutoff_nm"],
        exclude_sequence_separation=intra_params["exclude_sequence_separation"],
        long_range_separation=parameters["long_range_separation"],
        reference_scale=parameters["reference_scale"],
    )

    static_units = dict(
        intra_attract_contacts="pairs/chain",
        intra_short_range_attract_contacts="pairs/chain",
        intra_long_range_attract_contacts="pairs/chain",
        intra_repulse_contacts="pairs/chain",
        intra_neutral_contacts="pairs/chain",
        intra_attract_fraction="dimensionless",
        inter_attract_partition="dimensionless",
        pos_intra_degree="contacts/residue",
        neg_intra_degree="contacts/residue",
        neutral_intra_degree="contacts/residue",
        pos_intra_attract_degree="contacts/residue",
        neg_intra_attract_degree="contacts/residue",
        pos_mode_intra_degree="contacts/residue",
        neg_mode_intra_degree="contacts/residue",
        pos_mode_intra_attract_degree="contacts/residue",
        neg_mode_intra_attract_degree="contacts/residue",
    )
    records = [
        files.scalar_record(
            f"coupling.{key}",
            decomp_row[key],
            unit,
            object=object_name,
            definition=definition,
        )
        for key, unit in static_units.items()
    ]

    summary_df = pd.read_csv(stage / "summary.csv")
    first_row = summary_df.iloc[0]
    records.append(
        files.scalar_record(
            "coupling.inter_attract_to_intra_ratio",
            first_row["inter_attract_to_intra_ratio_mean_frames"],
            "dimensionless",
            object=object_name,
            definition=definition,
            sem=first_row.get("inter_attract_to_intra_ratio_sem_blocks"),
        )
    )

    threshold_units = dict(
        max_partner_attract_share="dimensionless",
        partner_per_attract_contact="chains/pair",
        partner_per_neg_residue="chains/residue",
    )
    records.extend(
        files.records_from_summary(
            stage / "summary.csv",
            "coupling",
            threshold_units,
            object=object_name,
            conditions=("threshold",),
            definition=definition,
        )
    )
    return records


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    chosen = store.validate_objects(OBJECTS if objects is None else objects)
    parent = Path(results_root) / simulation["id"] / "03_ams/02_coupling"
    for obj in chosen:
        target = files.reset_output(parent, obj)
        io.validate_parameters(PARAMETERS | dict(chunk_frames=1))
        inputs, decomp_row = compute(
            simulation, centered_root, target, PARAMETERS, obj, io, store, results_root
        )
        files.finish(
            target,
            simulation,
            PARAMETERS,
            inputs,
            scalars=summary_records(target, obj, inputs, decomp_row, PARAMETERS),
        )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.02_ams.02_coupling.run").main())
