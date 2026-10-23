"""Dilute-phase chain-equivalent concentration outside the upstream 2w boundary."""
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
PARAMETERS = dict(dilute_coefficient=2.0, chunk_frames=25, block_frames=50)


def concentration_summary(frame, block_frames, output):
    """All frames enter the mean; only full, equal-duration blocks enter SEM."""
    blocks = frame.groupby(frame.frame // block_frames).csat_mM.agg(["mean", "count"])
    blocks.index.name = "block"
    blocks["used_for_sem"] = blocks["count"] == block_frames
    complete = blocks.loc[blocks.used_for_sem, "mean"]
    sem = (
        float(complete.std(ddof=1) / np.sqrt(len(complete)))
        if len(complete) >= 2
        else None
    )
    summary = dict(
        csat_mM=float(frame.csat_mM.mean()),
        sem_blocks_mM=sem,
        n_frames=len(frame),
        n_full_blocks=len(complete),
        mean_left_mM=float(frame.dilute_left_mM.mean()),
        mean_right_mM=float(frame.dilute_right_mM.mean()),
    )
    blocks.to_csv(output / "block_means.csv")
    pd.DataFrame([summary]).to_csv(output / "summary.csv", index=False)
    return summary


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    stage = files.reset_output(
        Path(results_root) / simulation["id"] / "01_basic", "csat"
    )
    io.validate_parameters(PARAMETERS)
    coefficient = PARAMETERS["dilute_coefficient"]
    if not np.isfinite(coefficient) or coefficient < 1:
        raise ValueError("dilute_coefficient must be finite and >= 1")
    source = io.phase_source(simulation, centered_root, results_root)
    boundary = files.read_json(source / "boundary.json")
    metadata = files.read_json(source / "metadata.json")
    w = float(boundary["right_half_height_z_nm"])
    box0 = np.asarray(boundary["box_nm"], dtype=float)
    width = float(boundary["bin_width_nm"])
    profile = pd.read_csv(source / "density_profile.csv")
    z = profile.z_nm.to_numpy()
    edges = np.linspace(0, box0[2], len(z) + 1)
    if not np.allclose(
        z, (edges[:-1] + edges[1:]) / 2 - box0[2] / 2, atol=1e-5, rtol=0
    ):
        raise ValueError("Upstream density bins do not match the upstream box")
    if not np.isfinite(w) or w <= 0 or coefficient * w >= box0[2] / 2:
        raise ValueError("No dilute sampling region outside the requested boundary")
    left, right = z < -coefficient * w, z > coefficient * w
    if not left.any() or not right.any():
        raise ValueError("Both dilute sampling regions must contain bins")
    selected = left | right
    regions = dict(
        half_height_width_nm=w,
        dilute_start_nm=coefficient * w,
        bin_width_nm=width,
        selected_bin_centers_nm=z[selected].tolist(),
        dilute_volume_nm3=float(selected.sum() * width * np.prod(box0[:2])),
        rule="abs(bin_center_z_nm) > dilute_coefficient * right_half_height_z_nm",
    )
    files.write_json(stage / "regions.json", regions)
    rows = []
    for info, xyz, box, chains, *_ in io.centered_frames(
        simulation, centered_root, PARAMETERS["chunk_frames"]
    ):
        if not np.allclose(box, box0, atol=1e-5, rtol=0):
            raise ValueError("Csat requires the same fixed box as the upstream density")
        counts, _ = np.histogram(np.mod(xyz[:, 2], box[2]), bins=edges)
        conversion = (
            1e27 / 6.02214076e23 / (chains.shape[1] * np.prod(box0[:2]) * width)
        )
        values = counts * conversion
        rows.append(
            info
            | dict(
                csat_mM=float(values[selected].mean()),
                dilute_left_mM=float(values[left].mean()),
                dilute_right_mM=float(values[right].mean()),
                dilute_beads=int(counts[selected].sum()),
            )
        )
    if len(rows) != boundary["n_frames"]:
        raise ValueError("Csat frames do not match the upstream density calculation")
    frame = pd.DataFrame(rows)
    frame.to_csv(stage / "concentration_frames.csv", index=False)
    summary = concentration_summary(frame, PARAMETERS["block_frames"], stage)
    definition = dict(
        method="bead_density_divided_by_chain_length",
        region=regions["rule"],
        dilute_coefficient=coefficient,
        averaging="spatial_bin_mean_then_equal_frame_mean_including_zero_counts",
        sem="sample_SD_of_full_equal_block_means/sqrt(n_full_blocks)",
        block_frames=PARAMETERS["block_frames"],
    )
    files.finish(
        stage,
        simulation,
        PARAMETERS,
        dict(
            phase_boundary=metadata,
            centered=io.centered_reference(simulation, centered_root),
        ),
        scalars=[
            files.scalar_record(
                "csat",
                summary["csat_mM"],
                "mM",
                definition=definition,
                sem=summary["sem_blocks_mM"],
            )
        ],
    )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.01_basic.run").main(task="csat"))
