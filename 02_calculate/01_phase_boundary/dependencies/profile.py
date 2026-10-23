"""Time-mean chain-equivalent density and the legacy right half-height rule."""
import numpy as np
import pandas as pd
import result_files as files


def estimate_boundary(z, concentration, parameters):
    """Use a sampled bin crossing, not a tanh fit or 60%-height crossing.

    ana_v4 estimates the central peak from +/-N bins, then searches right
    for the first bin below half that peak. The left crossing is diagnostic.
    Invalid profiles fail explicitly instead of using ana_v4's fixed fallback.
    """
    z, concentration = np.asarray(z), np.asarray(concentration)
    if len(z) < 3 or concentration.shape != z.shape or not np.isfinite(concentration).all():
        raise ValueError("Density profile must contain finite, aligned bins")
    center = len(z) // 2
    halfwindow = parameters["peak_halfwidth_bins"]
    peak = float(concentration[max(0, center - halfwindow) : center + halfwindow + 1].mean())
    if peak < parameters["min_peak_concentration_mM"]:
        raise ValueError(f"Central density {peak:.6g} mM is below the minimum; no dense boundary")
    level = peak / 2
    right = np.flatnonzero(concentration[center:] < level)
    if not len(right):
        raise ValueError("No right half-height crossing; dense boundary is undefined (no fallback)")
    right_z = float(z[center + right[0]])
    if right_z <= 0:
        raise ValueError("Right half-height width must be positive")
    left = np.flatnonzero(concentration[:center][::-1] < level)
    left_z = float(z[center - 1 - left[0]]) if len(left) else None
    coefficient = parameters["dense_coefficient"]
    return dict(
        method="central_peak_half_height_right_first_bin",
        peak_chain_mM=peak,
        half_height_chain_mM=level,
        right_half_height_z_nm=right_z,
        left_half_height_z_nm=left_z,
        half_width_asymmetry_nm=right_z - abs(left_z) if left_z is not None else None,
        dense_coefficient=coefficient,
        dense_half_width_nm=right_z * coefficient,
        left_crossing_found=left_z is not None,
        fixed_width_fallback=False,
    )


def compute(simulation, centered_root, output, parameters, io):
    mean = m2 = box0 = centers = None
    n_frames = 0
    for info, xyz, box, chains, atom_chain, residue, com in io.centered_frames(
        simulation, centered_root, parameters["chunk_frames"]
    ):
        if box0 is None:
            box0 = box.astype(float)
            n_bins = round(box0[2] / parameters["bin_width_nm"])
            if n_bins < 3 or not np.isclose(
                n_bins * parameters["bin_width_nm"], box0[2], atol=1e-5, rtol=0
            ):
                raise ValueError("z box length must be a multiple of bin_width_nm (within 1e-5 nm)")
            edges = np.linspace(0, box0[2], n_bins + 1)
            centers = (edges[:-1] + edges[1:]) / 2 - box0[2] / 2
            mean = np.zeros(n_bins)
            m2 = np.zeros(n_bins)
            # Beads / chain length -> chain-equivalent number density; nm^3 -> mM.
            conversion = (
                1e27 / 6.02214076e23 / (chains.shape[1] * np.prod(box0[:2]) * (edges[1] - edges[0]))
            )
        if not np.allclose(box, box0, rtol=0, atol=1e-5):
            raise ValueError("Density boundary calculation requires a fixed box")
        counts, _ = np.histogram(np.mod(xyz[:, 2], box[2]), bins=edges)
        if counts.sum() != len(xyz):
            raise ValueError("Density histogram omitted beads")
        values = counts * conversion
        n_frames += 1
        delta = values - mean
        mean += delta / n_frames
        m2 += delta * (values - mean)
    if not n_frames:
        raise ValueError("No centered frames for density calculation")
    boundary = estimate_boundary(centers, mean, parameters)
    boundary.update(
        n_frames=n_frames,
        box_nm=box0.tolist(),
        bin_width_nm=float(edges[1] - edges[0]),
        peak_halfwidth_bins=parameters["peak_halfwidth_bins"],
    )
    pd.DataFrame(
        dict(
            z_nm=centers,
            mean_chain_mM=mean,
            std_chain_mM=np.sqrt(m2 / (n_frames - 1))
            if n_frames > 1
            else np.full_like(mean, np.nan),
        )
    ).to_csv(output / "density_profile.csv", index=False)
    files.write_json(output / "boundary.json", boundary)
