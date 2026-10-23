"""Object-specific intrachain counts, participation profiles and Rg pairing."""
import numpy as np
import pandas as pd


def compute(simulation, centered_root, output, parameters, object_name, io, store):
    source = store.prepared_source(simulation, centered_root, io)
    rg_source = store.rg_source(simulation, centered_root, io)
    rg = store.read_rg(rg_source)
    frames, profiles, distributions, chain_buffer = [], [], [], []
    probability = None
    valid_frames = 0
    for chunk in sorted((source / "edges").glob("*.npz")):
        with np.load(chunk) as saved:
            n = int(saved["n_residues"])
            if probability is None:
                probability = np.zeros((n, n))
            for index in range(len(saved["frame"])):
                info = store.frame_info(simulation, saved, index, parameters["block_frames"])
                selected = store.selected_chains(saved, index, object_name)
                counts = saved["intra_counts"][index, selected]
                degrees = saved["residue_degree"][index, selected]
                graph = saved["intra_graph_stats"][index, selected]
                radii = np.array(
                    [rg.loc[(info["source_frame"], int(c))] for c in selected], dtype=float
                )
                for k, chain in enumerate(selected):
                    chain_buffer.append(
                        info
                        | dict(
                            chain_id=int(chain),
                            contact_count=int(counts[k]),
                            rg_nm=float(radii[k]),
                            largest_component_fraction=float(graph[k, 0]),
                            global_efficiency=float(graph[k, 1]),
                            n_residues=n,
                        )
                    )
                row = info | dict(n_chains=len(selected))
                values = dict(
                    contact_count=counts,
                    rg_nm=radii,
                    largest_component_fraction=graph[:, 0],
                    global_efficiency=graph[:, 1],
                )
                frames.append(
                    row
                    | {k: float(v.mean()) if len(selected) else np.nan for k, v in values.items()}
                )
                # Probability of at least one contact, NOT the sum of pair probabilities.
                participation = (degrees > 0).mean(axis=0) if len(selected) else np.full(n, np.nan)
                profiles.extend(
                    row | dict(residue_index=i, fraction=float(value))
                    for i, value in enumerate(participation)
                )
                frequency = (
                    np.bincount(degrees.ravel(), minlength=n) / degrees.size
                    if len(selected)
                    else np.full(n, np.nan)
                )
                distributions.extend(
                    row | dict(degree=k, fraction=float(value)) for k, value in enumerate(frequency)
                )
                if len(selected):
                    edges = saved["intra_edges"]
                    keep = (edges[:, 0] == info["frame"]) & np.isin(edges[:, 1], selected)
                    i, j = edges[keep, 2], edges[keep, 3]
                    np.add.at(probability, (i, j), 1 / len(selected))
                    np.add.at(probability, (j, i), 1 / len(selected))
                    valid_frames += 1
        if chain_buffer:
            io.append_csv(chain_buffer, output / "chain_contacts.csv.gz")
            chain_buffer.clear()
    if probability is None:
        raise ValueError("No prepared contact frames")
    chain_file = output / "chain_contacts.csv.gz"
    if not chain_file.exists():
        pd.DataFrame(
            columns=[
                "simulation_id",
                "frame",
                "source_frame",
                "time_ps",
                "block",
                "chain_id",
                "contact_count",
                "rg_nm",
                "largest_component_fraction",
                "global_efficiency",
                "n_residues",
            ]
        ).to_csv(chain_file, index=False)
    frame = pd.DataFrame(frames)
    frame.to_csv(output / "frame_means.csv", index=False)
    io.summarize(
        frame,
        ["contact_count", "rg_nm", "largest_component_fraction", "global_efficiency"],
        ["simulation_id"],
        output,
    )
    for name, rows, key in [
        ("residue_probability", profiles, "residue_index"),
        ("degree_distribution", distributions, "degree"),
    ]:
        frame_distribution = pd.DataFrame(rows)
        frame_distribution.to_csv(output / (name + "_frames.csv.gz"), index=False)
        store.distribution_summary(frame_distribution, ["simulation_id", key], output / name, io)
    np.savez_compressed(
        output / "contact_map.npz",
        probability=probability / valid_frames
        if valid_frames
        else np.full_like(probability, np.nan),
        valid_frames=valid_frames,
        residue_index=np.arange(n),
    )
    return store.object_inputs(source, rg_source) | dict(object=object_name)
