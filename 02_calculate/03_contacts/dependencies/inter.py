"""Selected-chain contacts and partners; whole dense-chain network statistics."""
import numpy as np
import pandas as pd
import project

HERE = project.ROOT / "02_calculate/03_contacts/dependencies"
paths = project.load_module("inter_path_distribution", HERE / "path_distribution.py")


def graph_stats(adjacency):
    # Identical definitions to the existing intermolecular implementation.
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path, connected_components

    n = len(adjacency)
    if not n:
        return np.nan, np.nan
    if n == 1:
        return np.nan, 1.0
    graph = csr_matrix(adjacency)
    distances = shortest_path(graph, directed=False, unweighted=True)
    inverse = np.zeros_like(distances)
    reachable = np.isfinite(distances) & (distances > 0)
    inverse[reachable] = 1 / distances[reachable]
    _, labels = connected_components(graph, directed=False)
    return float(inverse.sum() / (n * (n - 1))), float(np.bincount(labels).max() / n)


def compute(simulation, centered_root, output, parameters, object_name, io, store):
    io.validate_parameters(parameters | dict(chunk_frames=1))
    source = store.prepared_source(simulation, centered_root, io)
    frame_rows, networks, distribution_rows, chain_buffer = [], [], [], []
    for chunk in sorted((source / "edges").glob("*.npz")):
        with np.load(chunk) as saved:
            n = saved["dense_chains"].shape[1]
            pair_i, pair_j = saved["pair_i"], saved["pair_j"]
            for index in range(len(saved["frame"])):
                info = store.frame_info(simulation, saved, index, parameters["block_frames"])
                selected = store.selected_chains(saved, index, object_name)
                dense = saved["dense_chains"][index]
                counts = np.zeros((n, n), dtype=np.int32)
                counts[pair_i, pair_j] = saved["contact_counts"][index]
                counts[pair_j, pair_i] = saved["contact_counts"][index]
                counts[~dense, :] = 0
                counts[:, ~dense] = 0
                totals = counts[selected].sum(axis=1)
                for threshold in parameters["partner_thresholds"]:
                    partners = (counts[selected] >= threshold).sum(axis=1)
                    ratios = np.divide(
                        totals, partners, out=np.full(len(selected), np.nan), where=partners > 0
                    )
                    for k, chain in enumerate(selected):
                        chain_buffer.append(
                            info
                            | dict(
                                chain_id=int(chain),
                                threshold=threshold,
                                total_contacts=int(totals[k]),
                                partner_count=int(partners[k]),
                                contacts_per_partner=float(ratios[k]),
                            )
                        )
                    frame_rows.append(
                        info
                        | dict(
                            threshold=threshold,
                            n_chains=len(selected),
                            n_dense=int(dense.sum()),
                            total_contacts=float(totals.mean()) if len(selected) else np.nan,
                            partner_count=float(partners.mean()) if len(selected) else np.nan,
                            contacts_per_partner=float(np.nanmean(ratios))
                            if np.isfinite(ratios).any()
                            else np.nan,
                            n_valid_contacts_per_partner=int(np.isfinite(ratios).sum()),
                        )
                    )
                    frequency = (
                        np.bincount(partners, minlength=n) / len(selected)
                        if len(selected)
                        else np.full(n, np.nan)
                    )
                    distribution_rows.extend(
                        info | dict(threshold=threshold, partner_count=k, fraction=float(v))
                        for k, v in enumerate(frequency)
                    )
                if object_name == "dense_chains":
                    matrix = counts[np.ix_(dense, dense)]
                    for threshold in parameters["network_thresholds"]:
                        efficiency, largest = graph_stats(matrix >= threshold)
                        networks.append(
                            info
                            | dict(
                                threshold=threshold,
                                n_dense=int(dense.sum()),
                                global_efficiency=efficiency,
                                largest_component_fraction=largest,
                            )
                        )
        if chain_buffer:
            io.append_csv(chain_buffer, output / "chain_partners.csv.gz")
            chain_buffer.clear()
    if not frame_rows:
        raise ValueError("No prepared contact frames")
    if not (output / "chain_partners.csv.gz").exists():
        pd.DataFrame(
            columns=[
                "simulation_id",
                "frame",
                "source_frame",
                "time_ps",
                "block",
                "chain_id",
                "threshold",
                "total_contacts",
                "partner_count",
                "contacts_per_partner",
            ]
        ).to_csv(output / "chain_partners.csv.gz", index=False)
    frame = pd.DataFrame(frame_rows)
    frame.to_csv(output / "frame_means.csv", index=False)
    io.summarize(
        frame,
        ["total_contacts", "partner_count", "contacts_per_partner"],
        ["simulation_id", "threshold"],
        output,
    )
    distribution = pd.DataFrame(distribution_rows)
    distribution.to_csv(output / "partner_distribution_frames.csv.gz", index=False)
    destination = output / "partner_distribution"
    store.distribution_summary(
        distribution, ["simulation_id", "threshold", "partner_count"], destination, io
    )
    overall = pd.read_csv(destination / "summary.csv").rename(
        columns={"fraction_mean_frames": "mean", "fraction_count_frames": "count"}
    )
    overall.to_csv(output / "partner_distribution.csv", index=False)
    pd.read_csv(destination / "block_means.csv").to_csv(
        output / "partner_distribution_blocks.csv", index=False
    )
    if object_name == "dense_chains":
        network_frame = pd.DataFrame(networks)
        network_frame.to_csv(output / "network_frames.csv", index=False)
        (output / "network").mkdir()
        io.summarize(
            network_frame,
            ["global_efficiency", "largest_component_fraction"],
            ["simulation_id", "threshold"],
            output / "network",
        )
        path_summary, path_blocks = paths.summarize_paths(
            source, parameters["network_thresholds"], parameters["block_frames"]
        )
        path_summary.insert(0, "simulation_id", simulation["id"])
        path_blocks.insert(0, "simulation_id", simulation["id"])
        path_summary.to_csv(output / "path_distribution.csv", index=False)
        path_blocks.to_csv(output / "path_distribution_blocks.csv", index=False)
    return store.object_inputs(source) | dict(
        object=object_name, contact_partners="other dense chains"
    )
