from importlib import import_module

"""Interchain contacts, partners, networks and paths."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

PARAMETERS = dict(
    contact_cutoff_nm=0.8,
    chunk_frames=25,
    block_frames=50,
    partner_thresholds=[1, 5, 10],
    network_thresholds=[1, 5, 10, 20, 30, 35],
)
OBJECTS = ["central_chain", "dense_chains"]
io = import_module("02_calculate.common.analysis_io")
store = import_module("02_calculate.02_contacts.common.contact_io")


paths = import_module("02_calculate.02_contacts.02_inter.path_distribution")


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


def compute(
    simulation, centered_root, output, parameters, object_name, io, store, source, results_root
):
    io.validate_parameters(parameters | dict(chunk_frames=1))
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
        raise ValueError("No calculated contact frames")
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


def scan(simulation, centered_root, output, parameters, io, results_root):
    (output / "edges").mkdir()
    buffers, pairs_buffer, membership = [], [], []
    chunk_id = 0

    def flush():
        nonlocal chunk_id
        if not buffers:
            return
        np.savez_compressed(
            output / "edges" / f"{chunk_id:06d}.npz",
            frame=np.array([r["frame"] for r in buffers]),
            source_frame=np.array([r["source_frame"] for r in buffers]),
            time_ps=np.array([r["time_ps"] for r in buffers]),
            dense_chains=np.array([r["dense"] for r in buffers]),
            central_chain=np.array([r["central"] for r in buffers]),
            contact_counts=np.array([r["inter"] for r in buffers], dtype=np.int32),
            pair_i=pair_i,
            pair_j=pair_j,
            n_residues=n_residues,
        )
        io.append_csv(membership, output / "chain_selection.csv.gz")
        buffers.clear()
        pairs_buffer.clear()
        membership.clear()
        chunk_id += 1

    for info, xyz, box, chains, atom_chain, residue, com, dense in io.trajectory(
        simulation, centered_root, parameters["chunk_frames"], results_root
    ):
        n_chains, n_residues = chains.shape
        pair_i, pair_j = np.triu_indices(n_chains, 1)
        pairs = io.close_pairs(xyz, box, parameters["contact_cutoff_nm"])
        ci, cj = atom_chain[pairs[:, 0]], atom_chain[pairs[:, 1]]
        same = ci == cj
        a, b = np.minimum(ci[~same], cj[~same]), np.maximum(ci[~same], cj[~same])
        upper = np.bincount(a * n_chains + b, minlength=n_chains * n_chains).reshape(
            n_chains, n_chains
        )
        central, z = store.chain_selection(com, box, dense)
        buffers.append(
            info
            | dict(
                dense=dense,
                central=central,
                inter=upper[pair_i, pair_j],
            )
        )
        membership.extend(
            info
            | dict(
                chain_id=c, dense=bool(dense[c]), central=c == central, z_center_nm=float(z[c])
            )
            for c in range(n_chains)
        )
        if len(buffers) == parameters["chunk_frames"]:
            flush()
            print(f'  {simulation["id"]}: {info["frame"]+1} inter contact frames', flush=True)
    flush()
    if chunk_id == 0:
        raise ValueError("No centered frames")
    phase = io.phase_source(simulation, centered_root, results_root)
    return dict(
        centered=io.centered_reference(simulation, centered_root),
        phase=files.read_json(phase / "metadata.json"),
        n_chains=n_chains,
        n_residues=n_residues,
        selection="central: closest geometric center among dense chains; tie: lowest chain_id",
    )


def summary_records(stage, object_name, parameters):
    stage = Path(stage)
    definition = dict(
        contact_cutoff_nm=parameters["contact_cutoff_nm"], partners="other_dense_chains"
    )
    records = files.records_from_summary(
        stage / "summary.csv",
        "inter",
        {"partner_count": "chains", "contacts_per_partner": "pairs/partner"},
        object=object_name,
        conditions=("threshold",),
        definition=definition,
    )
    # Total contacts include all dense partners and do not depend on partner threshold.
    row = pd.read_csv(stage / "summary.csv").iloc[0]
    records.insert(
        0,
        files.scalar_record(
            "inter.total_contacts",
            row.total_contacts_mean_frames,
            "pairs/chain",
            object=object_name,
            definition=definition,
            sem=row.total_contacts_sem_blocks,
        ),
    )
    if object_name == "dense_chains":
        records += files.records_from_summary(
            stage / "network/summary.csv",
            "network",
            {
                "global_efficiency": "dimensionless",
                "largest_component_fraction": "dimensionless",
            },
            object=object_name,
            conditions=("threshold",),
            definition=definition,
        )
    return records


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    chosen = store.validate_objects(OBJECTS if objects is None else objects)
    parent = store.result_path(simulation, "02_inter", results_root)
    # The common contact scan is replaced, so old object products must not stay complete.
    for obj in ("central_chain", "dense_chains"):
        (parent / obj / "COMPLETE").unlink(missing_ok=True)
    targets = {obj: files.reset_output(parent, obj) for obj in chosen}
    source = files.reset_output(parent, "common")
    io.validate_parameters(PARAMETERS)
    inputs = scan(simulation, centered_root, source, PARAMETERS, io, results_root)
    files.finish(source, simulation, PARAMETERS, inputs)
    for obj, target in targets.items():
        inputs = compute(
            simulation, centered_root, target, PARAMETERS, obj, io, store, source, results_root
        )
        files.finish(
            target,
            simulation,
            PARAMETERS,
            inputs,
            scalars=summary_records(target, obj, PARAMETERS),
        )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.02_contacts.02_inter.run").main())
