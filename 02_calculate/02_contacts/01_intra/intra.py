from importlib import import_module

"""Intrachain contacts, object statistics, loop size distributions, residue profiles and Rg pairing."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

PARAMETERS = dict(
    contact_cutoff_nm=0.8, exclude_sequence_separation=2, chunk_frames=25, block_frames=50
)
OBJECTS = ["central_chain", "dense_chains"]
io = import_module("02_calculate.common.analysis_io")
store = import_module("02_calculate.02_contacts.common.contact_io")
network = import_module("02_calculate.02_contacts.common.network")


def compute(
    simulation, centered_root, output, parameters, object_name, io, store, source, results_root
):
    rg_source = store.rg_source(simulation, centered_root, io, results_root)
    rg = store.read_rg(rg_source)
    frames, profiles, distributions, chain_buffer = [], [], [], []
    probability = None
    valid_frames = 0
    
    # 新增: 序列间距分布与残基接触度
    loop_dist = None
    res_degrees = None
    total_intra_contacts = 0
    total_selected_chains = 0

    for chunk in sorted((source / "edges").glob("*.npz")):
        with np.load(chunk) as saved:
            n = int(saved["n_residues"])
            if probability is None:
                probability = np.zeros((n, n))
                loop_dist = np.zeros(n, dtype=np.float64)
                res_degrees = np.zeros(n, dtype=np.float64)

            for index in range(len(saved["frame"])):
                info = store.frame_info(simulation, saved, index, parameters["block_frames"])
                selected = store.selected_chains(saved, index, object_name)
                counts = saved["intra_counts"][index, selected]
                degrees = saved["residue_degree"][index, selected]
                graph = saved["intra_graph_stats"][index, selected]
                radii = np.array(
                    [rg.loc[(info["source_frame"], int(c))] for c in selected], dtype=float
                )

                # 收集当前帧所选对象的 intra_edges 并统计短/长程
                edges = saved["intra_edges"]
                frame_edges = edges[edges[:, 0] == info["frame"]]
                sel_edges = frame_edges[np.isin(frame_edges[:, 1], selected)]
                
                # 为每条所选链预先统计短/中/长程
                chain_edge_map = {}
                if len(sel_edges):
                    c_ids = sel_edges[:, 1]
                    deltas = np.abs(sel_edges[:, 2] - sel_edges[:, 3])
                    # 累加到全局成环谱与残基接触度
                    binc = np.bincount(deltas, minlength=n)
                    loop_dist += binc
                    total_intra_contacts += len(sel_edges)
                    np.add.at(res_degrees, sel_edges[:, 2], 1)
                    np.add.at(res_degrees, sel_edges[:, 3], 1)
                    
                    for u_cid in np.unique(c_ids):
                        c_deltas = deltas[c_ids == u_cid]
                        chain_edge_map[u_cid] = (
                            int(np.sum(c_deltas <= 15)),
                            int(np.sum((c_deltas > 15) & (c_deltas <= 30))),
                            int(np.sum(c_deltas > 30)),
                        )

                total_selected_chains += len(selected)

                for k, chain in enumerate(selected):
                    c_id = int(chain)
                    sh, med, lg = chain_edge_map.get(c_id, (0, 0, 0))
                    chain_buffer.append(
                        info
                        | dict(
                            chain_id=c_id,
                            contact_count=int(counts[k]),
                            short_range_contacts=sh,
                            medium_range_contacts=med,
                            long_range_contacts=lg,
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
                    | {
                        k: float(v.mean()) if len(selected) else np.nan
                        for k, v in values.items()
                    }
                )
                participation = (
                    (degrees > 0).mean(axis=0) if len(selected) else np.full(n, np.nan)
                )
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
                    row | dict(degree=k, fraction=float(value))
                    for k, value in enumerate(frequency)
                )
                if len(selected):
                    keep = (edges[:, 0] == info["frame"]) & np.isin(edges[:, 1], selected)
                    i, j = edges[keep, 2], edges[keep, 3]
                    np.add.at(probability, (i, j), 1 / len(selected))
                    np.add.at(probability, (j, i), 1 / len(selected))
                    valid_frames += 1
        if chain_buffer:
            io.append_csv(chain_buffer, output / "chain_contacts.csv.gz")
            chain_buffer.clear()

    if probability is None:
        raise ValueError("No calculated contact frames")

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
                "short_range_contacts",
                "medium_range_contacts",
                "long_range_contacts",
                "rg_nm",
                "largest_component_fraction",
                "global_efficiency",
                "n_residues",
            ]
        ).to_csv(chain_file, compression="gzip", index=False)

    probability /= max(valid_frames, 1)
    np.save(output / "contact_probability.npy", probability)

    # 保存成环尺度谱 P(|i-j|) 与残基接触度轮廓
    p_loop = loop_dist / max(total_intra_contacts, 1)
    mean_res_degree = res_degrees / max(total_selected_chains, 1)
    np.save(output / "loop_size_distribution.npy", p_loop)
    np.save(output / "residue_degree_profile.npy", mean_res_degree)

    frame_data = pd.DataFrame(frames)
    frame_data.to_csv(output / "frame_means.csv", index=False)
    metrics = [
        "contact_count",
        "rg_nm",
        "largest_component_fraction",
        "global_efficiency",
    ]
    io.summarize(frame_data, metrics, ["simulation_id"], output)
    pd.DataFrame(profiles).to_csv(output / "residue_profile.csv.gz", index=False)
    pd.DataFrame(distributions).to_csv(output / "degree_distribution.csv.gz", index=False)
    return dict(
        object=object_name,
        contacts=store.reference(source),
        rg=files.read_json(rg_source / "metadata.json"),
    )


def scan(simulation, centered_root, output, parameters, io, results_root):
    # 底层轨迹扫描部分保持 100% 相同
    edges_dir = output / "edges"
    edges_dir.mkdir(parents=True, exist_ok=True)
    buffers, pairs_buffer, membership = [], [], []
    chunk_id = 0

    def flush():
        nonlocal chunk_id
        if not buffers:
            return
        np.savez_compressed(
            edges_dir / f"{chunk_id:06d}.npz",
            frame=np.array([r["frame"] for r in buffers], dtype=np.int32),
            source_frame=np.array([r["source_frame"] for r in buffers], dtype=np.int32),
            time_ps=np.array([r["time_ps"] for r in buffers]),
            dense_chains=np.array([r["dense"] for r in buffers]),
            central_chain=np.array([r["central"] for r in buffers]),
            intra_counts=np.array([r["intra"] for r in buffers], dtype=np.int32),
            residue_degree=np.array([r["degree"] for r in buffers], dtype=np.int32),
            intra_graph_stats=np.array([r["graph"] for r in buffers]),
            intra_edges=np.concatenate(pairs_buffer),
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
        separated = (
            np.abs(residue[pairs[:, 0]] - residue[pairs[:, 1]])
            > parameters["exclude_sequence_separation"]
        )
        intra = pairs[same & separated]
        chain_ids = atom_chain[intra[:, 0]]
        i, j = residue[intra[:, 0]], residue[intra[:, 1]]
        counts = np.bincount(chain_ids, minlength=n_chains)
        degrees = np.zeros((n_chains, n_residues), dtype=np.int32)
        np.add.at(degrees, (chain_ids, i), 1)
        np.add.at(degrees, (chain_ids, j), 1)
        stats, _ = network.chain_statistics(chain_ids, i, j, n_chains, n_residues)
        central, z = store.chain_selection(com, box, dense)
        buffers.append(
            info
            | dict(
                dense=dense,
                central=central,
                intra=counts,
                degree=degrees,
                graph=stats,
            )
        )
        pairs_buffer.append(
            np.column_stack((np.full(len(i), info["frame"]), chain_ids, i, j)).astype(np.int32)
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
            print(f'  {simulation["id"]}: {info["frame"]+1} intra contact frames', flush=True)
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
    definition = {
        key: parameters[key] for key in ("contact_cutoff_nm", "exclude_sequence_separation")
    }
    
    # 提取原有标量记录
    records = files.records_from_summary(
        Path(stage) / "summary.csv",
        "intra",
        {
            "contact_count": "pairs/chain",
            "rg_nm": "nm",
            "largest_component_fraction": "dimensionless",
            "global_efficiency": "dimensionless",
        },
        object=object_name,
        definition=definition,
    )
    
    # 新增: 长程闭环占比标量 (Long-range loop fraction >30)
    p_loop_file = Path(stage) / "loop_size_distribution.npy"
    if p_loop_file.is_file():
        p_loop = np.load(p_loop_file)
        f_long = float(p_loop[31:].sum())
        records.append(
            files.scalar_record(
                "intra.long_range_fraction",
                f_long,
                "dimensionless",
                object=object_name,
                definition=definition | dict(long_range_threshold=30),
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
    parent = store.result_path(simulation, "01_intra", results_root)
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
    raise SystemExit(import_module("02_calculate.02_contacts.01_intra.run").main())
