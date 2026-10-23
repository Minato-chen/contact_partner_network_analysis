"""One trajectory scan for shared intrachain and interchain contacts."""
import numpy as np
import pandas as pd
import result_files as files


def compute(simulation, centered_root, output, parameters, io, network):
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
        simulation, centered_root, parameters["chunk_frames"]
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
        a, b = np.minimum(ci[~same], cj[~same]), np.maximum(ci[~same], cj[~same])
        upper = np.bincount(a * n_chains + b, minlength=n_chains * n_chains).reshape(
            n_chains, n_chains
        )
        eligible = np.flatnonzero(dense)
        z = io.minimum_image(com[:, 2] - box[2] / 2, box[2])
        central = int(eligible[np.argmin(np.abs(z[eligible]))]) if len(eligible) else -1
        buffers.append(
            info
            | dict(
                dense=dense,
                central=central,
                inter=upper[pair_i, pair_j],
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
            | dict(chain_id=c, dense=bool(dense[c]), central=c == central, z_center_nm=float(z[c]))
            for c in range(n_chains)
        )
        if len(buffers) == parameters["chunk_frames"]:
            flush()
            print(f'  {simulation["id"]}: {info["frame"]+1} shared contact frames', flush=True)
    flush()
    if chunk_id == 0:
        raise ValueError("No centered frames")
    phase = io.phase_source(simulation, centered_root)
    return dict(
        centered=io.centered_reference(simulation, centered_root),
        phase=files.read_json(phase / "metadata.json"),
        n_chains=n_chains,
        n_residues=n_residues,
        selection="central: closest geometric center among dense chains; tie: lowest chain_id",
    )
