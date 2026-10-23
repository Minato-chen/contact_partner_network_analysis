"""Intrachain contacts, instantaneous networks and frame-equal distributions."""
import numpy as np
import pandas as pd


def compute(simulation, centered_root, output, parameters, io, network):
    frame_rows = []
    maps = None
    valid_frames = dict(all=0, dense=0)
    buffer = []
    distribution_rows = []
    edge_buffer, edge_frames, dense_frames = [], [], []
    edge_directory = output / 'edges'
    edge_directory.mkdir()
    chunk = 0

    def flush_edges():
        nonlocal chunk
        if not edge_frames:
            return
        edges = np.concatenate(edge_buffer, axis=0)
        np.savez_compressed(edge_directory / f'{chunk:06d}.npz',
            frame=np.array(edge_frames), dense_chains=np.array(dense_frames),
            edges=edges, n_residues=n_residues)
        chunk += 1
        edge_buffer.clear()
        edge_frames.clear()
        dense_frames.clear()

    for info, xyz, box, chains, atom_chain, residue, com, dense in io.trajectory(
            simulation, centered_root, parameters['chunk_frames']):
        n_chains, n_residues = chains.shape
        if maps is None:
            maps = {scope: np.zeros((n_residues, n_residues), dtype=float) for scope in valid_frames}
        pairs = io.close_pairs(xyz, box, parameters['contact_cutoff_nm'])
        same = atom_chain[pairs[:, 0]] == atom_chain[pairs[:, 1]]
        separated = np.abs(residue[pairs[:, 0]] - residue[pairs[:, 1]]) > parameters['exclude_sequence_separation']
        pairs = pairs[same & separated]
        chain_ids = atom_chain[pairs[:, 0]]
        i, j = residue[pairs[:, 0]], residue[pairs[:, 1]]
        counts = np.bincount(chain_ids, minlength=n_chains)
        metrics, histograms = network.chain_statistics(chain_ids, i, j, n_chains, n_residues)
        edge_buffer.append(np.column_stack((np.full(len(i), info['frame']), chain_ids, i, j)).astype(np.int32))
        edge_frames.append(info['frame'])
        dense_frames.append(dense.copy())
        block = info['frame'] // parameters['block_frames']
        for chain in range(n_chains):
            buffer.append(info | dict(block=block, chain=chain, dense=bool(dense[chain]),
                                      contact_count=int(counts[chain]),
                                      largest_component_fraction=metrics[chain, 0],
                                      global_efficiency=metrics[chain, 1]))
        for scope, selected in [('all', np.ones(n_chains, bool)), ('dense', dense)]:
            n_selected = int(selected.sum())
            frame_rows.append(info | dict(block=block, scope=scope, n_chains=n_selected,
                contact_count=float(counts[selected].mean()) if n_selected else np.nan,
                largest_component_fraction=float(metrics[selected, 0].mean()) if n_selected else np.nan,
                global_efficiency=float(metrics[selected, 1].mean()) if n_selected else np.nan))
            frequency = histograms[selected].mean(axis=0) if n_selected else np.full(n_residues, np.nan)
            for degree, fraction in enumerate(frequency):
                distribution_rows.append(info | dict(block=block, scope=scope,
                    n_chains=n_selected, degree=degree, fraction=fraction))
            if n_selected:
                keep = selected[chain_ids]
                np.add.at(maps[scope], (i[keep], j[keep]), 1 / n_selected)
                np.add.at(maps[scope], (j[keep], i[keep]), 1 / n_selected)
                valid_frames[scope] += 1
        if (info['frame'] + 1) % parameters['chunk_frames'] == 0:
            io.append_csv(buffer, output / 'chain_contacts.csv.gz')
            buffer.clear()
            flush_edges()
            print(f'  {simulation["id"]}: {info["frame"] + 1} intramolecular frames', flush=True)
    if buffer:
        io.append_csv(buffer, output / 'chain_contacts.csv.gz')
    flush_edges()
    if maps is None:
        raise ValueError('No frames')
    frame = pd.DataFrame(frame_rows)
    frame.to_csv(output / 'frame_means.csv', index=False)
    io.summarize(frame, ['contact_count', 'largest_component_fraction', 'global_efficiency'], ['simulation_id', 'scope'],
                 output)
    arrays = {}
    for scope in maps:
        arrays[f'{scope}_probability'] = (maps[scope] / valid_frames[scope] if valid_frames[scope]
                                          else np.full_like(maps[scope], np.nan))
        arrays[f'{scope}_valid_frames'] = valid_frames[scope]
    arrays['residue_index'] = np.arange(n_residues)
    np.savez_compressed(output / 'contact_map.npz', **arrays)

    distribution = pd.DataFrame(distribution_rows)
    distribution.to_csv(output / 'degree_distribution_frames.csv.gz', index=False)
    destination = output / 'degree_distribution'
    destination.mkdir()
    io.summarize(distribution, ['fraction'], ['simulation_id', 'scope', 'degree'], destination)
