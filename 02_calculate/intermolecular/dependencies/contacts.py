"""Interchain contacts, partner distributions and network metrics; no plotting."""
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path


def graph_stats(adjacency):
    n = len(adjacency)
    if not n:
        return np.nan, np.nan
    if n == 1:
        return np.nan, 1.0
    distances = shortest_path(csr_matrix(adjacency), directed=False, unweighted=True)
    inverse = np.zeros_like(distances)
    valid = np.isfinite(distances) & (distances > 0)
    inverse[valid] = 1 / distances[valid]
    _, labels = connected_components(csr_matrix(adjacency), directed=False)
    return float(inverse.sum() / (n * (n - 1))), float(np.bincount(labels).max() / n)


def compute(simulation, centered_root, output, parameters, io):
    frame_rows, networks = [], []
    distribution_buffer = []
    distribution_sums = {}
    chain_buffer, contact_buffer, dense_buffer, frame_buffer = [], [], [], []
    chunk_number = 0
    pair_i = pair_j = None
    (output / 'edges').mkdir()

    def flush_edges():
        nonlocal chunk_number
        if contact_buffer:
            np.savez_compressed(output / 'edges' / f'{chunk_number:06d}.npz',
                frame=np.array(frame_buffer), contact_counts=np.array(contact_buffer, dtype=np.int32),
                dense_chains=np.array(dense_buffer), pair_i=pair_i, pair_j=pair_j)
            chunk_number += 1
            contact_buffer.clear()
            dense_buffer.clear()
            frame_buffer.clear()

    for info, xyz, box, chains, atom_chain, residue, com, dense in io.trajectory(
            simulation, centered_root, parameters['chunk_frames']):
        n = len(chains)
        if pair_i is None:
            pair_i, pair_j = np.triu_indices(n, 1)
        pairs = io.close_pairs(xyz, box, parameters['contact_cutoff_nm'])
        ci, cj = atom_chain[pairs[:, 0]], atom_chain[pairs[:, 1]]
        keep = ci != cj
        lo, hi = np.minimum(ci[keep], cj[keep]), np.maximum(ci[keep], cj[keep])
        upper = np.bincount(lo * n + hi, minlength=n * n).reshape(n, n)
        counts = upper + upper.T
        contact_buffer.append(counts[pair_i, pair_j])
        dense_buffer.append(dense.copy())
        frame_buffer.append(info['frame'])
        counts[~dense, :] = 0
        counts[:, ~dense] = 0
        selected = np.flatnonzero(dense)
        totals = counts[selected].sum(axis=1)
        block = info['frame'] // parameters['block_frames']
        for threshold in parameters['partner_thresholds']:
            partners = (counts[selected] >= threshold).sum(axis=1)
            ratios = np.divide(totals, partners, out=np.full(len(partners), np.nan), where=partners > 0)
            for k, chain in enumerate(selected):
                chain_buffer.append(info | dict(block=block, chain=int(chain), threshold=threshold,
                    total_contacts=int(totals[k]), partner_count=int(partners[k]),
                    contacts_per_partner=ratios[k]))
            frame_rows.append(info | dict(block=block, threshold=threshold, n_dense=len(selected),
                n_valid_contacts_per_partner=int(np.isfinite(ratios).sum()),
                total_contacts=float(totals.mean()) if len(selected) else np.nan,
                partner_count=float(partners.mean()) if len(selected) else np.nan,
                contacts_per_partner=float(np.nanmean(ratios)) if np.isfinite(ratios).any() else np.nan))
            frequency = np.bincount(partners, minlength=n) if len(selected) else None
            # Empty dense frames carry NaN fractions, never a zero distribution.
            for k in range(n):
                fraction = float(frequency[k] / len(selected)) if len(selected) else np.nan
                distribution_buffer.append(info | dict(block=block, threshold=threshold,
                                                        partner_count=k, fraction=fraction))
                key = (threshold, k, block)
                total, valid = distribution_sums.get(key, (0.0, 0))
                distribution_sums[key] = (total + fraction, valid + 1) if len(selected) else (total, valid)
        dense_counts = counts[np.ix_(selected, selected)]
        for threshold in parameters['network_thresholds']:
            efficiency, giant = graph_stats(dense_counts >= threshold)
            networks.append(info | dict(block=block, threshold=threshold, n_dense=len(selected),
                global_efficiency=efficiency, largest_component_fraction=giant))
        if (info['frame'] + 1) % parameters['chunk_frames'] == 0:
            if chain_buffer:
                io.append_csv(chain_buffer, output / 'chain_partners.csv.gz')
                chain_buffer.clear()
            io.append_csv(distribution_buffer, output / 'partner_distribution_frames.csv.gz')
            distribution_buffer.clear()
            flush_edges()
            print(f'  {simulation["id"]}: {info["frame"] + 1} intermolecular frames', flush=True)
    if not frame_rows:
        raise ValueError('No frames')
    if chain_buffer:
        io.append_csv(chain_buffer, output / 'chain_partners.csv.gz')
    else:
        # All-empty dense trajectories still have a readable chain-level schema.
        if not (output / 'chain_partners.csv.gz').exists():
            pd.DataFrame(columns=['simulation_id', 'frame', 'source_frame', 'time_ps', 'block', 'chain',
                'threshold', 'total_contacts', 'partner_count', 'contacts_per_partner']).to_csv(output / 'chain_partners.csv.gz', index=False)
    flush_edges()
    frame = pd.DataFrame(frame_rows)
    frame.to_csv(output / 'frame_means.csv', index=False)
    io.summarize(frame, ['total_contacts', 'partner_count', 'contacts_per_partner'],
                 ['simulation_id', 'threshold'], output)
    if distribution_buffer:
        io.append_csv(distribution_buffer, output / 'partner_distribution_frames.csv.gz')
    distribution = pd.DataFrame([
        dict(simulation_id=simulation['id'], threshold=t, partner_count=k, block=b,
             fraction=total / count if count else np.nan, valid_frames=count, total_fraction=total)
        for (t, k, b), (total, count) in distribution_sums.items()])
    distribution.drop(columns='total_fraction').to_csv(output / 'partner_distribution_blocks.csv', index=False)
    keys = ['simulation_id', 'threshold', 'partner_count']
    overall = distribution.groupby(keys)[['total_fraction', 'valid_frames']].sum()
    overall['mean'] = overall.total_fraction / overall.valid_frames.replace(0, np.nan)
    overall.rename(columns={'valid_frames': 'count'}).drop(columns='total_fraction').reset_index().to_csv(
        output / 'partner_distribution.csv', index=False)
    network = pd.DataFrame(networks)
    network.to_csv(output / 'network_frames.csv', index=False)
    destination = output / 'network'
    destination.mkdir()
    io.summarize(network, ['global_efficiency', 'largest_component_fraction'],
                 ['simulation_id', 'threshold'], destination)
