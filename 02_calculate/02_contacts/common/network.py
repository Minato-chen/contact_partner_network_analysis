"""Contact-only residue networks; backbone bonds are deliberately absent."""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path


def chain_statistics(chain_ids, residue_i, residue_j, n_chains, n_residues):
    """Return [largest cluster fraction, efficiency] and normalized P(k) per chain.

    All residues are nodes, including isolates. Input edges are unique unordered
    residue pairs. Distances count contact edges, not spatial or sequence length.
    """
    statistics = np.empty((n_chains, 2), dtype=float)
    histograms = np.empty((n_chains, n_residues), dtype=float)
    for chain in range(n_chains):
        selected = chain_ids == chain
        i, j = residue_i[selected], residue_j[selected]
        degree = np.bincount(np.concatenate((i, j)), minlength=n_residues)
        histograms[chain] = np.bincount(degree, minlength=n_residues) / n_residues
        if not len(i):
            statistics[chain] = (1 / n_residues, 0.0 if n_residues > 1 else np.nan)
            continue
        graph = csr_matrix(
            (np.ones(2 * len(i)), (np.concatenate((i, j)), np.concatenate((j, i)))),
            shape=(n_residues, n_residues),
        )
        _, labels = connected_components(graph, directed=False)
        largest = np.bincount(labels).max() / n_residues
        distances = shortest_path(graph, directed=False, unweighted=True)
        inverse = np.zeros_like(distances)
        reachable = np.isfinite(distances) & (distances > 0)
        inverse[reachable] = 1 / distances[reachable]
        efficiency = inverse.sum() / (n_residues * (n_residues - 1))
        statistics[chain] = largest, efficiency
    return statistics, histograms
