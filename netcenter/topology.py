"""Turn segment arrays into the sparse graph used by SciPy shortest paths.

This module intentionally has no GIS dependencies. Keeping topology separate
from file I/O makes the core easy to test and prevents subtle differences
between production code and synthetic test networks.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


def adjacency_from_edges(u, w, length, n_nodes: int) -> csr_matrix:
    """Build an undirected sparse adjacency matrix from road segments.

    Parallel edges are collapsed by keeping the *shortest* segment between a
    pair of nodes. This is essential: converting a COO matrix with duplicate
    coordinates directly to CSR makes SciPy add the duplicate values, so two
    500 m roads can silently become a 1,000 m connection.

    Parameters
    ----------
    u, w
        Integer arrays containing the node at each end of every segment.
    length
        Positive finite segment lengths, in metres.
    n_nodes
        Total number of nodes in the network.
    """
    if (
        isinstance(n_nodes, (bool, np.bool_))
        or not isinstance(n_nodes, (int, np.integer))
        or n_nodes < 0
    ):
        raise ValueError("n_nodes must be a non-negative integer")
    n_nodes = int(n_nodes)

    u_raw = np.asarray(u)
    w_raw = np.asarray(w)
    if not np.issubdtype(u_raw.dtype, np.integer) or not np.issubdtype(
        w_raw.dtype, np.integer
    ):
        raise ValueError("u and w must contain integer node indices")
    u = u_raw.astype(np.int64, copy=False).reshape(-1)
    w = w_raw.astype(np.int64, copy=False).reshape(-1)
    length = np.asarray(length, dtype=np.float64).reshape(-1)

    if not (len(u) == len(w) == len(length)):
        raise ValueError("u, w and length must all be the same length")
    if len(u) == 0:
        return csr_matrix((n_nodes, n_nodes), dtype=np.float64)
    if n_nodes == 0:
        raise ValueError("edges were supplied but n_nodes is zero")
    if u.min() < 0 or w.min() < 0 or u.max() >= n_nodes or w.max() >= n_nodes:
        raise ValueError("edge endpoints fall outside the range of node indices")
    if not np.isfinite(length).all():
        raise ValueError("edge lengths must be finite")
    if (length <= 0).any():
        raise ValueError("edge lengths must be strictly positive")

    # Canonicalise each undirected pair so 5->9 and 9->5 are the same key.
    lo = np.minimum(u, w)
    hi = np.maximum(u, w)

    # Self-loops never shorten a route. Closed road rings are handled earlier
    # in graph.py by splitting them into ordinary segments before this point.
    keep = lo != hi
    lo, hi, length = lo[keep], hi[keep], length[keep]
    if len(lo) == 0:
        return csr_matrix((n_nodes, n_nodes), dtype=np.float64)

    # Sort primarily by node pair and secondarily by length. The first entry in
    # each run is therefore the shortest parallel segment for that pair.
    key = lo * np.int64(n_nodes) + hi
    order = np.lexsort((length, key))
    key_sorted = key[order]
    first = np.ones(len(key_sorted), dtype=bool)
    first[1:] = key_sorted[1:] != key_sorted[:-1]

    lo_u = lo[order][first]
    hi_u = hi[order][first]
    len_u = length[order][first]

    rows = np.concatenate([lo_u, hi_u])
    cols = np.concatenate([hi_u, lo_u])
    data = np.concatenate([len_u, len_u])
    return coo_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes)).tocsr()
