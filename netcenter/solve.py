"""High-level orchestration from a built network to solved map coordinates.

The mathematical routines operate only on arrays. This module is the small
bridge that chooses demand nodes, computes shortest paths, solves the three
objectives, and attaches x/y coordinates to the results.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from netcenter.center import (
    DEFAULT_MAX_CELLS,
    CenterResult,
    absolute_center,
    vertex_center,
    weighted_median,
)
from netcenter.distances import DEFAULT_TEMP_MB, distance_matrix

if TYPE_CHECKING:
    from netcenter.graph import Network


def _locate(net: Network, res: CenterResult) -> CenterResult:
    """Attach map coordinates to a solver result."""
    if res.node is not None:
        x, y = net.node_xy[res.node]
        res.xy = (float(x), float(y))
    elif res.edge is not None:
        if res.t is None:
            raise ValueError("edge result is missing its offset t")
        res.xy = net.interpolate(res.edge, res.t)
    return res


def _coalesce_demands(demand_nodes, weights, n_nodes: int):
    """Collapse repeated snapped nodes without changing either objective.

    Why this is exact:

    * For the minimax centre, repeating the same location does not change a
      maximum at all.
    * For the median, repetitions *do* matter, so their counts (or supplied
      weights) are summed into one weight for that unique node.

    The result can greatly reduce Dijkstra work when many centroids or points
    snap to the same road junction.
    """
    node_raw = np.asarray(demand_nodes)
    if not np.issubdtype(node_raw.dtype, np.integer):
        raise ValueError("demand_nodes must contain integer node indices")
    nodes = node_raw.astype(np.int64, copy=False).reshape(-1)
    if nodes.size == 0:
        raise ValueError("no demand points; every location was dropped or filtered out")
    if nodes.min() < 0 or nodes.max() >= n_nodes:
        raise ValueError("demand node index is outside the network")

    unique, inverse = np.unique(nodes, return_inverse=True)

    if weights is None:
        median_weights = np.bincount(inverse).astype(np.float64)
    else:
        w = np.asarray(weights, dtype=np.float64).reshape(-1)
        if len(w) != len(nodes):
            raise ValueError(
                f"got {len(w)} weights for {len(nodes)} demand points; "
                "they must correspond one to one"
            )
        if not np.isfinite(w).all():
            raise ValueError("weights contain NaN or infinity")
        if (w < 0).any():
            raise ValueError("weights must be non-negative")
        if not (w > 0).any():
            raise ValueError("at least one demand weight must be positive")
        median_weights = np.bincount(inverse, weights=w).astype(np.float64)

    # Avoid a needless matrix multiply when every unique node counts exactly
    # once. weighted_median's unweighted sum is slightly simpler and faster.
    if np.all(median_weights == 1.0):
        median_weights = None
    return unique, median_weights


def solve(
    net: Network,
    demand_nodes=None,
    weights=None,
    n_jobs: int = 1,
    backend: str = "loky",
    sweep_backend: str = "threading",
    dtype=np.float64,
    max_cells: int = DEFAULT_MAX_CELLS,
    max_temp_mb: float = DEFAULT_TEMP_MB,
    include_absolute: bool = True,
) -> dict[str, CenterResult]:
    """Compute the weighted 1-median, vertex centre, and absolute centre.

    ``demand_nodes`` defaults to every network node. That is mathematically
    valid but can be expensive; for applied facility-location work, snapped
    settlement/facility/population points are usually more meaningful.

    ``weights`` affects the median only. The centre implemented here is the
    unweighted minimax problem over the chosen demand locations.

    ``n_jobs=1`` is deliberately conservative.  If more workers are requested,
    ``backend`` controls shortest-path scheduling while ``sweep_backend``
    controls the absolute-centre sweep.  Separate settings are exposed because
    process and thread trade-offs differ between these two stages and across
    machines; benchmark on the target system before assuming either is faster.
    """
    if net.n_nodes == 0:
        raise ValueError("network has no nodes")

    if demand_nodes is None:
        demand_nodes = np.arange(net.n_nodes, dtype=np.int64)

    unique_nodes, median_weights = _coalesce_demands(demand_nodes, weights, net.n_nodes)

    D = distance_matrix(
        net.csr,
        unique_nodes,
        n_jobs=n_jobs,
        backend=backend,
        dtype=dtype,
        max_temp_mb=max_temp_mb,
    )

    # A disconnected demand/network pair appears as infinity. All objectives
    # become meaningless in that case, so fail before returning a plausible map.
    if not np.isfinite(D.max()):
        raise ValueError(
            "some demand points cannot reach some network nodes by road. "
            "Rebuild with keep_largest_component=True or use a connected study "
            "network."
        )

    out = {
        "median": _locate(net, weighted_median(D, median_weights)),
        "vertex_center": _locate(net, vertex_center(D)),
    }
    if include_absolute:
        out["absolute_center"] = _locate(
            net,
            absolute_center(
                D,
                net.edge_u,
                net.edge_w,
                net.edge_len,
                n_jobs=n_jobs,
                backend=sweep_backend,
                max_cells=max_cells,
            ),
        )
    return out
