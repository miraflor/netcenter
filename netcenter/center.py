"""The two location problems, and the algorithms that solve them.

READING THIS FILE WITHOUT A PROGRAMMING BACKGROUND
--------------------------------------------------
Three ideas cover almost everything below.

*Arrays.* An array is a grid of numbers. ``D`` is our main one: row ``v``,
column ``j`` holds the road distance from demand point ``v`` to intersection
``j``. So ``D`` has one row per demand and one column per intersection.

*Vectorisation.* Rather than looping over rows one at a time, we ask numpy to
do the same arithmetic to every row simultaneously. ``a + b`` where both are
grids adds them cell by cell, in optimised C code. This is why the code below
looks like algebra rather than like step-by-step instructions -- each statement
is performing thousands of small calculations at once.

*Axis.* Operations take an ``axis`` telling them which direction to work in.
``axis=0`` collapses down the columns, ``axis=1`` collapses across the rows.
``D.max(axis=0)`` gives, for each intersection, the distance to its FARTHEST
demand point.

WHAT THE TWO PROBLEMS ARE
-------------------------
Both ask "where is the middle of this road network?" and they disagree.

The MEDIAN minimises the total distance everyone travels: the sum of
``weight * distance``. It is pulled towards where the people are. Hakimi proved
in 1964 that a best median always sits exactly on an intersection, so finding
it is just "score every intersection, keep the best" -- one pass over ``D``.

The CENTRE minimises the distance for the worst-off person: the largest
distance to anyone. It is pulled towards the geographic extremes and is decided
by the two or three most remote places. Crucially Hakimi's result does NOT hold
here: the best point can sit partway along a road, between two intersections.
Finding it needs the sweep described next.

HOW THE CENTRE SWEEP WORKS
--------------------------
Walk along one road segment from intersection ``u`` to intersection ``w``, and
call your position ``t`` metres from ``u``. The segment is ``L`` metres long.

To reach some demand point ``v`` you either walk back to ``u`` and travel from
there, or forward to ``w`` and travel from there. So your distance is

    d(t, v) = min(a_v + t,  b_v + L - t)

where ``a_v`` is the road distance from ``u`` to ``v`` and ``b_v`` from ``w``
to ``v``. As ``t`` grows this rises at first (going via ``u`` is still the
better route) and later falls (going via ``w`` takes over). Plotted against
``t`` it is a tent shape.

Your eccentricity at ``t`` is the worst of these over all demands -- the
highest tent above you. With many tents that upper outline is a jagged ridge of
peaks and valleys, and we want its lowest valley.

The trick that makes this fast: tent ``v`` is still rising exactly while
``t <= t*_v = (b_v - a_v + L) / 2``. Sort the demands by that switchover point.
Between two consecutive switchovers the set of still-rising tents is a fixed
block at one end of the sorted list and the set of already-falling tents a
fixed block at the other, so across that whole stretch the ridge simplifies to

    ecc(t) = max(A + t,  B + L - t)

where ``A`` is the largest ``a`` among the risers and ``B`` the largest ``b``
among the fallers -- two numbers that do not change within the stretch. That is
a simple V shape, and the bottom of a V is found by setting its two lines
equal. So each stretch is solved by arithmetic, with no searching or sampling.

Sorting costs ``k log k`` per segment for ``k`` demands. Everything else is
running maximums, which numpy computes in a single pass.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# How many array cells to process at once. Peak memory in the sweep is roughly
# 15 arrays of this size in 8-byte floats, so 1,000,000 is about 120 MB per
# worker. Lower it on a memory-constrained machine; it never changes results.
DEFAULT_MAX_CELLS = 1_000_000

# Below this much work, worker startup is usually not worth attempting.
MIN_PARALLEL_WORK = 5_000_000


@dataclass
class CenterResult:
    """One solved location and its objective value.

    A vertex result has ``node`` set. A continuous-network result has ``edge``
    and ``t`` set, where ``t`` is distance from that edge's ``u`` endpoint.
    ``xy`` is attached later by :mod:`netcenter.solve` because the mathematical
    core intentionally knows nothing about map coordinates.
    """

    kind: str
    objective: float
    node: int | None = None
    edge: int | None = None
    t: float | None = None
    xy: tuple[float, float] | None = None


def _validate_D(D: np.ndarray) -> np.ndarray:
    """Validate the demand-by-node distance matrix used by every solver."""
    D = np.asarray(D)
    if D.ndim != 2:
        raise ValueError("D must be a 2-D array shaped (demands, nodes)")
    if D.shape[0] == 0:
        raise ValueError("D contains no demand points")
    if D.shape[1] == 0:
        raise ValueError("D contains no network nodes")
    if not np.issubdtype(D.dtype, np.number):
        raise ValueError("D must contain numeric distances")
    if not np.isfinite(D).all():
        raise ValueError("D contains NaN or infinite distances")
    if (D < 0).any():
        raise ValueError("D contains negative distances")
    return D


def _validate_edges(D, edge_u, edge_w, edge_len):
    """Normalise and validate the segment arrays used by the centre sweep."""
    edge_u_raw = np.asarray(edge_u)
    edge_w_raw = np.asarray(edge_w)
    if not np.issubdtype(edge_u_raw.dtype, np.integer) or not np.issubdtype(
        edge_w_raw.dtype, np.integer
    ):
        raise ValueError("edge_u and edge_w must contain integer node indices")
    edge_u = edge_u_raw.astype(np.int64, copy=False).reshape(-1)
    edge_w = edge_w_raw.astype(np.int64, copy=False).reshape(-1)
    edge_len = np.asarray(edge_len, dtype=np.float64).reshape(-1)
    if not (len(edge_u) == len(edge_w) == len(edge_len)):
        raise ValueError("edge_u, edge_w and edge_len must be the same length")
    if len(edge_u):
        if edge_u.min() < 0 or edge_w.min() < 0:
            raise ValueError("edge endpoints must be non-negative node indices")
        if edge_u.max() >= D.shape[1] or edge_w.max() >= D.shape[1]:
            raise ValueError("edge endpoint index is outside D's node columns")
        if not np.isfinite(edge_len).all():
            raise ValueError("edge lengths must be finite")
        if (edge_len <= 0).any():
            raise ValueError("edge lengths must be strictly positive")
    return edge_u, edge_w, edge_len


def vertex_eccentricity(D: np.ndarray) -> np.ndarray:
    """For every node, return its distance to the farthest demand point."""
    D = _validate_D(D)
    return D.max(axis=0)


def vertex_center(D: np.ndarray) -> CenterResult:
    """Best *node* under the minimax objective (the Jordan centre)."""
    ecc = vertex_eccentricity(D)
    node = int(np.argmin(ecc))
    return CenterResult("vertex_center", float(ecc[node]), node=node)


def weighted_median(D: np.ndarray, weights=None) -> CenterResult:
    """Best node under total weighted distance.

    Hakimi's vertex-optimality result means an optimum for vertex demand can be
    chosen at a network vertex, so the continuous edge interiors do not need a
    second search for this objective.
    """
    D = _validate_D(D)
    if weights is None:
        obj = D.sum(axis=0)
    else:
        w = np.asarray(weights, dtype=np.float64).reshape(-1)
        if w.shape[0] != D.shape[0]:
            raise ValueError(
                f"got {w.shape[0]} weights for {D.shape[0]} demand points; "
                "they must correspond one to one"
            )
        if not np.isfinite(w).all():
            raise ValueError("weights contain NaN or infinity")
        if (w < 0).any():
            raise ValueError("weights must be non-negative")
        if not (w > 0).any():
            raise ValueError("at least one demand weight must be positive")
        obj = w @ D
    node = int(np.argmin(obj))
    return CenterResult("median", float(obj[node]), node=node)


def _validate_jobs(n_jobs: int) -> int:
    """Require a genuine positive integer worker count."""
    if (
        isinstance(n_jobs, (bool, np.bool_))
        or not isinstance(n_jobs, (int, np.integer))
        or n_jobs < 1
    ):
        raise ValueError("n_jobs must be a positive integer")
    return int(n_jobs)


def _lower_bounds(D, edge_u, edge_w, max_cells) -> np.ndarray:
    """Cheap segment-level lower bound used to prune the exact sweep."""
    if (
        isinstance(max_cells, (bool, np.bool_))
        or not isinstance(max_cells, (int, np.integer))
        or max_cells <= 0
    ):
        raise ValueError("max_cells must be a positive integer")
    m = len(edge_u)
    step = max(1, int(max_cells // max(1, D.shape[0])))
    out = np.empty(m, dtype=np.float64)
    for s in range(0, m, step):
        e = min(s + step, m)
        out[s:e] = np.minimum(D[:, edge_u[s:e]], D[:, edge_w[s:e]]).max(axis=0)
    return out


def _sweep(D, edge_u, edge_w, edge_len, max_cells):
    """Solve the exact minimax point on every supplied segment.

    Arrays inside each block have shape ``(segments, demands)``. The block size
    is derived from ``max_cells`` so the temporary arrays have a predictable
    upper bound independent of the total number of road segments.
    """
    if (
        isinstance(max_cells, (bool, np.bool_))
        or not isinstance(max_cells, (int, np.integer))
        or max_cells <= 0
    ):
        raise ValueError("max_cells must be a positive integer")
    m = len(edge_u)
    k = D.shape[0]
    step = max(1, int(max_cells // max(1, k)))
    best_val = np.empty(m, dtype=np.float64)
    best_t = np.empty(m, dtype=np.float64)

    for s in range(0, m, step):
        e = min(s + step, m)

        # Distances from both endpoints to every demand. Making these arrays
        # contiguous pays for itself because the next operation sorts each row.
        a = np.ascontiguousarray(D[:, edge_u[s:e]].T, dtype=np.float64)
        b = np.ascontiguousarray(D[:, edge_w[s:e]].T, dtype=np.float64)
        L = edge_len[s:e][:, None]
        n_e = a.shape[0]

        # A demand switches from the u-route to the w-route at t*. A valid
        # shortest-path matrix should already put t* in [0,L]; clipping is a
        # defensive guard against tiny numerical violations.
        tstar = np.clip(0.5 * (b - a + L), 0.0, L)
        order = np.argsort(tstar, axis=1, kind="stable")
        ts = np.take_along_axis(tstar, order, axis=1)
        a_s = np.take_along_axis(a, order, axis=1)
        b_s = np.take_along_axis(b, order, axis=1)
        del a, b, order, tstar

        # k switch points create k+1 intervals. On each interval the upper
        # envelope reduces to max(A+t, B+L-t), a V whose minimum is analytic.
        neg = np.full((n_e, 1), -np.inf)
        lo = np.concatenate([np.zeros((n_e, 1)), ts], axis=1)
        hi = np.concatenate([ts, L], axis=1)
        A = np.concatenate(
            [np.maximum.accumulate(a_s[:, ::-1], axis=1)[:, ::-1], neg], axis=1
        )
        B = np.concatenate([neg, np.maximum.accumulate(b_s, axis=1)], axis=1)
        del a_s, b_s, ts

        tc = np.clip(0.5 * (B + L - A), lo, hi)
        val = np.maximum(A + tc, B + L - tc)
        del A, B, lo, hi

        j = np.argmin(val, axis=1)
        rows = np.arange(n_e)
        best_val[s:e] = val[rows, j]
        best_t[s:e] = tc[rows, j]

    return best_val, best_t


def absolute_center(
    D: np.ndarray,
    edge_u: np.ndarray,
    edge_w: np.ndarray,
    edge_len: np.ndarray,
    n_jobs: int = 1,
    backend: str = "threading",
    max_cells: int = DEFAULT_MAX_CELLS,
    min_parallel_work: int = MIN_PARALLEL_WORK,
    tol: float = 1e-9,
) -> CenterResult:
    """Return the exact absolute 1-centre for vertex demand.

    The best vertex gives an incumbent radius. A mathematically safe lower
    bound then discards segments that cannot improve it, and only the survivors
    are swept exactly. This means the expensive sorting is normally performed
    on a tiny fraction of the road network.
    """
    D = _validate_D(D)
    edge_u, edge_w, edge_len = _validate_edges(D, edge_u, edge_w, edge_len)
    if (
        isinstance(max_cells, (bool, np.bool_))
        or not isinstance(max_cells, (int, np.integer))
        or max_cells <= 0
    ):
        raise ValueError("max_cells must be a positive integer")
    if tol < 0 or not np.isfinite(tol):
        raise ValueError("tol must be a finite non-negative number")
    if not isinstance(min_parallel_work, (int, np.integer)) or min_parallel_work < 0:
        raise ValueError("min_parallel_work must be a non-negative integer")
    n_jobs = _validate_jobs(n_jobs)
    if backend not in {"threading", "loky"}:
        raise ValueError("backend must be 'threading' or 'loky'")

    # D may deliberately be stored as float32. The centre is then exact for the
    # stored distances, not for the unrounded float64 values SciPy originally
    # produced. More importantly, pruning must never discard an edge merely
    # because rounding made its lower bound a hair too large. Add a scale-aware
    # numerical margin and keep borderline edges rather than pruning them.
    if np.issubdtype(D.dtype, np.floating):
        roundoff = 8.0 * np.finfo(D.dtype).eps * max(1.0, float(D.max()))
    else:
        roundoff = 0.0
    effective_tol = max(float(tol), roundoff)

    best = vertex_center(D)
    if len(edge_u) == 0:
        return best

    lb = _lower_bounds(D, edge_u, edge_w, max_cells)
    keep = np.flatnonzero(lb < best.objective + effective_tol)
    if keep.size == 0:
        return best

    too_small = keep.size * D.shape[0] < min_parallel_work
    if n_jobs == 1 or keep.size < 2 * n_jobs or too_small:
        vals, offs = _sweep(D, edge_u[keep], edge_w[keep], edge_len[keep], max_cells)
    else:
        from joblib import Parallel, delayed

        blocks = [blk for blk in np.array_split(keep, n_jobs) if blk.size]
        parts = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(_sweep)(D, edge_u[b], edge_w[b], edge_len[b], max_cells) for b in blocks
        )
        vals = np.concatenate([p[0] for p in parts])
        offs = np.concatenate([p[1] for p in parts])

    i = int(np.argmin(vals))
    if vals[i] >= best.objective - effective_tol:
        return best
    return CenterResult(
        "absolute_center",
        float(vals[i]),
        edge=int(keep[i]),
        t=float(offs[i]),
    )
