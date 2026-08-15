"""Convert line geometry into the compact undirected network used by netcenter.

A road layer is not automatically a routing graph. The important distinction is
between *geometry* (lines that happen to cross on a map) and *topology* (places
where travel can actually move from one road to another).

Safe default
------------
``build_network(..., node=False)`` preserves source topology: it splits at
vertices already shared by distinct input lines, but does not invent junctions
at geometric crossings. This is the safer default for routing-quality data,
including data that distinguishes bridges, tunnels, and grade-separated roads.

Optional planar noding
----------------------
``node=True`` splits every two-dimensional line crossing and therefore assumes
that every drawn crossing is a valid junction. That can be useful for simple
street drawings, but it is wrong for overpasses and underpasses. The option is
explicit because silently inventing a turn is worse than failing loudly.

Other protections in this module handle empty geometries, closed rings,
near-identical endpoints, duplicate parallel edges, metric CRS units, and
stranded connected components.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import shapely
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from netcenter.topology import adjacency_from_edges

_LINESTRING = 1


@dataclass(slots=True)
class Network:
    """Road network represented as numeric arrays plus original edge geometry."""

    node_xy: np.ndarray
    edge_u: np.ndarray
    edge_w: np.ndarray
    edge_len: np.ndarray
    edge_geom: np.ndarray
    csr: csr_matrix
    crs: object = None

    @property
    def n_nodes(self) -> int:
        return len(self.node_xy)

    @property
    def n_edges(self) -> int:
        return len(self.edge_u)

    def interpolate(self, edge: int, t: float) -> tuple[float, float]:
        """Return map coordinates ``t`` metres from the u-end of one edge."""
        if edge < 0 or edge >= self.n_edges:
            raise IndexError("edge index is outside the network")
        if not np.isfinite(t):
            raise ValueError("t must be finite")
        if t < -1e-9 or t > self.edge_len[edge] + 1e-9:
            raise ValueError("t falls outside the selected edge")
        # Clamp tiny floating-point excursions at exactly 0 or L.
        t = float(np.clip(t, 0.0, self.edge_len[edge]))
        pt = shapely.line_interpolate_point(self.edge_geom[edge], t)
        return (float(shapely.get_x(pt)), float(shapely.get_y(pt)))


def _explode_lines(geoms: np.ndarray) -> np.ndarray:
    """Flatten multipart geometry into non-empty LineStrings only."""
    parts = shapely.get_parts(geoms)
    parts = parts[shapely.get_type_id(parts) == _LINESTRING]
    return parts[~shapely.is_empty(parts)]


def _node_and_merge(geoms: np.ndarray) -> np.ndarray:
    """Planarise all crossings, then merge degree-2 shape segments.

    This treats *every* 2-D crossing as connected. Call it only when that is a
    valid assumption for the input data.
    """
    merged = shapely.line_merge(shapely.union_all(geoms))
    return _explode_lines(np.asarray([merged], dtype=object))


def _split_at_shared_vertices(parts: np.ndarray, snap: float) -> np.ndarray:
    """Split lines at vertices that are shared by *different* input lines.

    This is the topology-preserving default for OSM-shaped linework.  A real
    junction is often encoded as one coordinate that belongs to two or more
    ways, but that coordinate need not be an endpoint of either way.  Splitting
    only where a side street ends on a through-road therefore recovers many
    T-junctions but still misses X-junctions where both ways continue.

    We instead identify coordinate locations used by at least two distinct
    LineStrings and split every participating line at that location when it is
    an interior vertex.  Crucially, this does *not* invent a coordinate at a
    mere geometric crossing.  A bridge crossing a road below shares no source
    vertex, so the two lines remain disconnected unless the caller explicitly
    asks for full planar noding with ``node=True``.

    Coordinates are compared on the same ``snap`` grid later used for endpoint
    identity.  This keeps the meaning of "same network location" consistent
    throughout graph construction.
    """
    if len(parts) < 2:
        return parts

    xy, line_id = shapely.get_coordinates(parts, return_index=True)
    if len(xy) == 0:
        return parts

    starts = np.searchsorted(line_id, np.arange(len(parts)), side="left")
    ends = np.searchsorted(line_id, np.arange(len(parts)), side="right") - 1

    # Quantise all source vertices once.  Locations shared by two distinct
    # lines are genuine topology supplied by the data; unlike planar noding, no
    # new intersection coordinate is created here.
    grid = _quantise_xy(xy, snap)
    _, location_id = np.unique(grid, axis=0, return_inverse=True)
    location_id = location_id.astype(np.int64, copy=False)

    # A location is shared when the minimum and maximum contributing line IDs
    # differ.  Sorting one integer vector is materially lighter than building a
    # second (location, line) matrix and running another 2-D unique operation.
    # Repeated visits by the *same* line therefore do not create a false junction.
    order = np.argsort(location_id, kind="stable")
    loc_sorted = location_id[order]
    line_sorted = line_id[order]
    first = np.r_[True, loc_sorted[1:] != loc_sorted[:-1]]
    starts_of_groups = np.flatnonzero(first)
    min_line = np.minimum.reduceat(line_sorted, starts_of_groups)
    max_line = np.maximum.reduceat(line_sorted, starts_of_groups)
    shared_ids = loc_sorted[starts_of_groups][min_line != max_line]
    shared_location = np.zeros(int(location_id.max()) + 1, dtype=bool)
    shared_location[shared_ids] = True

    interior = np.ones(len(xy), dtype=bool)
    interior[starts] = False
    interior[ends] = False
    cut_here = interior & shared_location[location_id]
    if not cut_here.any():
        return parts

    needs_cut = np.zeros(len(parts), dtype=bool)
    np.logical_or.at(needs_cut, line_id[cut_here], True)

    out: list[object] = []
    for p in range(len(parts)):
        if not needs_cut[p]:
            out.append(parts[p])
            continue

        lo, hi = starts[p], ends[p] + 1
        coords = xy[lo:hi]
        local_cuts = np.flatnonzero(cut_here[lo:hi])
        positions = [0, *local_cuts.tolist(), len(coords) - 1]

        # Consecutive repeated cut positions are harmless but produce zero-size
        # pieces.  Deduplicate positions before rebuilding the line.
        positions = np.unique(np.asarray(positions, dtype=np.int64))
        for a, b in zip(positions[:-1], positions[1:], strict=True):
            piece = coords[a : b + 1]
            if len(piece) >= 2:
                out.append(shapely.linestrings(piece))

    return np.asarray(out, dtype=object)


def _quantise_xy(xy: np.ndarray, snap: float) -> np.ndarray:
    """Map coordinates to an integer snap grid with an overflow guard."""
    xy = np.asarray(xy, dtype=np.float64)
    scaled = xy / float(snap)
    if not np.isfinite(scaled).all():
        raise ValueError("coordinates are too large relative to the snap grid spacing")

    limit = np.iinfo(np.int64).max - 1
    if np.abs(scaled).max(initial=0.0) > limit:
        raise ValueError(
            "snap is too small for the coordinate magnitude; choose a larger "
            "snap grid spacing"
        )
    return np.rint(scaled).astype(np.int64)


def _split_closed_rings(parts: np.ndarray, tolerance: float) -> np.ndarray:
    """Cut closed/nearly-closed lines before endpoint snapping deletes them.

    Real GIS rings are not always *bit-for-bit* closed: the final coordinate can
    differ from the first by floating-point dust. We therefore use the same
    endpoint tolerance that will later define node identity.
    """
    if len(parts) == 0:
        return parts
    p0, p1 = _endpoints(parts)
    nearly_closed = np.linalg.norm(p0 - p1, axis=1) <= tolerance
    is_ring = np.asarray(shapely.is_closed(parts)) | nearly_closed
    if not is_ring.any():
        return parts

    from shapely.ops import substring

    out = list(parts[~is_ring])
    for geom in parts[is_ring]:
        total = float(geom.length)
        if total <= 0:
            continue
        out.append(substring(geom, 0.0, total / 2.0))
        out.append(substring(geom, total / 2.0, total))
    return np.asarray(out, dtype=object)


def _endpoints(parts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised extraction of the first and last coordinate of each line."""
    xy, idx = shapely.get_coordinates(parts, return_index=True)
    starts = np.searchsorted(idx, np.arange(len(parts)), side="left")
    ends = np.searchsorted(idx, np.arange(len(parts)), side="right") - 1
    return xy[starts], xy[ends]


def _is_metre_crs(crs) -> bool:
    """True when both horizontal axes use metres as their linear unit."""
    from pyproj import CRS

    crs = CRS.from_user_input(crs)
    if not crs.is_projected:
        return False
    axes = crs.axis_info[:2]
    if len(axes) < 2:
        return False
    return all(np.isclose(axis.unit_conversion_factor, 1.0) for axis in axes)


def _working_crs(gdf, target_crs=None):
    """Choose a projected metre-based CRS and reject ambiguous unit mistakes."""
    from pyproj import CRS

    source_crs = CRS.from_user_input(gdf.crs)

    if target_crs is not None:
        chosen = CRS.from_user_input(target_crs)
        if not chosen.is_projected:
            raise ValueError("target_crs must be projected, not latitude/longitude")
        if not _is_metre_crs(chosen):
            raise ValueError(
                "target_crs must use metres because all netcenter distances and "
                "snap grid spacing is defined in metres"
            )
        return chosen

    if source_crs.is_projected and _is_metre_crs(source_crs):
        # Do not reproject good projected input: changing CRS needlessly changes
        # every measured edge length a little.
        return source_crs

    # Geographic input, or projected input in feet/another unit: choose a local
    # UTM CRS so lengths and tolerances have the documented metre interpretation.
    chosen = gdf.estimate_utm_crs()
    if chosen is None:
        raise ValueError(
            "could not infer a local metre-based CRS; pass target_crs explicitly"
        )

    if source_crs.is_geographic:
        minx, miny, maxx, maxy = gdf.total_bounds
        if maxx - minx > 6.0:
            warnings.warn(
                "input spans more than 6 degrees of longitude; one inferred UTM "
                "zone may distort a very large study area. Pass an explicit "
                "metre-based projected CRS appropriate to the full extent.",
                UserWarning,
                stacklevel=3,
            )
    else:
        warnings.warn(
            "input CRS is projected but not metre-based; reprojecting to an "
            "estimated local UTM CRS so reported distances remain metres.",
            UserWarning,
            stacklevel=3,
        )
    return CRS.from_user_input(chosen)


def _compact_nodes(
    node_xy: np.ndarray, u: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Drop nodes no surviving edge uses and remap endpoints densely.

    Geometry filtering can remove a sliver or snapped self-loop after node IDs
    have already been assigned.  Leaving those now-orphaned nodes in the CSR
    creates isolated components that can later make an otherwise connected
    network appear unreachable.
    """
    used = np.unique(np.concatenate([u, w]))
    if len(used) == len(node_xy):
        return node_xy, u, w

    remap = np.full(len(node_xy), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    return node_xy[used], remap[u], remap[w]


def _largest_component_by_length(labels, u, w, lengths) -> int:
    """Choose the connected component containing the most road length."""
    n_components = int(labels.max()) + 1
    totals = np.zeros(n_components, dtype=np.float64)
    # Every retained edge has endpoints in one component. Count its length once.
    np.add.at(totals, labels[u], lengths)
    return int(np.argmax(totals))


def build_network(
    source,
    layer: str | None = None,
    target_crs=None,
    snap: float = 0.1,
    node: bool = False,
    split_shared_vertices: bool = True,
    min_length: float = 1e-9,
    keep_largest_component: bool = True,
) -> Network:
    """Build a :class:`Network` from a vector file or GeoDataFrame.

    Parameters
    ----------
    source
        Shapefile/GeoPackage/etc. path or a GeoDataFrame containing linework.
    layer
        Optional layer name for multi-layer files.
    target_crs
        Optional projected CRS whose horizontal units are metres. Geographic
        input is otherwise moved to an inferred local UTM CRS automatically.
    snap
        Node-identity grid spacing in metres. Endpoints and shared source
        vertices are quantised to this grid; arbitrary points are never
        projected onto roads. Closed-ring detection also uses this value as a
        small Euclidean closure threshold.
    node
        If true, split every geometric line crossing. Use only for a planar
        network where every 2-D crossing is a genuine junction. This also
        connects bridges to whatever passes beneath them.
    split_shared_vertices
        If true (the default), split lines at vertices already shared by two or
        more distinct input LineStrings. This recovers OSM-style T- and
        X-junctions without inventing a junction at a bridge/underpass crossing.
    min_length
        Segments at or below this length are discarded as slivers.
    keep_largest_component
        If true, retain the connected component with the greatest road length.
    """
    import geopandas as gpd

    if not np.isfinite(snap) or snap <= 0:
        raise ValueError("snap must be a finite positive distance in metres")
    if not np.isfinite(min_length) or min_length < 0:
        raise ValueError("min_length must be a finite non-negative distance")

    gdf = (
        source
        if isinstance(source, gpd.GeoDataFrame)
        else gpd.read_file(source, layer=layer)
    )

    # Remove missing/empty geometry before endpoint indexing; an empty geometry
    # contributes no coordinates and can otherwise shift subsequent endpoints.
    gdf = gdf[~gdf.geometry.is_empty & ~gdf.geometry.isna()]
    gdf = gdf[gdf.geom_type.isin(["LineString", "MultiLineString"])]
    if len(gdf) == 0:
        raise ValueError("no usable line geometries found in source")
    if gdf.crs is None:
        raise ValueError(
            "source has no coordinate reference system; set one before building"
        )

    crs = _working_crs(gdf, target_crs)
    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)

    geoms = np.asarray(gdf.geometry.values, dtype=object)
    if node:
        warnings.warn(
            "planar noding is enabled: every geometric crossing will become a "
            "junction, including bridges/underpasses unless the input has already "
            "been separated geometrically.",
            UserWarning,
            stacklevel=2,
        )
        parts = _node_and_merge(geoms)
    else:
        parts = _explode_lines(geoms)
        if split_shared_vertices:
            parts = _split_at_shared_vertices(parts, snap)

    parts = _split_closed_rings(parts, tolerance=snap)
    if len(parts) == 0:
        raise ValueError("no line segments survived preprocessing")

    p0, p1 = _endpoints(parts)
    lengths = np.asarray(shapely.length(parts), dtype=np.float64)

    # Quantise endpoints to integer bins. This avoids unreliable direct equality
    # tests on floating-point map coordinates while keeping the original drawn
    # coordinate as the representative node location.
    stacked = np.vstack([p0, p1])
    quantised = _quantise_xy(stacked, snap)
    _, first_seen, inverse = np.unique(
        quantised, axis=0, return_index=True, return_inverse=True
    )
    inverse = inverse.ravel()
    node_xy = stacked[first_seen]

    m = len(parts)
    u, w = inverse[:m], inverse[m:]
    usable = (u != w) & np.isfinite(lengths) & (lengths > min_length)
    u, w = u[usable], w[usable]
    lengths, parts = lengths[usable], parts[usable]
    if len(u) == 0:
        raise ValueError("no usable segments after removing loops and slivers")

    # Node IDs were created before sliver/self-loop filtering.  Compact them now
    # so filtered geometry cannot leave isolated orphan nodes in the graph.
    node_xy, u, w = _compact_nodes(node_xy, u, w)
    csr = adjacency_from_edges(u, w, lengths, len(node_xy))

    if keep_largest_component:
        n_components, labels = connected_components(csr, directed=False)
        if n_components > 1:
            keep_label = _largest_component_by_length(labels, u, w, lengths)
            node_keep = labels == keep_label
            remap = np.full(len(node_xy), -1, dtype=np.int64)
            remap[node_keep] = np.arange(int(node_keep.sum()))
            edge_keep = node_keep[u] & node_keep[w]

            kept_length = float(lengths[edge_keep].sum())
            total_length = float(lengths.sum())
            warnings.warn(
                f"network split into {n_components} disconnected pieces; keeping "
                f"the component with the most road length "
                f"({kept_length:,.1f} of {total_length:,.1f} m).",
                UserWarning,
                stacklevel=2,
            )

            u, w = remap[u[edge_keep]], remap[w[edge_keep]]
            lengths, parts = lengths[edge_keep], parts[edge_keep]
            node_xy = node_xy[node_keep]
            csr = adjacency_from_edges(u, w, lengths, len(node_xy))

    return Network(
        node_xy=np.asarray(node_xy, dtype=np.float64),
        edge_u=np.asarray(u, dtype=np.int64),
        edge_w=np.asarray(w, dtype=np.int64),
        edge_len=np.asarray(lengths, dtype=np.float64),
        edge_geom=np.asarray(parts, dtype=object),
        csr=csr,
        crs=crs,
    )


def snap_points(net: Network, xy: np.ndarray, max_dist: float | None = None):
    """Attach point demand to the nearest *network node*.

    This is intentionally node snapping, not nearest-edge snapping. The returned
    distance tells the caller how far each demand point was moved, so a sensible
    ``max_dist`` can prevent long simplified edges from causing hidden shifts.
    """
    xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
    if not np.isfinite(xy).all():
        raise ValueError("demand coordinates contain NaN or infinity")
    if max_dist is not None and (not np.isfinite(max_dist) or max_dist < 0):
        raise ValueError("max_dist must be a finite non-negative distance")
    if net.n_nodes == 0:
        raise ValueError("network has no nodes")

    dist, idx = cKDTree(net.node_xy).query(xy, k=1)
    idx = idx.astype(np.int64)
    if max_dist is not None:
        idx = np.where(dist <= max_dist, idx, -1)
    return idx, dist
