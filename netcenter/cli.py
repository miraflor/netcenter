"""Command-line interface for users who do not want to write Python code."""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from netcenter import __version__
from netcenter.distances import DEFAULT_TEMP_MB, estimate_distance_matrix_mb
from netcenter.solve import DEFAULT_MAX_CELLS, solve


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="netcenter",
        description=(
            "Locate the weighted 1-median, vertex centre, and exact absolute "
            "1-centre of an undirected road network."
        ),
    )
    p.add_argument("roads", help="road line file (Shapefile, GeoPackage, etc.)")
    p.add_argument("--version", action="version", version=f"netcenter {__version__}")
    p.add_argument("--layer", default=None, help="road layer name if needed")
    p.add_argument(
        "--crs",
        default=None,
        help="metre-based projected CRS; otherwise geographic/non-metre input "
        "is moved to an inferred local UTM CRS",
    )
    p.add_argument(
        "--snap",
        type=float,
        default=0.1,
        help="node-identity grid spacing in metres (default: 0.1)",
    )
    p.add_argument(
        "--node-crossings",
        action="store_true",
        help="PLANAR networks only: turn every geometric crossing into a junction; "
        "do not use on bridges/underpasses",
    )
    # Compatibility with 0.1 commands. Since 0.2 already defaults to no planar
    # noding, this flag is a harmless no-op rather than an abrupt CLI break.
    p.add_argument("--no-node", action="store_true", help=argparse.SUPPRESS)
    p.add_argument(
        "--no-shared-vertex-noding",
        action="store_true",
        help=(
            "skip safe splitting at vertices shared by distinct input lines; "
            "use only when the source is already segmented at every real junction"
        ),
    )

    demand = p.add_argument_group("demand")
    demand.add_argument(
        "--demand",
        default=None,
        help="point/polygon demand layer; default is every network node",
    )
    demand.add_argument("--demand-layer", default=None)
    demand.add_argument(
        "--weight-field", default=None, help="numeric field weighting the median"
    )
    demand.add_argument(
        "--max-snap",
        type=float,
        default=None,
        help="drop demand points farther than this from the nearest NETWORK NODE",
    )

    perf = p.add_argument_group("performance")
    perf.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="maximum workers; conservative default is 1",
    )
    perf.add_argument(
        "--backend",
        default="loky",
        choices=["threading", "loky"],
        help=(
            "shortest-path backend when --jobs > 1: loky uses processes; "
            "threading shares memory (default: loky)"
        ),
    )
    perf.add_argument(
        "--float32",
        action="store_true",
        help="store the distance matrix in float32 to roughly halve its size",
    )
    perf.add_argument(
        "--max-temp-mb",
        type=float,
        default=DEFAULT_TEMP_MB,
        help=(
            "target max float64 Dijkstra block size per worker in MiB "
            f"(default: {DEFAULT_TEMP_MB})"
        ),
    )
    perf.add_argument(
        "--max-cells",
        type=int,
        default=DEFAULT_MAX_CELLS,
        help="maximum segment-demand cells per centre-sweep block",
    )

    p.add_argument(
        "--skip-absolute",
        action="store_true",
        help="solve only junction-restricted median and vertex centre",
    )
    p.add_argument("--out", default=None, help="optional point output file")
    p.add_argument("--quiet", action="store_true")
    return p


def _load_demand(args, net, log):
    """Read demand, convert polygons to inner points, then snap to network nodes."""
    import geopandas as gpd

    pts = gpd.read_file(args.demand, layer=args.demand_layer)
    pts = pts[~pts.geometry.isna() & ~pts.geometry.is_empty].copy()
    if len(pts) == 0:
        raise ValueError(f"demand layer {args.demand!r} contains no usable geometry")

    if not (pts.geom_type == "Point").all():
        log("note: non-point demand detected; using each geometry's representative point")
        pts = pts.set_geometry(pts.geometry.representative_point())

    if pts.crs is None:
        raise ValueError("demand layer has no coordinate reference system")
    pts = pts.to_crs(net.crs)

    if args.weight_field and args.weight_field not in pts.columns:
        fields = ", ".join(c for c in pts.columns if c != "geometry")
        raise ValueError(
            f"--weight-field {args.weight_field!r} is not in the demand layer; "
            f"available fields: {fields}"
        )

    from netcenter.graph import snap_points

    xy = np.column_stack([pts.geometry.x.to_numpy(), pts.geometry.y.to_numpy()])
    idx, dist = snap_points(net, xy, max_dist=args.max_snap)
    ok = idx >= 0
    if not ok.any():
        raise ValueError(
            "every demand point was beyond --max-snap; check study area and CRS"
        )
    if not ok.all():
        log(f"warning: dropped {int((~ok).sum()):,} demand points beyond --max-snap")

    weights = None
    if args.weight_field:
        try:
            weights = pts[args.weight_field].to_numpy(dtype=float)[ok]
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"--weight-field {args.weight_field!r} must be numeric"
            ) from exc

    farthest = float(dist[ok].max())
    unique = int(np.unique(idx[ok]).size)
    log(
        f"demand: {int(ok.sum()):,} points -> {unique:,} unique network nodes; "
        f"farthest snap {farthest:,.1f} m"
    )
    if args.max_snap is None and farthest > 1000:
        log(
            "warning: at least one demand point moved more than 1 km to a network "
            "node; consider --max-snap or a less aggressively simplified network"
        )
    return idx[ok], weights


def _write(path, results, crs, log):
    """Write solved locations as a geospatial point file."""
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import Point

    gdf = gpd.GeoDataFrame(
        {
            "kind": list(results),
            "objective": [r.objective for r in results.values()],
            "node": pd.array([r.node for r in results.values()], dtype="Int64"),
            "edge": pd.array([r.edge for r in results.values()], dtype="Int64"),
            "t_m": [r.t for r in results.values()],
        },
        geometry=[Point(*r.xy) for r in results.values()],
        crs=crs,
    )
    gdf.to_file(path)
    log(f"wrote {path}")


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    log = (lambda *a: None) if args.quiet else (lambda *a: print(*a, file=sys.stderr))

    try:
        from netcenter.graph import build_network

        t0 = time.perf_counter()
        net = build_network(
            args.roads,
            layer=args.layer,
            target_crs=args.crs,
            snap=args.snap,
            node=bool(args.node_crossings and not args.no_node),
            split_shared_vertices=not args.no_shared_vertex_noding,
        )
        log(
            f"graph: {net.n_nodes:,} nodes, {net.n_edges:,} segments, "
            f"crs={net.crs.to_string() if net.crs else '?'} "
            f"[{time.perf_counter() - t0:.1f}s]"
        )

        demand_nodes, weights = (None, None)
        if args.demand:
            demand_nodes, weights = _load_demand(args, net, log)

        k = net.n_nodes if demand_nodes is None else int(np.unique(demand_nodes).size)
        dtype = np.float32 if args.float32 else np.float64
        matrix_mb = estimate_distance_matrix_mb(k, net.n_nodes, dtype)
        log(f"distance matrix: about {matrix_mb:,.1f} MiB stored as {np.dtype(dtype).name}")

        t1 = time.perf_counter()
        results = solve(
            net,
            demand_nodes=demand_nodes,
            weights=weights,
            n_jobs=args.jobs,
            backend=args.backend,
            dtype=dtype,
            max_cells=args.max_cells,
            max_temp_mb=args.max_temp_mb,
            include_absolute=not args.skip_absolute,
        )
        log(f"solved in {time.perf_counter() - t1:.1f}s")
    except (ImportError, ValueError, OSError) as exc:
        raise SystemExit(f"error: {exc}") from exc

    n_demand = len(demand_nodes) if demand_nodes is not None else net.n_nodes
    total_weight = float(np.sum(weights)) if weights is not None else float(n_demand)
    for name, r in results.items():
        where = f"node {r.node}" if r.node is not None else f"edge {r.edge} @ {r.t:.1f} m"
        extra = f"  mean {r.objective / total_weight:,.1f} m" if name == "median" else ""
        print(
            f"{name:16s}  {r.objective:14,.3f}  {where:24s}  "
            f"({r.xy[0]:.3f}, {r.xy[1]:.3f}){extra}"
        )

    if args.out:
        _write(args.out, results, net.crs, log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
