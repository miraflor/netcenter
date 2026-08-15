"""Regression tests for the GIS/topology boundary."""

import numpy as np
import pytest
from scipy.sparse.csgraph import connected_components

gpd = pytest.importorskip("geopandas")
shapely_geom = pytest.importorskip("shapely.geometry")
LineString = shapely_geom.LineString

from netcenter.graph import build_network, snap_points  # noqa: E402

UTM51N = "EPSG:32651"


def gdf(geoms, crs=UTM51N):
    return gpd.GeoDataFrame(geometry=list(geoms), crs=crs)


def test_crossings_are_not_connected_by_default():
    roads = gdf(
        [
            LineString([(0, 50), (100, 50)]),
            LineString([(50, 0), (50, 100)]),
        ]
    )
    net = build_network(roads, keep_largest_component=False)
    n_components, _ = connected_components(net.csr, directed=False)
    assert net.n_nodes == 4
    assert net.n_edges == 2
    assert n_components == 2


def test_planar_noding_is_explicit_and_splits_crossings():
    roads = gdf(
        [
            LineString([(0, 50), (100, 50)]),
            LineString([(50, 0), (50, 100)]),
        ]
    )
    with pytest.warns(UserWarning, match="planar noding"):
        net = build_network(roads, node=True)
    assert net.n_nodes == 5
    assert net.n_edges == 4
    assert net.edge_len.sum() == pytest.approx(200.0)


def test_shape_points_do_not_become_nodes():
    curve = LineString([(x, 10 * np.sin(x / 40)) for x in range(0, 400, 5)])
    net = build_network(gdf([curve]))
    assert net.n_edges == 1
    assert net.n_nodes == 2


def test_closed_loop_is_kept_not_deleted():
    ring_pts = [
        (1000 + 300 * np.cos(a), 300 * np.sin(a)) for a in np.linspace(0, 2 * np.pi, 33)
    ]
    expected = 1300.0 + LineString(ring_pts).length
    net = build_network(gdf([LineString([(0, 0), (1300, 0)]), LineString(ring_pts)]))
    assert net.edge_len.sum() == pytest.approx(expected, abs=0.01)
    assert (net.edge_u != net.edge_w).all()


def test_empty_geometry_does_not_shift_endpoints():
    net = build_network(
        gdf(
            [
                LineString([(0, 0), (100, 0)]),
                LineString(),
                LineString([(100, 0), (200, 0)]),
            ]
        )
    )
    assert net.edge_len == pytest.approx([100.0, 100.0])


def test_projected_metre_input_is_not_reprojected():
    net = build_network(gdf([LineString([(0, 0), (1300, 0)])]))
    assert net.crs.to_string() == UTM51N
    assert net.edge_len[0] == pytest.approx(1300.0)


def test_latlon_input_is_reprojected_to_metres():
    net = build_network(gdf([LineString([(121.0, 14.6), (121.01, 14.6)])], crs="EPSG:4326"))
    assert not net.crs.is_geographic
    assert 1000 < net.edge_len[0] < 1200


def test_explicit_geographic_target_crs_is_rejected():
    with pytest.raises(ValueError, match="projected"):
        build_network(gdf([LineString([(0, 0), (100, 0)])]), target_crs="EPSG:4326")


def test_explicit_non_metre_target_crs_is_rejected():
    with pytest.raises(ValueError, match="metres"):
        build_network(gdf([LineString([(0, 0), (100, 0)])]), target_crs="EPSG:2263")


def test_largest_component_is_chosen_by_road_length_not_node_count():
    # Dense but tiny component: 10 one-metre links.
    tiny = [LineString([(i, 0), (i + 1, 0)]) for i in range(10)]
    # Sparse but substantively larger component: one 100-metre road.
    long = [LineString([(1000, 0), (1100, 0)])]
    with pytest.warns(UserWarning, match="most road length"):
        net = build_network(gdf(tiny + long))
    assert net.n_nodes == 2
    assert net.n_edges == 1
    assert net.edge_len.sum() == pytest.approx(100.0)


def test_missing_crs_is_rejected():
    with pytest.raises(ValueError, match="coordinate reference system"):
        build_network(gdf([LineString([(0, 0), (1, 0)])], crs=None))


def test_invalid_snap_is_rejected():
    with pytest.raises(ValueError, match="snap"):
        build_network(gdf([LineString([(0, 0), (1, 0)])]), snap=0)


def test_snap_points_respects_max_dist():
    net = build_network(gdf([LineString([(0, 0), (100, 0)])]))
    idx, dist = snap_points(net, np.array([[0.0, 5.0], [0.0, 5000.0]]), max_dist=50.0)
    assert idx[0] >= 0 and idx[1] == -1
    assert dist[0] == pytest.approx(5.0)


def test_snap_points_rejects_nonfinite_coordinates():
    net = build_network(gdf([LineString([(0, 0), (100, 0)])]))
    with pytest.raises(ValueError, match="NaN or infinity"):
        snap_points(net, np.array([[np.nan, 5.0]]))


def test_end_to_end_solve_finds_known_midpoint():
    from netcenter.solve import solve

    net = build_network(gdf([LineString([(0, 0), (3, 0)]), LineString([(3, 0), (10, 0)])]))
    ends, _ = snap_points(net, np.array([[0.0, 0.0], [10.0, 0.0]]))
    res = solve(net, demand_nodes=ends)
    assert res["vertex_center"].objective == pytest.approx(7.0)
    assert res["absolute_center"].objective == pytest.approx(5.0)
    assert res["absolute_center"].xy[0] == pytest.approx(5.0)


def test_shared_vertex_split_recovers_osm_style_junctions():
    """OSM stores a through-road as one way whose interior vertices ARE the
    junctions, with side streets ending on them. Endpoint-only matching misses
    most of those connections in this synthetic network."""
    lines, n = [], 8
    for i in range(n):
        lines.append(LineString([(x * 100, i * 100) for x in range(n)]))
    for j in range(1, n - 1):
        for i in range(n - 1):
            lines.append(LineString([(j * 100, i * 100), (j * 100, (i + 1) * 100)]))
    g = gdf(lines)
    total = sum(line.length for line in lines)

    kept = build_network(g).edge_len.sum()
    assert kept == pytest.approx(total)

    with pytest.warns(UserWarning, match="disconnected"):
        shattered = build_network(g, split_shared_vertices=False).edge_len.sum()
    assert shattered < 0.2 * total


def test_shared_vertex_split_does_not_weld_a_bridge():
    """An overpass shares no vertex with the road beneath it, so it must stay
    separate. This is what planar noding gets wrong."""
    from scipy.sparse.csgraph import connected_components

    bridge = gdf(
        [
            LineString([(0, 50), (100, 50)]),
            LineString([(50, 0), (50, 100)]),
        ]
    )
    net = build_network(bridge, keep_largest_component=False)
    assert connected_components(net.csr, directed=False)[0] == 2

    with pytest.warns(UserWarning, match="planar noding"):
        planar = build_network(bridge, node=True, keep_largest_component=False)
    assert connected_components(planar.csr, directed=False)[0] == 1


def test_shared_vertex_split_leaves_untouched_lines_alone():
    net = build_network(gdf([LineString([(0, 0), (50, 0), (100, 0)])]))
    assert net.n_edges == 1
    assert net.edge_len[0] == pytest.approx(100.0)


def test_shared_interior_vertices_are_split_on_both_lines():
    """Two continuing ways can share a real OSM junction at interior vertices."""
    roads = gdf(
        [
            LineString([(0, 50), (50, 50), (100, 50)]),
            LineString([(50, 0), (50, 50), (50, 100)]),
        ]
    )
    net = build_network(roads, keep_largest_component=False)
    n_components, _ = connected_components(net.csr, directed=False)
    assert n_components == 1
    assert net.n_nodes == 5
    assert net.n_edges == 4
    assert net.edge_len.sum() == pytest.approx(200.0)


def test_filtered_sliver_does_not_leave_an_orphan_node():
    roads = gdf(
        [
            LineString([(0, 0), (100, 0)]),
            LineString([(1000, 0), (1000.00000001, 0)]),
        ]
    )
    net = build_network(roads, min_length=1e-3, keep_largest_component=False)
    assert net.n_nodes == 2
    assert net.n_edges == 1
    assert net.csr.shape == (2, 2)


def test_absurdly_small_snap_is_rejected_before_integer_overflow():
    with pytest.raises(ValueError, match="snap is too small"):
        build_network(gdf([LineString([(500000, 0), (500100, 0)])]), snap=1e-20)
