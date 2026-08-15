import pytest
from shapely.geometry import LineString

pytest.importorskip("geopandas")
import geopandas as gpd

from netcenter.graph import build_network
from netcenter.solve import _coalesce_demands, solve


def make_path():
    gdf = gpd.GeoDataFrame(
        geometry=[LineString([(0, 0), (1, 0)]), LineString([(1, 0), (2, 0)])],
        crs="EPSG:32651",
    )
    return build_network(gdf)


def test_duplicate_demands_are_coalesced_with_counts():
    nodes, weights = _coalesce_demands([0, 0, 2], None, 3)
    assert nodes.tolist() == [0, 2]
    assert weights.tolist() == [2.0, 1.0]


def test_duplicate_weighted_demands_are_summed():
    nodes, weights = _coalesce_demands([0, 0, 2], [2.0, 3.0, 4.0], 3)
    assert nodes.tolist() == [0, 2]
    assert weights.tolist() == [5.0, 4.0]


def test_coalescing_preserves_median_multiplicity():
    net = make_path()
    res = solve(net, demand_nodes=[0, 0, 2])
    assert res["median"].node == 0
    assert res["median"].objective == pytest.approx(2.0)
    assert res["absolute_center"].objective == pytest.approx(1.0)


def test_weight_length_is_checked_before_shortest_paths():
    net = make_path()
    with pytest.raises(ValueError, match="weights for"):
        solve(net, demand_nodes=[0, 2], weights=[1.0])


def test_public_solve_symbol_is_callable():
    from netcenter import solve as public_solve

    assert callable(public_solve)


def test_fractional_demand_node_indices_are_rejected_instead_of_truncated():
    net = make_path()
    with pytest.raises(ValueError, match="integer node indices"):
        solve(net, demand_nodes=[0.5, 2.0])
