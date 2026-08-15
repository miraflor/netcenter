"""Small, self-contained example that does not require a GIS file on disk."""

import geopandas as gpd
from shapely.geometry import LineString

from netcenter import build_network, solve

# A simple 10 m path with one junction at x=3.
roads = gpd.GeoDataFrame(
    geometry=[
        LineString([(0, 0), (3, 0)]),
        LineString([(3, 0), (10, 0)]),
    ],
    crs="EPSG:32651",  # projected CRS with metre units
)

net = build_network(roads)

# Demand at the two end nodes. The best vertex centre is x=3 with radius 7 m,
# while the absolute centre is x=5 with radius 5 m.
results = solve(net, demand_nodes=[0, 2])

for name, result in results.items():
    print(name, result)
