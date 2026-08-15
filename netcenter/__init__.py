"""Exact 1-centre and weighted 1-median tools for undirected road networks.

The numerical core depends only on NumPy/SciPy/joblib. GIS functionality is
loaded lazily so users with their own graph and distance matrix do not need the
GeoPandas/Shapely stack.
"""

from importlib.metadata import PackageNotFoundError, version

from netcenter.center import (
    CenterResult,
    absolute_center,
    vertex_center,
    vertex_eccentricity,
    weighted_median,
)
from netcenter.distances import distance_matrix, estimate_distance_matrix_mb
from netcenter.solve import solve

try:
    __version__ = version("netcenter")
except PackageNotFoundError:  # Running directly from an uninstalled source tree.
    __version__ = "0+unknown"

_LAZY = {
    "Network": "netcenter.graph",
    "build_network": "netcenter.graph",
    "snap_points": "netcenter.graph",
}

__all__ = [
    "CenterResult",
    "Network",
    "absolute_center",
    "build_network",
    "distance_matrix",
    "estimate_distance_matrix_mb",
    "snap_points",
    "solve",
    "vertex_center",
    "vertex_eccentricity",
    "weighted_median",
]


def __getattr__(name):
    if name in _LAZY:
        import importlib

        return getattr(importlib.import_module(_LAZY[name]), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
