"""Small packaging-level regressions that catch release metadata drift."""

from importlib.metadata import PackageNotFoundError, version

import pytest

import netcenter


def test_exported_version_matches_installed_distribution():
    # A fresh clone has no installed distribution, so this would otherwise fail
    # for a contributor who simply runs pytest. CI installs first, where the
    # check is meaningful and does run.
    try:
        installed = version("netcenter")
    except PackageNotFoundError:
        pytest.skip("netcenter is not installed; run `pip install -e .` first")
    assert netcenter.__version__ == installed


def test_public_solve_is_callable():
    assert callable(netcenter.solve)
