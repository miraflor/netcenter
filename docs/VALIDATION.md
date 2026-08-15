# Validation strategy

The dangerous failures in a network-location package are usually plausible
wrong answers rather than crashes. The test suite therefore mixes mathematical
checks, GIS/topology regressions, memory/execution checks, and packaging checks.

## Mathematical checks

- hand-computable path and triangle examples;
- weighted-median movement under changed weights;
- twelve random networks compared against independent brute-force sampling at
  4,001 locations on every segment;
- the absolute centre may never be worse than the vertex centre;
- serial/parallel and large/small sweep blocks must agree;
- float32 storage is checked on a known case;
- invalid matrices, edges, weights, and worker settings fail loudly.

## GIS/topology regressions

- a plain geometric crossing is disconnected by default;
- explicit planar noding connects that crossing;
- a shared endpoint-to-interior T-junction is recovered;
- a shared interior-to-interior X-junction is recovered;
- a bridge-like crossing with no shared source vertex remains disconnected;
- ordinary shape vertices on a single road do not become graph nodes;
- closed and nearly closed rings survive preprocessing;
- empty geometries do not shift endpoint indexing;
- filtered slivers do not leave orphaned graph nodes;
- largest-component selection follows total road length;
- metric CRS behavior and invalid target CRSs are tested;
- non-finite demand coordinates and invalid snap grid settings are rejected;
- extreme snap scaling is rejected before integer overflow.

## Distance/memory checks

- bounded Dijkstra chunking matches effectively unbounded calculation;
- float32 storage is supported;
- matrix-size estimates are tested;
- threaded parallel chunks match serial output;
- invalid temporary-memory and worker settings are rejected;
- duplicate demand nodes are consolidated while preserving median counts and
  supplied weights.

## Packaging checks

- the public `solve` symbol remains callable;
- exported `netcenter.__version__` matches installed package metadata;
- CI lints, compiles, tests, exercises the CLI, and builds the package.

## Run locally

```bash
python -m pytest -q
```

A clean v0.3.0 tree reports **69 passed** with the development dependencies
installed.
