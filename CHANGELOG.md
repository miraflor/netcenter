# Changelog

All notable changes to this repository are recorded here.

## 0.3.0 — 2026-08-15

### Fixed

- Junction recovery now splits at any snapped coordinate shared by two or more
  distinct input LineStrings, not only where one line's endpoint meets another
  line's interior vertex. The narrower endpoint-only rule missed X-junctions
  where both ways continue through the shared node; on a grid of overshooting
  ways it retained 6.2% of the network versus 100% now. Bridge safety is
  unchanged: lines that merely cross in 2-D share no source vertex and stay
  disconnected.
- Node IDs are compacted immediately after edge filtering, so a node left
  behind by a discarded sliver no longer survives as an isolated matrix row
  producing spurious infinite distances when `keep_largest_component=False`.
- The parallel shortest-path path preallocates its output and writes each block
  into place, instead of holding every block plus a second stacked copy.
- `n_jobs` values of zero or below are rejected rather than silently coerced.
- Snap-grid integer conversion checks for overflow before casting.

### Changed

- Shortest paths and the centre sweep take separate backend settings
  (`backend`, `sweep_backend`), since the two stages differ in memory profile
  and GIL behaviour. Both `loky` and `threading` are documented as options to
  benchmark rather than one being assumed faster.
- `n_jobs` defaults to 1. Raise it explicitly after benchmarking on the target
  machine; the work-size thresholds keep small jobs serial regardless.
- `__version__` is read from installed package metadata.

### Added

- `docs/TECHNICAL_NOTE.tex` / `.pdf`: 14-page derivation, complexity analysis,
  and an explicit audit of where v0.3.0 changes graph construction rather than
  the facility-location mathematics.
- `docs/REVIEW_NOTES.md`, `docs/NETWORK_ASSUMPTIONS.md`, `docs/VALIDATION.md`.

## 0.2.0 — 2026-08-15

### Correctness and modelling safety

- Changed road-building default from planar noding to topology-preserving input:
  geometric crossings are no longer assumed to be valid turns.
- Added explicit `--node-crossings` for genuinely planar linework.
- Added projected/metre CRS validation.
- Geographic and non-metre projected inputs are converted to an inferred local
  UTM CRS when no explicit working CRS is supplied.
- Added warning for geographic study areas spanning more than one typical UTM
  zone width.
- Changed largest-component selection from node count to total road length.
- Fixed nearly closed rings whose first and last coordinates differ only by
  floating-point noise.
- Added validation for invalid edge lengths, distance matrices, weights, node
  indices, snap tolerances, and non-finite demand coordinates.

### Performance and memory

- Repeated snapped demand nodes are coalesced before shortest-path calculation.
- Median multiplicity/weights are preserved exactly during coalescing.
- Shortest-path calculation now uses bounded Dijkstra row chunks.
- Added `max_temp_mb` / `--max-temp-mb`.
- Added distance-matrix memory estimation.
- Changed default worker count to 1 for predictable laptop-safe behaviour.
- Changed default parallel backend to `threading` so workers can share the graph
  and large arrays rather than duplicating them across processes.

### Repository structure

- Added `CHANGELOG.md`.
- Added `GITHUB_SETUP.md` with replacement-repo instructions.
- Moved detailed algorithm notes under `docs/`.
- Added architecture, network-assumption, and validation documentation.
- Added a minimal Python example.
- Expanded CI and tests.

### Tests

- Expanded from 33 to 52 passing tests.
- Added regression coverage for false planar intersections, CRS misuse,
  component selection, nearly closed rings, demand coalescing, Dijkstra
  chunking, and additional invalid-input cases.

## 0.1.0

- Initial packaged implementation of weighted 1-median, vertex centre, and
  exact absolute 1-centre.
- Initial GeoPandas/Shapely road-network builder and command-line interface.
