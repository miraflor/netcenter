# Changelog

All notable changes to this repository are recorded here.

## 0.1.0 — 2026-09-24

First public release of `netcenter`.

### Facility-location solvers

- Weighted network 1-median for vertex demand.
- Unweighted vertex 1-center.
- Exact unweighted absolute 1-center, allowing the optimum to lie inside a network edge.
- Exact coalescing of repeated snapped demand nodes before shortest-path calculation.
- Safe lower-bound pruning before the continuous absolute-center edge sweep.
- Continuous edgewise breakpoint/envelope solution based on the classical network-location formulation of Hakimi (1964) and Kariv & Hakimi (1979).

### Network construction and topology

- GeoPandas/Shapely linework to sparse undirected network conversion.
- Shared-source-vertex noding as the default topology rule.
- Recovery of both T-junctions and shared interior-interior junctions already encoded in the source linework.
- Optional full planar noding with `--node-crossings`.
- Strict endpoint-only mode with `--no-shared-vertex-noding` for networks already segmented at all true junctions.
- Geometric crossings are not assumed to be connected by default, reducing false junctions at bridges, flyovers, tunnels, and other grade-separated crossings.
- Parallel edges between the same node pair are reduced using the minimum direct edge length for shortest-path calculations.
- Closed and nearly closed rings are preserved rather than discarded as self-loops.
- Sliver and self-loop filtering is followed by node-ID compaction so discarded edges do not leave orphan graph nodes.
- Optional connected-component filtering based on total retained road length.

### Demand handling

- Point and polygon demand layers.
- Representative points are used for non-point demand geometries.
- Demand points are transformed into the working network CRS before snapping.
- Nearest-network-node snapping.
- Optional maximum snapping distance with `--max-snap`.
- Exact aggregation of repeated observations snapped to the same network node.
- Optional non-negative demand weights for the 1-median objective.
- Duplicate demand locations are removed from the minimax center calculation because multiplicity does not affect a maximum.

### Coordinate reference systems and units

- All reported network distances, snapping distances, edge lengths, and edge offsets are expressed in metres.
- Geographic input is reprojected to an inferred local UTM CRS unless an explicit working CRS is supplied.
- Projected input already using metres is preserved.
- Projected input using other units is reprojected.
- Explicit target CRSs are validated to ensure they are projected and metre-based.

### Shortest-path computation

- Sparse shortest-path distances computed using SciPy's Dijkstra implementation.
- Demand-by-network shortest-path matrix reused across the median and center solvers.
- Bounded Dijkstra source blocks to control temporary memory use.
- Optional float32 storage to approximately halve distance-matrix memory requirements.
- Parallel shortest-path blocks are streamed directly into a preallocated result matrix instead of being accumulated and concatenated afterward.
- Worker counts and memory-control parameters are validated explicitly.

### Absolute-center computation

- Exact continuous search along candidate network edges rather than restricting the minimax solution to existing vertices.
- For edge `(u, w)` of length `L`, distance from an interior point at offset `t` to demand `v` is evaluated as

  ```text
  min(d(u, v) + t, d(w, v) + L - t)
  ```

- Demand-specific route-switch breakpoints are sorted along each candidate edge.
- Prefix and suffix running maxima reduce each breakpoint interval to the minimum of a two-line upper envelope.
- Segment-level lower bounds prune edges that cannot improve the incumbent vertex-center solution.
- Scale-aware numerical margins prevent float32 rounding from incorrectly eliminating borderline candidate edges.
- Continuous results can be returned as an edge index plus an offset in metres and are interpolated back to map coordinates.

### Performance and memory

- Conservative default of one worker.
- Optional process-based and thread-based parallel backends.
- Separate backend controls for shortest-path computation and absolute-center edge sweeps.
- Bounded-memory segment-demand sweep blocks through `--max-cells`.
- Bounded temporary Dijkstra blocks through `--max-temp-mb`.
- Distance-matrix memory estimation.
- Duplicate-demand coalescing reduces redundant shortest-path computation.

### Command-line interface

- Road-network input from standard GIS line formats supported by GeoPandas.
- Optional road layer selection.
- Optional demand layer and demand layer-name selection.
- Optional weight field.
- Optional output to a spatial point file.
- Optional suppression of the continuous absolute-center calculation with `--skip-absolute`.
- Explicit topology controls for shared-vertex, endpoint-only, and planar interpretations.
- Quiet mode for reduced console output.
- Package version reporting through `--version`.

### Python API

- `build_network()` for network construction.
- `snap_points()` for assigning external coordinates to network nodes.
- `solve()` for computing the weighted 1-median, vertex 1-center, and absolute 1-center.
- Package version read from installed package metadata.

### Validation and testing

- Hand-computable weighted median cases.
- Hand-computable vertex-center cases.
- Edge-interior absolute-center cases.
- Independent brute-force edge sampling used to validate the continuous center sweep on random networks.
- Serial versus parallel invariance tests.
- Block-size invariance tests.
- Float32 behavior tests.
- T-junction topology regression tests.
- Shared interior-interior junction regression tests.
- Bridge versus true-junction regression tests.
- Parallel-edge regression tests.
- Closed-ring and nearly closed-ring regression tests.
- CRS and unit-validation tests.
- Empty-geometry and sliver handling tests.
- Disconnected-network tests.
- Worker-count, memory-control, and package-metadata tests.

### Documentation

- `README.md` with installation, usage, mathematical interpretation, and scope.
- `docs/ALGORITHMS.md` with algorithmic derivations and literature attribution.
- `docs/ARCHITECTURE.md` describing the package structure.
- `docs/NETWORK_ASSUMPTIONS.md` documenting topology and modeling assumptions.
- `docs/VALIDATION.md` documenting the testing and validation strategy.
- `docs/REVIEW_NOTES.md` recording pre-release review and implementation hardening.
- `docs/TECHNICAL_NOTE.tex` and rendered PDF containing the detailed mathematical and implementation discussion.

### Algorithmic lineage

The facility-location mathematics implemented by `netcenter` is classical.

Principal references include:

- Hakimi, S. L. (1964). “Optimum Locations of Switching Centers and the Absolute Centers and Medians of a Graph.” *Operations Research*, 12(3), 450–459. https://doi.org/10.1287/opre.12.3.450
- Hakimi, S. L. (1965). “Optimum Distribution of Switching Centers in a Communication Network and Some Related Graph Theoretic Problems.” *Operations Research*, 13(3), 462–475. https://doi.org/10.1287/opre.13.3.462
- Kariv, O., & Hakimi, S. L. (1979). “An Algorithmic Approach to Network Location Problems. I: The p-Centers.” *SIAM Journal on Applied Mathematics*, 37(3), 513–538. https://doi.org/10.1137/0137040
- Kariv, O., & Hakimi, S. L. (1979). “An Algorithmic Approach to Network Location Problems. II: The p-Medians.” *SIAM Journal on Applied Mathematics*, 37(3), 539–560. https://doi.org/10.1137/0137041
- Handler, G. Y., & Mirchandani, P. B. (1979). *Location on Networks: Theory and Algorithms*. MIT Press.
- Daskin, M. S. (2013). *Network and Discrete Location: Models, Algorithms, and Applications* (2nd ed.). Wiley.

The package does not claim a new network-location algorithm. Its additional work is primarily in GIS topology construction, numerical safeguards, vectorization, memory management, pruning, parallel execution, and validation.