# netcenter

**Exact 1-centre and weighted 1-median locations on an undirected road network.**

`netcenter` answers three related facility-location questions:

| Result | Objective | Candidate locations |
|---|---|---|
| weighted 1-median | minimise total weighted road distance | network nodes |
| vertex centre | minimise worst road distance | network nodes |
| absolute 1-centre | minimise worst road distance | **any point on the network** |

The last case is the distinctive one. The optimum can lie in the interior of a
road segment, so checking junctions alone can return the wrong minimax location.

## Why v0.3.0

This release is primarily a **topology and memory-safety release**.

The important change is a three-level distinction between ways of interpreting
road linework:

1. **Default — shared-vertex noding.** If two input LineStrings already contain
   the same source vertex, that vertex is treated as a real junction and both
   lines are split there. This recovers common OSM-style T- and X-junctions.
2. **Strict endpoint-only mode.** Use `split_shared_vertices=False` (or
   `--no-shared-vertex-noding`) only when the source is already segmented at
   every real junction.
3. **Explicit planar noding.** Use `node=True` / `--node-crossings` only when
   every geometric crossing is genuinely connected. This can incorrectly weld
   a bridge to the road beneath it.

The default therefore preserves topology already encoded by the data without
inventing turns at arbitrary map crossings.

Other v0.3.0 changes:

- genuine shared **interior-interior** junctions are now recovered, not only
  side-street endpoints touching a through-road;
- nodes left orphaned after sliver/self-loop filtering are removed before the
  sparse graph is built;
- parallel shortest-path blocks stream directly into the final matrix instead
  of being accumulated and copied with `numpy.vstack`;
- worker counts and memory controls are validated instead of silently coerced;
- package version metadata, CLI defaults, documentation, and tests are aligned;
- one worker remains the conservative default; parallel backends are opt-in.

See [`CHANGELOG.md`](CHANGELOG.md) and
[`docs/REVIEW_NOTES.md`](docs/REVIEW_NOTES.md) for the detailed review.

---

## Repository map

```text
netcenter/
├─ .github/
│  └─ workflows/
│     └─ ci.yml
├─ docs/
│  ├─ ALGORITHMS.md
│  ├─ ARCHITECTURE.md
│  ├─ NETWORK_ASSUMPTIONS.md
│  ├─ REVIEW_NOTES.md
│  ├─ TECHNICAL_NOTE.tex
│  ├─ TECHNICAL_NOTE.pdf
│  └─ VALIDATION.md
├─ examples/
│  └─ quickstart.py
├─ netcenter/
│  ├─ __init__.py
│  ├─ center.py
│  ├─ cli.py
│  ├─ distances.py
│  ├─ graph.py
│  ├─ solve.py
│  └─ topology.py
├─ tests/
│  ├─ test_center.py
│  ├─ test_distances.py
│  ├─ test_graph.py
│  ├─ test_package.py
│  └─ test_solve.py
├─ .gitignore
├─ CHANGELOG.md
├─ GITHUB_SETUP.md
├─ LICENSE
├─ README.md
└─ pyproject.toml
```

If you want to understand the code rather than merely run it, read in this
order:

```text
solve.py -> center.py -> distances.py -> graph.py -> topology.py
```

---

## Installation

### Full GIS installation

```bash
python -m pip install -e ".[gis]"
```

### Development installation

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
```

### Numerical core only

If you already have a sparse graph and do not need GeoPandas/Shapely file I/O:

```bash
python -m pip install -e .
```

---

## Command-line use

Minimal:

```bash
netcenter roads.gpkg
```

With demand points and weights:

```bash
netcenter roads.gpkg \
  --layer roads \
  --demand barangays.gpkg \
  --weight-field population \
  --max-snap 1000 \
  --float32 \
  --out centers.gpkg
```

### Topology options

Default, recommended for OSM-shaped/routing-quality linework:

```bash
netcenter roads.gpkg
```

Disable shared-vertex splitting only if the source is already segmented at all
true junctions:

```bash
netcenter roads.gpkg --no-shared-vertex-noding
```

Treat every drawn crossing as connected only for genuinely planar data:

```bash
netcenter roads.gpkg --node-crossings
```

### Parallel execution

The safe default is one worker:

```bash
netcenter roads.gpkg
```

For a large problem, try more workers deliberately:

```bash
netcenter roads.gpkg --jobs 4 --backend loky
```

Or, on a memory-constrained system, compare the shared-memory backend:

```bash
netcenter roads.gpkg --jobs 4 --backend threading
```

Parallel speedups vary by graph size, number of demand nodes, SciPy build,
operating system, and available RAM. Benchmark the actual workload rather than
assuming a backend is universally faster.

---

## Python use

```python
from netcenter import build_network, snap_points, solve

net = build_network("roads.gpkg")

# Coordinates must be in net.crs.
demand_xy = [
    [305000.0, 1615000.0],
    [306250.0, 1616200.0],
]

nodes, snap_distance = snap_points(net, demand_xy, max_dist=1000)

results = solve(
    net,
    demand_nodes=nodes,
    weights=[1200, 800],
    dtype="float32",
)

print(results["median"])
print(results["vertex_center"])
print(results["absolute_center"])
```

The absolute-centre result is either a node or an edge plus an offset `t` in
metres. `solve()` attaches the corresponding map coordinate in `.xy`.

---

## Demand semantics

Demand is represented at **network nodes**.

`snap_points()` finds the nearest node, not the nearest arbitrary point along an
edge. This makes the mathematical problem explicit and allows repeated demand
locations to be consolidated exactly, but it also means a simplified road
network can move a demand point farther than expected.

Inspect the reported snap distances and use `--max-snap` in production work.

Repeated demand observations are handled exactly:

- multiplicity does not affect a maximum, so the centre uses each unique demand
  node once;
- multiplicity does affect the median, so counts or supplied weights are summed
  at each unique snapped node.

`weights` affect the **median only**. The centre objective is currently
unweighted over the chosen demand locations.

---

## CRS and distance units

All reported network distances, edge lengths, snap grid spacing, and offsets are
in **metres**.

- projected metre-based input is preserved;
- geographic input is moved to an inferred local UTM CRS unless you supply a
  target CRS;
- projected input using feet or another unit is also reprojected;
- an explicit `target_crs` must be projected and metre-based.

For a large multi-zone study area, supply a projection appropriate to the full
extent instead of relying on one inferred UTM zone.

---

## Memory model

For `k` unique demand nodes and `n` network nodes, the stored shortest-path
matrix is approximately:

```text
float64: 8 x k x n bytes
float32: 4 x k x n bytes
```

Example: 2,000 demand nodes by 12,000 network nodes is about 183 MiB in float64
or 92 MiB in float32.

`distance_matrix()` computes SciPy Dijkstra output in bounded row chunks. The
`max_temp_mb` / `--max-temp-mb` setting controls the target float64 result size
of **one chunk per worker**. With several workers, several chunks can exist at
once, so reduce the budget when RAM is the bottleneck.

Parallel results are streamed into a preallocated final matrix. The previous
implementation retained all returned blocks and then used `numpy.vstack`, which
could transiently duplicate most of the distance matrix.

`float32` is a memory/precision trade-off. The centre solver widens its pruning
margin according to the stored dtype so rounding cannot incorrectly eliminate a
borderline candidate edge, but the final objective is still limited by the
precision of the stored matrix.

---

## Absolute 1-centre in one paragraph

For a point `t` metres along edge `(u, w)` of length `L`, distance to demand
node `v` is

```text
d(t, v) = min(d(u, v) + t, d(w, v) + L - t).
```

Each demand therefore contributes a tent-shaped function along the edge. The
worst-demand distance is the upper envelope of those functions. Once the switch
points are sorted, every interval reduces to the maximum of one rising and one
falling line, whose minimum is analytic. A lower bound first eliminates edges
that cannot improve the best vertex centre, so only surviving edges need the
exact sweep.

See [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) for the compact derivation and references, or [`docs/TECHNICAL_NOTE.tex`](docs/TECHNICAL_NOTE.tex) for the full LaTeX discussion and explicit audit of how v0.3.0 differs from the initial prototype.

---

## Scope and non-goals

Implemented:

- undirected non-negative road costs;
- vertex demand;
- weighted 1-median;
- unweighted vertex centre;
- exact unweighted absolute 1-centre;
- GIS linework -> sparse network conversion;
- shared-source-vertex topology recovery;
- optional explicit planar noding.

Not implemented:

- one-way streets or directed graphs;
- turn restrictions;
- time-dependent or asymmetric travel costs;
- weighted minimax centre;
- demand placed continuously along edges;
- automatic inference of bridge/tunnel semantics from attributes;
- out-of-core storage for a distance matrix too large for RAM.

---

## Validation

The development suite contains **69 tests** in this release. It includes:

- hand-computable centre/median cases;
- random-network comparison with independent brute-force edge sampling;
- serial/parallel and block-size invariance checks;
- bridge versus true-junction topology regressions;
- T-junction and interior-interior shared-vertex regressions;
- CRS/unit, closed-ring, empty-geometry, sliver, and component regressions;
- distance chunking, float32, worker validation, and package-metadata checks.

Run:

```bash
python -m pytest -q
```

For the validation philosophy, see [`docs/VALIDATION.md`](docs/VALIDATION.md).

---

## License

MIT. See [`LICENSE`](LICENSE).
