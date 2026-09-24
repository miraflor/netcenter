# Review notes for v0.1.0

This file records the pre-release review and hardening work that led to the
first public release of `netcenter`, v0.1.0.

The codebase went through several internal development states before public
release. Those states are treated here as **pre-release development history**,
not as a public semantic-version sequence.

The review focused on two different questions:

1. **Are the facility-location algorithms mathematically correct?**
2. **Does the GIS/network construction produce the graph those algorithms are
   supposed to solve?**

Those are separate failure modes. A correct network-center algorithm applied to
the wrong graph can still return a precise but substantively wrong answer.

---

## 1. Mathematical scope remained stable

The mathematical core of the package remained the same throughout pre-release
development:

- weighted network 1-median;
- vertex-restricted 1-center;
- exact absolute 1-center with edge-interior facility location;
- Dijkstra shortest-path distance on an undirected non-negative graph;
- breakpoint-based continuous edge search for the absolute center.

The relevant classical lineage is documented in `docs/ALGORITHMS.md` and
`docs/TECHNICAL_NOTE.tex`, principally:

- Hakimi (1964);
- Hakimi (1965);
- Kariv & Hakimi (1979);
- Handler & Mirchandani (1979);
- Daskin (2013).

The release does **not** claim a new facility-location algorithm.

Most review-driven changes affected graph construction, numerical safety,
memory behavior, execution policy, validation, and documentation.

---

## 2. Topology review

### 2.1 Why topology was the highest-risk area

GIS linework is not automatically a routing graph.

Two roads can cross geometrically without being connected. Conversely, a true
junction can be encoded as a coordinate in the interior of one or more source
LineStrings.

The review therefore treated graph construction as part of the model rather
than as a neutral preprocessing step.

### 2.2 Default topology rule

The final v0.1.0 default is **shared-source-vertex noding**.

A snapped coordinate is treated as a junction when it occurs as a source vertex
in at least two distinct input LineStrings.

This recovers:

- ordinary T-junctions;
- shared interior-interior junctions;
- common OSM-style junction encoding;

without automatically inventing connectivity at every two-dimensional crossing.

### 2.3 Why full planar noding is not the default

Treating every geometric crossing as connected can incorrectly join:

- bridges to roads below them;
- flyovers to surface streets;
- tunnel crossings;
- other grade-separated structures.

Full planar noding remains available explicitly through `--node-crossings`.

That mode is appropriate only when every planar crossing genuinely represents a
routable connection.

### 2.4 Strict endpoint-only mode

`--no-shared-vertex-noding` is retained for already segmented networks where
every true junction is encoded as a line endpoint.

It should not be treated as the general default.

---

## 3. Important correctness bugs found and fixed

Several pre-release failures could produce plausible wrong answers rather than
obvious crashes. These are the most important findings from review.

### 3.1 Parallel sparse edges could be summed

A direct sparse COO-to-CSR conversion can add duplicate entries.

If two physical segments join the same pair of nodes, summing their lengths
creates an artificial direct edge with the combined cost.

That is wrong for shortest-path use and can violate assumptions used by the
absolute-center sweep.

**Fix:** parallel node pairs are collapsed using the **minimum direct edge
length**.

---

### 3.2 Closed loops could disappear

A ring road can merge into a closed LineString whose start and end coordinates
coincide.

If this becomes a self-loop and is filtered, the entire physical ring can
disappear.

**Fix:** closed and nearly closed rings are detected and split before graph
construction.

---

### 3.3 Nearly closed rings could evade exact closure tests

Some GIS rings are visually closed but have start/end coordinates differing by
floating-point noise.

Exact equality is not reliable.

**Fix:** ring closure uses the same tolerance logic as endpoint/node identity.

---

### 3.4 Empty geometries could corrupt endpoint indexing

Unusable geometries can create empty coordinate runs that make downstream
endpoint lookups refer to the wrong data.

**Fix:** empty and unusable geometries are filtered before endpoint extraction.

---

### 3.5 Shared-vertex recovery was initially too narrow

An intermediate rule handled T-junctions where one line ended on another line's
interior vertex.

It still missed cases where two source LineStrings both continued through the
same shared interior coordinate.

**Fix:** a junction is defined by participation of at least two distinct source
LineStrings at the snapped coordinate.

This covers both T and X-style shared-source junctions without planarizing
arbitrary crossings.

---

### 3.6 Filtered slivers could leave orphan graph nodes

Node IDs were originally assigned before all edge filtering completed.

A discarded sliver or snapped self-loop could therefore leave an unused graph
row/column.

That could later create infinite distances in an otherwise connected surviving
network.

**Fix:** node IDs are compacted immediately after edge filtering.

---

### 3.7 CRS handling was too aggressive

Automatically moving all data into an inferred local UTM CRS is appropriate for
geographic coordinates, but unnecessary and potentially undesirable for input
already projected in valid metre units.

**Fix:**

- preserve projected metre-based input;
- reproject geographic input;
- reproject projected input whose units are not metres;
- require explicit target CRSs to be projected and metre-based.

---

## 4. Demand handling review

### 4.1 Demand is vertex-based

External demand observations are snapped to the nearest network node.

This is intentional and documented.

The package does **not** currently solve the continuous-demand-on-edge problem.

### 4.2 Duplicate snapped demand is coalesced exactly

If multiple observations snap to the same node:

- duplicate copies do not matter for the minimax center;
- counts or user-supplied weights do matter for the median.

The release therefore solves shortest paths once per unique snapped node and
aggregates median weights exactly.

This is both faster and mathematically equivalent to the unreduced formulation.

### 4.3 Snap distance should be inspected

A simplified network can move an observation farther than expected when demand
is snapped to nodes rather than arbitrary edge points.

For that reason:

- snap distances are reported;
- `--max-snap` can reject distant observations;
- large snaps should be treated as a modeling warning, not merely a technical
  detail.

---

## 5. Shortest-path review

Shortest paths are computed through SciPy's sparse Dijkstra implementation.

The review added or strengthened:

- square/symmetric graph validation;
- non-negative edge-cost assumptions;
- bounded source blocks;
- explicit temporary-memory control;
- streaming assembly into a preallocated output matrix;
- optional float32 storage.

The dominant practical memory cost is often the `k × n` shortest-path matrix,
not the final center sweep.

---

## 6. Absolute-center review

### 6.1 Edgewise distance formula

For a point `t` metres along edge `(u, w)` of length `L`, and demand `v`:

```text
d(t, v) = min(d(u, v) + t, d(w, v) + L - t)
```

Each demand contributes a tent-shaped piecewise-linear function.

The edge eccentricity is the upper envelope of those tents.

### 6.2 Breakpoint structure

Each demand switches preferred endpoint at

```text
t* = (d(w, v) - d(u, v) + L) / 2
```

Sorting those breakpoints fixes the active rising/falling sets on each interval.

The interval objective reduces to:

```text
max(A + t, B + L - t)
```

which has an analytic minimum.

### 6.3 Prefix/suffix implementation

The release computes:

- suffix maxima for the `u`-side terms;
- prefix maxima for the `w`-side terms.

This is a vectorized implementation choice built around the classical
Kariv–Hakimi edge structure.

### 6.4 Safe lower-bound pruning

For edge `(u, w)`:

```text
LB = max_i min(D[i,u], D[i,w])
```

is a valid lower bound on eccentricity anywhere on the edge.

If that lower bound cannot improve the current best vertex-center radius, the
edge is skipped.

This is an exact pruning argument in real arithmetic.

### 6.5 Numerical tolerance

With finite-precision distance storage, especially float32, pruning must not
discard a candidate merely because of rounding.

The final release therefore uses a dtype- and scale-aware margin.

---

## 7. Memory and parallelism review

### 7.1 Conservative default

`n_jobs=1` is the default.

The package does not assume that more workers are always faster.

### 7.2 Separate backends

Shortest-path computation and the absolute-center sweep expose separate backend
choices because their memory and execution characteristics differ.

### 7.3 Bounded Dijkstra blocks

SciPy Dijkstra output is generated in bounded source chunks.

This prevents one full float64 temporary matrix from coexisting unnecessarily
with the final output matrix.

### 7.4 Streaming parallel output

Parallel blocks are written directly into their final rows.

The implementation avoids holding all completed blocks and then creating a
second full matrix through `numpy.vstack`.

### 7.5 Bounded center-sweep blocks

The absolute-center sweep also processes bounded segment-demand cell blocks.

`max_cells` changes memory usage, not the mathematical result.

---

## 8. Optimization tried and rejected

A best-first bound-tightening strategy was implemented during pre-release work.

The idea was:

1. rank promising edges;
2. sweep the best candidates first;
3. tighten the incumbent center radius;
4. rerun pruning on the remaining edges.

This is standard branch-and-bound practice.

It was removed after measurement because the first lower-bound pass already
reduced the edge set enough that the second pass did not justify the extra
branching and code complexity.

The result is intentionally documented so future optimization work starts from
measurement rather than repeating the same assumption.

---

## 9. AI-assisted implementation work

Some engineering work was developed in AI-assisted programming sessions.

This includes:

- vectorized prefix/suffix envelope evaluation;
- bounded sweep blocks;
- contiguous sort inputs;
- streaming parallel matrix assembly;
- duplicate-demand reduction;
- pruning and tolerance hardening;
- several review-driven regression fixes.

These are implementation details around published algorithms.

They are not presented as new network-location theory and should not be cited as
such.

---

## 10. Validation review

### 10.1 Hand-computable cases

Tests include small networks with known:

- weighted medians;
- vertex centers;
- edge-interior absolute centers.

### 10.2 Independent brute-force comparison

Random test networks are independently sampled along edges.

The exact sweep must be at least as good as the sampled approximation and agree
within sampling resolution.

This is useful because the brute-force check does not share the breakpoint
implementation.

### 10.3 Important limitation of brute-force validation

Brute-force edge sampling validates the **optimizer conditional on the supplied
graph**.

It does not validate whether the graph represents the intended real transport
network.

For example, a numerical test cannot determine whether a bridge crossing should
be routable from geometry alone.

That remains a data/topology responsibility.

### 10.4 Regression coverage

Regression tests cover:

- false planar intersections;
- T-junction recovery;
- shared interior-interior junction recovery;
- parallel edges;
- rings and nearly closed rings;
- empty geometries;
- filtered slivers;
- disconnected graphs;
- CRS/unit misuse;
- shortest-path chunking;
- float32 behavior;
- serial/parallel invariance;
- block-size invariance;
- worker validation;
- package metadata and CLI version consistency.

---

## 11. What changed from the earliest prototype?

The key distinction is:

### Same optimization problem and mathematical method

The following remained stable:

- weighted 1-median objective;
- vertex minimax center;
- absolute 1-center objective;
- Dijkstra network metric;
- edgewise breakpoint/envelope structure;
- lower-bound pruning logic.

### Different end-to-end graph construction and execution behavior

The following were hardened or changed:

- default crossing semantics;
- shared-source-vertex recovery;
- parallel-edge handling;
- ring preservation;
- CRS handling;
- component selection;
- demand coalescing;
- shortest-path memory policy;
- numerical tolerance;
- parallel scheduling.

These changes can move the reported center because they change the graph or
stored metric being solved, not because they redefine the facility-location
objective.

---

## 12. Release audit table

| Aspect | Early prototype | v0.1.0 | Nature of change |
|---|---|---|---|
| Weighted 1-median | Weighted total shortest-path distance over vertices | Same | No mathematical change |
| Vertex 1-center | Minimax over vertices | Same | No mathematical change |
| Absolute 1-center | Continuous edge-interior minimax | Same | No mathematical change |
| Shortest paths | Dijkstra | Same, with validation/chunking | Engineering only |
| Lower-bound pruning | Endpoint-based bound | Same, safer tolerance | Engineering only |
| Demand duplicates | Could repeat shortest-path rows | Coalesced exactly | Exact reduction |
| Crossing semantics | Could planarize all crossings | Shared-source-vertex default | Preprocessing/model change |
| Shared interior junctions | Not fully recovered | Recovered | Preprocessing fix |
| Component selection | Could depend on node count | Uses total road length | Preprocessing policy |
| CRS policy | Weaker assumptions | Explicit metre-based handling | Measurement-policy change |
| Parallel edges | Duplicate entries could sum | Minimum direct edge retained | Correctness fix |
| Closed rings | Could disappear | Preserved | Correctness fix |
| Orphan nodes | Could survive filtering | Compacted | Correctness fix |
| Float32 | Limited support | Optional with safe tolerance | Numerical option |
| Parallelism | More aggressive experimentation | Conservative default | Execution-policy change |

---

## 13. Scope of v0.1.0

Implemented:

- undirected non-negative networks;
- weighted vertex 1-median;
- unweighted vertex 1-center;
- exact unweighted absolute 1-center;
- GIS linework to sparse network conversion;
- shared-source-vertex topology recovery;
- optional explicit planar noding;
- point/polygon demand input;
- nearest-node demand snapping;
- duplicate-demand coalescing;
- bounded-memory shortest-path calculation;
- optional float32 storage;
- optional parallel execution;
- GIS output.

Not implemented:

- directed roads;
- one-way restrictions;
- turn restrictions;
- asymmetric costs;
- time-dependent costs;
- weighted minimax center;
- continuous edge-based demand;
- `p > 1` center or median;
- automatic interpretation of all bridge/tunnel/layer semantics;
- out-of-core shortest-path matrices.

---

## 14. Final review conclusion

The final v0.1.0 release is best understood as a **classical network-location
solver with substantial GIS, numerical, and implementation hardening**.

The most important review lesson is that mathematical correctness and graph
correctness are different things.

The facility-location formulas are compact. The difficult practical work is
ensuring that:

- the intended junctions exist;
- false junctions do not;
- metric units are valid;
- graph cleaning does not silently delete roads;
- duplicate sparse edges do not corrupt direct distances;
- demand mapping is explicit;
- floating-point shortcuts do not invalidate pruning;
- memory optimizations preserve exact mathematical equivalence.

The public v0.1.0 release therefore retains the classical optimization
objectives while making the surrounding graph construction and numerical
behavior substantially more explicit and defensible.
