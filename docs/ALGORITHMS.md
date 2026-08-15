# Algorithms, and who they belong to

Every algorithm in this package was invented by someone else, most of them
decades ago. This document records what is used, who created it, and — kept
carefully separate at the end — which parts are implementation work done with
AI assistance rather than anything original.

Nothing here is presented as a new algorithm. The package implements the
classical absolute 1-centre directly because that continuous edge-interior case
was the gap motivating the original codebase.

## Notation

- `n` — junctions (graph vertices)
- `m` — road segments (graph edges)
- `k` — demand points; `k ≤ n`, and usually `k ≪ n`
- `d(x, y)` — shortest road distance
- `L` — the length of a particular segment

---

## Stage 1 — Turning a drawing into a network

### Recovering source topology without inventing crossings

A line layer records geometry, not necessarily road topology. Two lines may
cross on the page without a valid turn between them (for example, a bridge over
a road), while a genuine OSM-style junction can be stored as a coordinate in
the *interior* of one or both LineStrings.

The v0.3 default therefore uses a deliberately conservative rule: if a snapped
coordinate location is already present in at least two distinct input
LineStrings, it is treated as a shared source vertex and participating lines are
split there. This recovers T- and X-junctions already encoded by the source but
does not create a new coordinate at a mere geometric crossing.

This shared-vertex recovery is straightforward array bookkeeping rather than a
new graph algorithm; no novelty is claimed.

### Optional planar splitting of all geometric crossings

When the caller explicitly enables planar noding, every two-dimensional
crossing is treated as connected. That is appropriate only when the input is
genuinely planar.

The sweep-line method for reporting all intersections among a set of segments
is due to **Jon Bentley and Thomas Ottmann (1979)**. Producing output that
remains topologically consistent once coordinates are rounded to finite
precision — the practical difficulty — is **snap rounding**, introduced by
**Daniel Greene and F. Frances Yao (1986)** and made practical by
**John Hobby (1999)**.

This package does not implement those algorithms directly. Explicit planar
noding calls `shapely.union_all`, which binds to **GEOS**, itself a port of the
**JTS Topology Suite** written by **Martin Davis**.

### Merging away shape points

A curved road is drawn with many coordinates. Only the junctions matter for
routing; the rest describe shape. Merging chains of segments that meet
end-to-end at a point where nothing else joins can reduce `n` dramatically.
The distance matrix costs `O(kn)` memory and becomes `O(n²)` when every node is
used as demand, so unnecessary shape vertices can dominate the problem.
Implemented by `shapely.line_merge`, again GEOS/JTS (Martin Davis).

### Matching near-identical coordinates

Endpoints meant to coincide often differ in the last decimal places. They are
rounded to a tolerance and compared as integers. This is folklore rather than
anyone's algorithm.

### Finding the connected pieces

A network digitised from real data usually contains stranded fragments.
Identifying connected components in linear time is standard depth-first search,
in the form given by **John Hopcroft and Robert Tarjan (1973)**. Used via
`scipy.sparse.csgraph.connected_components`.

### Attaching demand points to the network

Demand points rarely sit exactly on a road, so each is matched to its nearest
junction using a **k-d tree**, invented by **Jon Bentley (1975)**. Used via
`scipy.spatial.cKDTree`.

---

## Stage 2 — Shortest paths

Distances come from **Dijkstra's algorithm**, published by
**Edsger W. Dijkstra in 1959** in a paper slightly over two pages long. It is
run once per demand point. Used via `scipy.sparse.csgraph.dijkstra`; not
reimplemented here.

Cost is `O(k · (m + n log n))` with a heap, and the result occupies
`O(k · n)` memory. On any network large enough to matter, the memory is what
stops you rather than the arithmetic.

---

## Stage 3 — The 1-median (minisum)

Find `x` minimising `Σ_v w_v · d(x, v)`.

The governing result is **Hakimi's vertex optimality theorem**, from
**S. L. Hakimi (1964)**: along any segment the objective is concave, so it
attains its minimum at an endpoint, and therefore an optimal median can always
be taken at a junction. Hakimi extended this to `p` medians in **(1965)**.

This collapses a continuous problem into a discrete one, and the implementation
is a single weighted sum over the distance matrix followed by an `argmin`.

Maximising closeness centrality — **Linton Freeman (1978)**, building on
**Alex Bavelas (1950)** — is the same computation with unit weights.

On a tree the 1-median is computable in linear time, shown by
**A. J. Goldman (1971)**.

---

## Stage 4 — The absolute 1-centre (minimax)

Find `x` anywhere on the network minimising `max_v d(x, v)`.

Vertex optimality **fails** here, and this is the whole difficulty. The optimum
can lie strictly inside a segment. Both the problem and the term *absolute
centre* are due to **Hakimi (1964)**, in the same paper as the median theorem.

The efficient algorithm is **Oded Kariv and S. L. Hakimi (1979)**, part I,
which solves the absolute `p`-centre and gives `O(mn log n)` for `p = 1`. The
same paper proves the `p`-centre problem NP-hard for general `p`; part II does
the same for the `p`-median. **Edward Minieka (1970)** treated the related
`m`-centre problem and the variant where every point of the network, not only
the junctions, counts as demand.

### The construction implemented here

For a segment from `u` to `w` of length `L`, and a position `t` metres from
`u`:

```
d(t, v) = min(a_v + t,  b_v + L − t)        a_v = d(u, v),  b_v = d(w, v)
```

a tent-shaped function of `t`. Eccentricity is the upper envelope of `k` tents:
piecewise linear with slopes ±1, not concave, so its minimum can be interior.

Tent `v` is on its rising branch exactly when `t ≤ t*_v = (b_v − a_v + L)/2`.
Sorting demands by `t*_v` makes the rising set a suffix and the falling set a
prefix of the sorted order, so between consecutive breakpoints

```
ecc(t) = max(A + t,  B + L − t)
```

with `A` a suffix maximum of `a` and `B` a prefix maximum of `b`, both constant
on the interval. That is a convex V, minimised where the two lines meet,
clamped to the interval. `O(k log k)` per segment, `O(mk log k)` overall.

The triangle inequality guarantees `t*_v ∈ [0, L]`, so no interval is degenerate.

This decomposition is the standard reading of Kariv and Hakimi's local-centre
computation, presented in textbook form by **Handler and Mirchandani (1979)**
and **Daskin (2013)**. **Nimrod Megiddo (1983)** later removed the logarithmic
factor using parametric search; that refinement is *not* implemented here,
because at these problem sizes the sort is not the bottleneck.

### Special cases in the literature, not implemented

On a tree the centre is far easier — **Camille Jordan (1869)** characterised it
as one vertex or two adjacent ones, in a paper about trees written a century
before anyone had a computer.

---

## Complexity summary

| Stage | Cost | Source |
|---|---|---|
| Noding | `O((s + i) log s)` for `s` segments, `i` crossings | Bentley & Ottmann 1979 |
| Shortest paths | `O(k(m + n log n))` | Dijkstra 1959 |
| 1-median | `O(kn)` | Hakimi 1964 |
| Jordan centre | `O(kn)` | — |
| Absolute centre | `O(mk log k)` worst case | Kariv & Hakimi 1979 |

In practice the absolute centre stage is far below its worst case, because of
the bound described in the next section.

---

## AI-assisted implementation work

The following were developed in a working session with Claude (Anthropic). They
are engineering on top of the published algorithms above — none is a new
algorithm, and none should be cited as a contribution to the literature. They
are listed because the person reviewing this code should know which lines came
from a textbook and which came from a conversation, and therefore where to look
hardest for mistakes.

### Optimisations

1. **Segment pruning by a lower bound.** From any point on a segment,
   `d(x, v) ≥ min(d(u, v), d(w, v))`, so `max_v min(d(u,v), d(w,v))` is a lower
   bound on that segment's eccentricity. Segments whose bound already exceeds
   the best known junction radius cannot contain a better centre and are not
   swept. This follows directly from the distance formula; no novelty is
   claimed.

2. **Batched vectorisation of the sweep.** The interval decomposition is
   Kariv–Hakimi's; expressing it as prefix/suffix running maxima over a
   `(segments × demands)` array so many segments are solved together is an
   implementation choice, as is the `-inf` padding that removes boundary
   special cases.

3. **Memory-budgeted sweep blocking.** Both the lower-bound calculation and the
   exact sweep process segments in cell-bounded blocks so temporary memory does
   not scale directly with the full edge count.

4. **Contiguous row-wise sort inputs.** Endpoint-distance columns are copied
   into row-contiguous arrays before the per-edge switch-point sort. This costs
   one controlled copy and makes the subsequent vectorised operations simpler
   and generally faster on large blocks.

5. **Conservative parallelism.** `n_jobs` is a ceiling, not a promise. Small
   jobs remain serial because worker startup can dominate useful work. v0.3
   defaults to one worker and exposes both process and thread backends for
   deliberate benchmarking rather than assuming one is universally superior.

6. **Bounded shortest-path chunks.** SciPy returns Dijkstra distances in
   float64. Source nodes are split into bounded row groups, each result is cast
   promptly to the requested storage dtype, and rows are written into the final
   preallocated matrix. This keeps a float32 request from requiring one full
   float64 matrix at the same time as the complete float32 matrix.

7. **Streaming parallel assembly.** Parallel Dijkstra blocks are consumed as an
   ordered joblib generator and written directly to final output rows. The code
   therefore avoids retaining all completed blocks and then constructing a
   second complete matrix with `numpy.vstack`.

8. **Exact demand coalescing.** Repeated observations snapped to the same node
   are collapsed before Dijkstra. Multiplicity is irrelevant to a maximum, so
   the centre is unchanged; counts or user weights are summed for the median,
   preserving the minisum objective exactly.

The centre pruning step also uses a dtype-and-scale-aware safety margin so
float32 rounding cannot discard a borderline candidate edge.

### An optimisation that was tried and rejected

**Best-first bound tightening** — sweeping the most promising segments first,
tightening the incumbent radius from the result, and re-pruning the remainder.
This is standard branch-and-bound practice and it was implemented and measured.
It was **removed** because the first lower-bound pass already made the second
pass too small to justify the added branching and code complexity on the tested
workloads. It is recorded here so future optimisation work starts from a
measurement rather than an assumption.

### Correctness bugs found and fixed

The failures below all had the potential to produce plausible wrong answers
rather than obvious crashes, which is the dangerous kind. They are kept here as
regression history.

1. **Parallel edges silently summed.** `coo_matrix(...).tocsr()` adds duplicate
   entries, so two roads joining the same pair of junctions became one road of
   their combined length. This made `d(u,w) > L`, breaking the triangle
   inequality the sweep depends on, which pushed breakpoints outside `[0, L]`,
   inverted intervals, and returned a radius **below** the true one. Fixed by
   collapsing duplicates with a minimum (`topology.adjacency_from_edges`), plus
   a defensive clamp in the sweep.

2. **Closed loops deleted.** A ring road merges into a single line whose ends
   coincide; the self-loop filter then removed it, deleting real road. On a test
   case a **1.9 km ring vanished**. Fixed by cutting such loops in half.

3. **Empty geometries corrupting endpoints.** A geometry with no coordinates
   left an empty run in the coordinate index, so the "last coordinate" lookup
   reached into the previous line's data and attached segments to the wrong
   junctions. Fixed by filtering empties before the lookup.

4. **Needless reprojection.** Always converting to the estimated local UTM zone
   is right for latitude/longitude input and wrong for already-projected
   metre-based input. Projected metre input is now left untouched, while
   explicit target CRSs are checked to ensure the documented metre semantics.

5. **Nearly closed rings deleted.** GIS rings are sometimes visually closed but
   have a final coordinate that differs from the first by floating-point dust.
   Exact `is_closed` testing can miss them; endpoint snapping then maps both ends
   to one node and the resulting self-loop is discarded. The hardened graph builder detects closure
   using the same endpoint tolerance used for node identity and splits the ring
   before that can happen.

6. **False junctions at grade-separated crossings.** Treating every geometric
   crossing as a routable junction can connect a bridge to the road below it.
   Planar noding is explicit rather than the default.

7. **Shared-vertex noding stopped at T-junctions.** An intermediate patch split
   a through-road when another line *ended* on one of its interior vertices.
   That still missed a genuine junction where two ways both continued through
   the same shared source vertex. v0.3 defines a junction by participation of
   at least two distinct input LineStrings, covering both T and X cases without
   inventing a node at an ordinary geometric crossing.

8. **Filtered slivers left orphan graph nodes.** Node IDs were assigned before
   tiny edges and snapped self-loops were removed. With component filtering
   disabled, a discarded edge could leave an isolated CSR row/column and later
   produce infinite distances in an otherwise connected surviving network.
   Node IDs are now compacted immediately after edge filtering.

Each has a regression test in `tests/test_graph.py`, `tests/test_distances.py`,
or `tests/test_center.py`.

### Verification

Correctness of the sweep is checked against brute force: twelve random networks
are sampled at 4,001 points along every segment, and the sweep's answer must be
no worse than the sampled minimum and within the sampling resolution of it.
Further tests confirm hand-computable cases, that the answer never exceeds the
Jordan centre, and that results are invariant to worker count and block size.

Brute force is a genuinely independent check — it shares no code with the sweep
— but it only verifies the sweep given a distance matrix. It does not verify
that the distance matrix describes the roads you meant. That remains the
reviewer's job, and stage 1 is where the real risk lives.

---

## References

Bavelas, A. (1950). Communication patterns in task-oriented groups.
*Journal of the Acoustical Society of America*, 22(6), 725–730.

Bentley, J. L. (1975). Multidimensional binary search trees used for
associative searching. *Communications of the ACM*, 18(9), 509–517.

Bentley, J. L., & Ottmann, T. A. (1979). Algorithms for reporting and counting
geometric intersections. *IEEE Transactions on Computers*, C-28(9), 643–647.

Daskin, M. S. (2013). *Network and Discrete Location: Models, Algorithms and
Applications* (2nd ed.). Wiley.

Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
*Numerische Mathematik*, 1, 269–271.

Freeman, L. C. (1978). Centrality in social networks: conceptual clarification.
*Social Networks*, 1(3), 215–239.

Goldman, A. J. (1971). Optimal center location in simple networks.
*Transportation Science*, 5(2), 212–221.

Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry.
*27th Annual Symposium on Foundations of Computer Science*, 143–152.

Hakimi, S. L. (1964). Optimum locations of switching centers and the absolute
centers and medians of a graph. *Operations Research*, 12(3), 450–459.

Hakimi, S. L. (1965). Optimum distribution of switching centers in a
communication network and some related graph theoretic problems.
*Operations Research*, 13(3), 462–475.

Handler, G. Y., & Mirchandani, P. B. (1979). *Location on Networks: Theory and
Algorithms*. MIT Press.

Hobby, J. D. (1999). Practical segment intersection with finite precision
output. *Computational Geometry*, 13(4), 199–214.

Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: efficient algorithms for
graph manipulation. *Communications of the ACM*, 16(6), 372–378.

Jordan, C. (1869). Sur les assemblages de lignes. *Journal für die reine und
angewandte Mathematik*, 70, 185–190.

Kariv, O., & Hakimi, S. L. (1979). An algorithmic approach to network location
problems. I: The p-centers. *SIAM Journal on Applied Mathematics*, 37(3),
513–538.

Kariv, O., & Hakimi, S. L. (1979). An algorithmic approach to network location
problems. II: The p-medians. *SIAM Journal on Applied Mathematics*, 37(3),
539–560.

Megiddo, N. (1983). Applying parallel computation algorithms in the design of
serial algorithms. *Journal of the ACM*, 30(4), 852–865.

Minieka, E. (1970). The m-center problem. *SIAM Review*, 12(1), 138–139.

Mirchandani, P. B., & Francis, R. L. (Eds.) (1990). *Discrete Location Theory*.
Wiley.

### Software

GEOS and the JTS Topology Suite (Martin Davis and contributors) — noding and
line merging, via Shapely.
SciPy (`sparse.csgraph`, `spatial`) — Dijkstra, connected components, k-d tree.
NumPy — array operations throughout.
joblib — optional threading/process parallel scheduling.
