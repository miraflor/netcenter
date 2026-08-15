# Review notes for v0.3.0

This file records the second-pass review that produced v0.3.0 from the supplied
v0.2.1 candidate and its patch files.

## What was worth keeping

### Shared-source-vertex splitting

The v0.2.1 patch correctly identified a major weakness in strict endpoint-only
network construction: OSM-style road ways commonly carry real junctions as
vertices inside a LineString. A side street can therefore end on the interior
vertex of a through-road.

Recovering topology from source vertices is substantially safer than full
planar noding because it does not create a node where two lines merely cross in
2-D.

### Separate parallel backends

Shortest-path work and the absolute-centre sweep have different memory and
scheduling characteristics. Keeping separate backend settings is reasonable.
The revised code does not, however, claim that one backend is universally
faster.

## What needed correction

### 1. Endpoint-to-interior was not enough

The supplied shared-vertex patch looked only for locations used as an endpoint
of some line. That recovers many T-junctions but misses a real junction where
two continuing ways both contain the shared node as an interior vertex.

v0.3.0 defines a shared junction as a snapped coordinate location present in at
least two **distinct input LineStrings**. Every participating line is split at
that location when the vertex lies in its interior.

This still preserves bridge safety: two lines crossing geometrically without a
shared source vertex remain disconnected.

### 2. Filtered geometry could leave orphan nodes

Node IDs were assigned before snapped self-loops and tiny slivers were removed.
With `keep_largest_component=False`, a node used only by a discarded edge could
remain as an isolated CSR row/column. Later shortest-path calculations would
then contain infinities even though the surviving road component itself was
connected.

v0.3.0 compacts node IDs immediately after edge filtering.

### 3. Parallel shortest paths duplicated output memory

The parallel path returned all chunk arrays as a list and then called
`numpy.vstack`. Near completion that can hold the individual blocks *and* a
second complete stacked copy.

v0.3.0 preallocates the final matrix and consumes joblib results as an ordered
generator, writing each block directly into its final rows.

### 4. Release metadata had drifted

The supplied tree declared version `0.2.1` in `pyproject.toml` but exported
`0.2.0` from `netcenter.__init__`. The README, changelog, validation count, CLI
backend default, and performance narrative also described different releases.

The public version now comes from installed package metadata, and the release
documentation is aligned around v0.3.0.

### 5. Worker settings were silently coerced

`n_jobs=0` and negative values were silently converted to one worker. That can
hide user mistakes. The numerical and distance layers now require a positive
integer explicitly.

### 6. The GIL explanation was too categorical

The submitted patch included a machine-specific claim that SciPy Dijkstra
necessarily makes threads useless. That is too brittle for library-level
API documentation: performance depends on the exact wrapper path, SciPy build,
problem shape, process overhead, and memory pressure.

v0.3.0 therefore keeps `n_jobs=1` as the safe default and documents both
`loky` and `threading` as explicit options to benchmark on the target system.

## Additional hardening

- snap-grid integer conversion now checks for overflow before casting;
- `max_temp_mb` must be finite and positive;
- negative matrix dimensions are rejected by the estimator;
- threaded chunk execution has a regression test that must match serial output;
- package metadata has its own regression test.

## Result

The supplied v0.2.1 candidate passed 55 tests. The revised v0.3.0 tree passes
69 tests in the review environment.
