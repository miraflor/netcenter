# Network assumptions

These assumptions matter more than any command-line switch. Read them before
applying `netcenter` to a new road dataset.

## 1. The graph is undirected

Every retained road segment can be traversed in both directions. One-way roads,
turn restrictions, and asymmetric travel costs are not represented.

## 2. Source topology is preferred over geometric appearance

The default does **not** assume that every line crossing is a junction.
Instead, it treats a coordinate already present in two or more distinct input
LineStrings as a real shared vertex and splits participating lines there.

That recovers common OSM-style junctions while preserving grade separation when
a bridge and the road below merely cross geometrically and do not share a source
vertex.

## 3. Three topology modes are available

### Default: shared-source-vertex noding

```python
build_network(...)
```

Recommended for OSM-shaped or routing-quality linework. Real junctions already
encoded as shared vertices are recovered, including cases where both ways
continue through the junction.

### Strict endpoint-only mode

```python
build_network(..., split_shared_vertices=False)
```

Use only when the source is already split into separate segments at every true
junction. This can be faster on already-normalized data, but it will miss
interior shared vertices.

CLI equivalent:

```bash
netcenter roads.gpkg --no-shared-vertex-noding
```

### Full planar noding

```python
build_network(..., node=True)
```

This creates a junction at every 2-D crossing. Use only for genuinely planar
linework where every crossing is connected.

CLI equivalent:

```bash
netcenter roads.gpkg --node-crossings
```

Do not use this blindly on bridges, tunnels, or grade-separated roads.

## 4. A shared vertex is evidence, not an attribute model

`netcenter` does not inspect OSM `bridge`, `tunnel`, `layer`, turn-restriction,
or one-way attributes. It only sees line geometry at graph-construction time.

If the source contains a geometrically shared vertex that should *not* permit a
turn, clean or separate that topology before building the network.

## 5. Distances are metres

Every solver result, edge length, snap grid spacing, and offset `t` is interpreted
as metres.

- metre-based projected input is preserved;
- geographic input is reprojected to an inferred local UTM CRS;
- non-metre projected input is converted when possible;
- an explicit `target_crs` must be projected and metre-based.

For a very large study area, choose a projection suitable for the full extent
instead of relying on one inferred UTM zone.

## 6. Demand is snapped to network nodes

Demand is not projected to an arbitrary point along the nearest edge. On a
heavily simplified network, the nearest junction can be much farther away than
the nearest road itself.

Always inspect snap distances. In production work, set a defensible
`--max-snap` threshold.

## 7. Only one connected component is retained by default

The component with the greatest total road length is kept. This avoids infinite
shortest-path distances between disconnected pieces.

If disconnected components are substantively meaningful, solve them separately
rather than interpreting unreachable travel as one facility-location problem.

## 8. Demand is vertex demand

The centre objective is the maximum distance to the chosen demand nodes. It is
not the continuous-network centre in which every point on every road is itself
a demand location.
