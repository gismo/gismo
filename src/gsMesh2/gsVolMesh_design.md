# The half-face structure of `gsMesh2`

Design note for `gsMeshTopology`, `gsVolMeshTopology` and `gsVolMesh`. The API itself is
documented in the doxygen comments; what follows is the reasoning behind the shape of the
data structure, and how it maps onto the literature.

## 1. Why `gsSurfMeshTopology` could not simply be extended

`gsSurfMeshTopology` is a half-edge structure whose invariants are index arithmetic:

```cpp
Edge     edge(Halfedge h)             const { return Edge(h.idx() >> 1); }
Halfedge halfedge(Edge e, unsigned i) const { return Halfedge((e.idx() << 1) + i); }
Halfedge opposite_halfedge(Halfedge h) const { return Halfedge(h.idx() ^ 1); }
```

That is the statement *an edge has exactly two halfedges*. In a volume mesh a geometric
edge is surrounded by *k* faces and therefore carries 2·*k* darts, so it cannot be a
`gsSurfMeshTopology::Edge`. The same arithmetic is, however, exactly right one level down:
inside a single cell, an edge is shared by exactly two of that cell's faces.

That observation is the whole design. It is also the one the references make:

* Feng, Wang, Weng & Tong, *Compact combinatorial maps: a volume mesh data structure*,
  Graphical Models **75** (2013) 149–156, §3.2.1: “Each 3-cell is treated locally as a
  2-manifold cell complex, which can be represented by a local half-edge structure, i.e. a
  2D combinatorial map.”
* CGoGN\_3 declares `struct CMap3 : public CMap2` and distinguishes `Vertex2`, `Edge2`,
  `Face2` (the cell-local cells) from `Vertex`, `Edge`, `Face`, `Volume`.

## 2. Class layout

```
gsMeshTopology                     oriented 2-map with property arrays: handles,
   |                               connectivity, iterators, circulators, add_face(),
   |                               deletion, garbage_collection()
   +-- gsSurfMeshTopology          the surface reading, plus the Euler operators
   |      |                        (triangulate, collapse, flip, split, insert_*)
   |      +-- gsSurfMesh<Scalar>   vertex positions and surface geometry
   |
   +-- gsVolMeshTopology           the 3-map: beta3, Cell, and the geometric 0/1/2-cells
          |
          +-- gsVolMesh<Scalar>    vertex positions and volume geometry
```

`gsMeshTopology` was factored out of `gsSurfMeshTopology` by a pure move; the latter's
public API is unchanged. The split is not cosmetic: the Euler operators of a surface would
silently violate the 3-map invariants, and because they live in the *sibling* class they
are structurally unreachable from `gsVolMeshTopology` rather than merely discouraged.

The inherited 2-map instance of `gsVolMeshTopology` holds the **cell-boundary complex**:
the disjoint union of the oriented boundary surfaces of all cells. Every cell boundary is
a closed oriented genus-0 surface and hence a legal component of the base class.

## 3. Entities

| tier | `gsVolMeshTopology` | is | CGoGN\_3 | per hexahedron |
|---|---|---|---|---|
| geometric | `Vertex` | a point of the mesh | `Vertex` | — |
| geometric | `Edge` | an edge shared by all incident cells | `Edge` | — |
| geometric | `Face` | a face shared by one or two cells | `Face` | — |
| volumetric | `Cell` | a polyhedral cell | `Volume` | 1 |
| cell-local | `Corner` = `gsMeshTopology::Vertex` | one use of a vertex by one cell | `Vertex2` | 8 |
| cell-local | `Edge_use` = `gsMeshTopology::Edge` | one use of an edge by one cell | `Edge2` | 12 |
| cell-local | `Halfedge` = `gsMeshTopology::Halfedge` | a **dart** | `HalfEdge` | 24 |
| cell-local | `Halfface` = `gsMeshTopology::Face` | an oriented face, seen from one cell | `Face2` | 6 |

The geometric names hide the inherited ones, so `position(Vertex)` and `vertices(Cell)`
read naturally; the cell-local tier is reached through the four aliases. The geometric
handles use a small CRTP base with **typed** comparison operators, because
`gsMeshTopology::Base_handle::operator==` accepts any handle and, with two tiers of
0/1/2-cells around, `Corner(3) == Vertex(3)` would otherwise compile and mean nothing.

The corner array is the paper's `Cv2V` table — the per-cell vertex list — in disguise.

## 4. The one new relation

| | operator | source |
|---|---|---|
| β₁ | `next_halfedge` / `prev_halfedge` | inherited |
| β₂ | `opposite_halfedge` (`h^1`), the dart on the adjacent face of the **same** cell | inherited |
| β₃ | `mate`, the dart on the same geometric edge in the opposite half-face | **new**: `Halfedge_property<Halfedge> hmate_` |

β₂ and β₃ each reverse the orientation of a dart, so their composition preserves it:

```
radial_next(h) = opposite_halfedge(mate(h))     // rotate around the edge, +1 cell
radial_prev(h) = mate(opposite_halfedge(h))
```

The radial orbit therefore visits **one dart per incident cell**, never repeating a cell.
It is a cycle for an interior edge and an open chain for a boundary edge; the circulator
rewinds to the free end in its constructor.

β₃ is partial — a dart with no mate is a boundary dart — mirroring the surface mesh, where
a halfedge with no face is a boundary halfedge. Materialising boundary cells to make β₃
total (CGoGN's `close()`, the paper's `B2D` table) is a possible later addition; it would
make every orbit a closed cycle at the price of maintaining the boundary incrementally.

## 5. Rings

Four singly linked circular rings, one `int` per entity:

* per `Cell`: its half-faces (`hfnext_`), its corners (`ccnext_`), its edge-uses (`ucnext_`);
* per `Vertex`: its corners (`cvnext_`).

With them `vertices(Cell)`, `edges(Cell)`, `faces(Cell)`, `cells(Vertex)`, `cells(Edge)`,
`vertices(Face)` and `edges(Face)` are genuine cycles and need no visited marks — CGoGN
needs a `DartMarker` for these. Only `edges(Vertex)` and `faces(Vertex)` require
de-duplication, and those two are deliberately `std::vector`-returning collect functions so
the cost is visible at the call site.

## 6. `add_cell`

`add_cell` takes the boundary faces as outward-oriented vertex loops and

1. **validates first**: every directed vertex pair must occur exactly once and its reverse
   exactly once (closed, orientable, manifold), and the boundary must have Euler
   characteristic 2. Nothing is written to the mesh before this passes, so the operation
   cannot fail half way;
2. allocates one `Corner` per distinct vertex and links the rings;
3. builds the half-faces directly rather than through `gsMeshTopology::add_face()`. The
   corners are private to the cell and the input is known to be a closed manifold surface,
   so every cell-local edge is used by exactly two loops and the wiring is unambiguous —
   there is no boundary to re-link and no complex vertex to guard against;
4. matches each cell-local edge against the existing geometry with `find_edge`, and each
   half-face with `find_free_halfface`, sewing dart by dart where a neighbour is found.

The searches walk the corner ring of one vertex and are O(vertex degree). No persistent
hash map is kept, so the result stays correct across garbage collection. `add_tet`,
`add_hex`, `add_prism` and `add_pyramid` are thin wrappers over static face templates in
VTK/CGNS ordering, in the spirit of the cell types of the paper's §3.2.1.

## 7. Invariants

`is_valid_topology(std::string*)` checks, and every unit test calls it:

1. β₃ is a fixed-point-free involution on the non-boundary darts;
2. β₃ reverses the half-face cycle: `mate(next(h)) == prev(mate(h))`, i.e. (β₁∘β₃)² = id;
3. a dart and its mate agree on their two vertices, on their geometric edge and on their
   geometric face;
4. a geometric face carries one or two half-faces, never more;
5. the four rings close and are consistent with the entities they list;
6. every cell boundary has Euler characteristic 2.

## 8. Garbage collection

`gsMeshTopology::garbage_collection()` compacts by disjoint transpositions, so the
resulting permutation is an involution and a single map serves both directions. It
destroys its index maps before returning, so a new protected hook was added:

```cpp
virtual void remap_handles(const Vertex_property<Vertex>&,
                           const Halfedge_property<Halfedge>&,
                           const Face_property<Face>&);
```

called while the maps are still alive. `gsVolMeshTopology::garbage_collection()` compacts
the four geometric arrays first, rewrites the cell-local → geometric maps, then delegates
to the base, whose callback fixes β₃, the four rings and the geometric → cell-local
pointers.

## 9. Cost

Storage was chosen to match the rest of `gsMesh2`: explicit per-dart connectivity in
property arrays, arbitrary polyhedra, cheap local edits. Measured resident-set growth
(release build, `int` indices, `double` coordinates, includes allocator slack):

| mesh | cells | darts | bytes/cell |
|---|---|---|---|
| 20×20×20 hexahedra | 8 000 | 192 000 | ≈ 1 850 |
| the same block, Kuhn-split into tetrahedra | 48 000 | 576 000 | ≈ 534 |

For comparison, Table 2 of the paper measures ≈ 465 bytes per tetrahedron for
OpenVolumeMesh, so this layout is in the same class. The paper's own compact scheme — dart
id = (cell, local dart), β₁/β₂ from per-cell-type lookup tables, only β₃ stored — reaches
≈ 8 integers per tetrahedron, roughly twelve times less, at the price of restricting cells
to declared types and making local topological edits hard. That trade was considered and
declined; the public API (β₁/β₂/β₃ accessors and orbit traversals, no index arithmetic
exposed) leaves room for such a back-end later.

## 10. Not done yet

* Reading volume meshes. `gsVolMesh::write` emits VTK `.vtu` (as `VTK_POLYHEDRON`, with the
  face stream) and delegates surface formats to `boundary_mesh()`; there is no reader yet.
* `close()` / boundary cells, i.e. a total β₃.
* Volumetric subdivision and trivariate patch extraction, which are the reason the
  structure exists.
