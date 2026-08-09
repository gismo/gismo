/** @file gsVolMesh_test.cpp

    @brief Tests gsVolMesh and its non-templated topology base gsVolMeshTopology.

    The structure is a combinatorial 3-map, so almost every bug shows up as a
    broken invariant rather than as a wrong number.  Every case therefore ends
    with is_valid_topology(), which checks that beta3 is a fixed-point-free
    involution, that it reverses the half-face cycles, that the cell-local and
    the geometric tier agree on the vertices, edges and faces of every dart, and
    that all four rings close.

    Beyond that the emphasis is on the three things that are easy to get subtly
    wrong:

      - gluing.  Two cells that share a face have to end up with one geometric
        face and one geometric edge per shared edge, not two of each.

      - the second tier.  Corners, edge-uses and half-faces are per cell, so
        their counts are a multiple of the cell count while the geometric counts
        are not; mixing the two up still satisfies Euler's formula.

      - value semantics and garbage collection, which have to renumber both
        tiers consistently and re-bind every property handle.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsVolMesh_test)
{

typedef gsVolMesh<real_t>          Mesh;
typedef gsVolMesh<real_t>::Point   Point;
typedef gsVolMesh<real_t>::Vertex  Vertex;
typedef gsVolMesh<real_t>::Edge    Edge;
typedef gsVolMesh<real_t>::Face    MFace;
typedef gsVolMesh<real_t>::Cell    Cell;

// V - E + F - C, i.e. 1 for a mesh whose union is a topological ball
inline index_t euler(const Mesh & m)
{
    return (index_t)m.n_vertices() - (index_t)m.n_edges()
         + (index_t)m.n_faces()    - (index_t)m.n_cells();
}

inline void checkValid(const Mesh & m)
{
    std::string msg;
    const bool sound = m.is_valid_topology(&msg);
    if (!sound) gsInfo << "is_valid_topology: " << msg << "\n";
    CHECK( sound );
}

// the unit tetrahedron, in VTK_TETRA ordering
inline Cell unitTet(Mesh & m)
{
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    return m.add_tet(v0,v1,v2,v3);
}

// the unit cube, in VTK_HEXAHEDRON ordering
inline Cell unitHex(Mesh & m)
{
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(1,1,0));
    const Vertex v3 = m.add_vertex(Point(0,1,0));
    const Vertex v4 = m.add_vertex(Point(0,0,1));
    const Vertex v5 = m.add_vertex(Point(1,0,1));
    const Vertex v6 = m.add_vertex(Point(1,1,1));
    const Vertex v7 = m.add_vertex(Point(0,1,1));
    return m.add_hex(v0,v1,v2,v3,v4,v5,v6,v7);
}

// an n x n x n block of unit hexahedra; ids is filled with the (n+1)^3 vertices
inline void hexBlock(Mesh & m, int n, std::vector<Vertex> & ids)
{
    const int s = n+1;
    ids.assign(s*s*s, Vertex());
    for (int i=0; i!=s; ++i)
        for (int j=0; j!=s; ++j)
            for (int k=0; k!=s; ++k)
                ids[(i*s+j)*s+k] = m.add_vertex(Point(i,j,k));

#   define ID(a,b,c) ids[(((a)*s+(b))*s+(c))]
    for (int i=0; i!=n; ++i)
        for (int j=0; j!=n; ++j)
            for (int k=0; k!=n; ++k)
                m.add_hex(ID(i,j,k),     ID(i+1,j,k),     ID(i+1,j+1,k),     ID(i,j+1,k),
                          ID(i,j,k+1),   ID(i+1,j,k+1),   ID(i+1,j+1,k+1),   ID(i,j+1,k+1));
#   undef ID
}


//---------------------------------------------------------------------------
// The four standard cell types
//---------------------------------------------------------------------------

TEST(SingleTet)
{
    Mesh m;
    const Cell c = unitTet(m);
    CHECK( c.is_valid() );
    checkValid(m);

    CHECK_EQUAL(4u, m.n_vertices());
    CHECK_EQUAL(6u, m.n_edges());
    CHECK_EQUAL(4u, m.n_faces());
    CHECK_EQUAL(1u, m.n_cells());
    CHECK_EQUAL(1,  euler(m));

    CHECK( m.is_tet(c) );
    CHECK( m.is_tet_mesh() );
    CHECK_EQUAL(4u, m.valence(c));

    // every face, edge and vertex of a lone cell is on the boundary
    for (auto f : m.faces()) { CHECK( m.is_boundary(f) ); CHECK_EQUAL(1u, m.valence(f)); }
    for (auto e : m.edges()) { CHECK( m.is_boundary(e) ); CHECK_EQUAL(1u, m.valence(e)); }
    for (auto v : m.vertices()) { CHECK( m.is_boundary(v) ); CHECK_EQUAL(1u, m.valence(v)); }

    // beta3 is empty, so no half-face has an opposite
    for (auto hf : m.halffaces()) CHECK( !m.opposite_halfface(hf).is_valid() );

    CHECK_CLOSE(1.0/6.0, m.volume(c), 1e-12);
}

// The second tier is per cell: a lone tetrahedron has as many corners as
// vertices, but four half-faces for four faces and twelve darts for six edges.
TEST(SingleTet_CellLocalTier)
{
    Mesh m;
    unitTet(m);

    CHECK_EQUAL(4u,  m.n_corners());
    CHECK_EQUAL(6u,  m.n_edge_uses());
    CHECK_EQUAL(4u,  m.n_halffaces());
    CHECK_EQUAL(12u, m.n_halfedges());

    // each corner is a use of exactly one vertex by exactly one cell
    for (auto cn : m.corners())
    {
        CHECK( m.vertex(cn).is_valid() );
        CHECK( m.cell(cn).is_valid() );
    }
}

TEST(SingleHex)
{
    Mesh m;
    const Cell c = unitHex(m);
    CHECK( c.is_valid() );
    checkValid(m);

    CHECK_EQUAL(8u,  m.n_vertices());
    CHECK_EQUAL(12u, m.n_edges());
    CHECK_EQUAL(6u,  m.n_faces());
    CHECK_EQUAL(1u,  m.n_cells());
    CHECK_EQUAL(8u,  m.n_corners());
    CHECK_EQUAL(12u, m.n_edge_uses());
    CHECK_EQUAL(6u,  m.n_halffaces());
    CHECK_EQUAL(24u, m.n_halfedges());

    CHECK( m.is_hex(c) );
    CHECK( m.is_hex_mesh() );
    CHECK_CLOSE(1.0, m.volume(c), 1e-12);

    // half-faces are oriented outwards, so every normal points away from the
    // barycenter of the cell
    const Point b = m.barycenter(c);
    for (auto hf : m.halffaces())
        CHECK( m.normal(hf).dot(m.barycenter(hf) - b) > 0 );
}

TEST(SinglePrismAndPyramid)
{
    {
        Mesh m;
        const Vertex v0 = m.add_vertex(Point(0,0,0));
        const Vertex v1 = m.add_vertex(Point(1,0,0));
        const Vertex v2 = m.add_vertex(Point(0,1,0));
        const Vertex v3 = m.add_vertex(Point(0,0,1));
        const Vertex v4 = m.add_vertex(Point(1,0,1));
        const Vertex v5 = m.add_vertex(Point(0,1,1));
        const Cell c = m.add_prism(v0,v1,v2,v3,v4,v5);
        CHECK( c.is_valid() );
        checkValid(m);
        CHECK_EQUAL(6u, m.n_vertices());
        CHECK_EQUAL(9u, m.n_edges());
        CHECK_EQUAL(5u, m.n_faces());
        CHECK( m.is_prism(c) );
        CHECK_CLOSE(0.5, m.volume(c), 1e-12);
    }
    {
        Mesh m;
        const Vertex v0 = m.add_vertex(Point(0,0,0));
        const Vertex v1 = m.add_vertex(Point(1,0,0));
        const Vertex v2 = m.add_vertex(Point(1,1,0));
        const Vertex v3 = m.add_vertex(Point(0,1,0));
        const Vertex v4 = m.add_vertex(Point(0.5,0.5,1));
        const Cell c = m.add_pyramid(v0,v1,v2,v3,v4);
        CHECK( c.is_valid() );
        checkValid(m);
        CHECK_EQUAL(5u, m.n_vertices());
        CHECK_EQUAL(8u, m.n_edges());
        CHECK_EQUAL(5u, m.n_faces());
        CHECK( m.is_pyramid(c) );
        CHECK_CLOSE(1.0/3.0, m.volume(c), 1e-12);
    }
}

// A cell whose boundary is not a closed orientable surface has to be refused
// before anything is written to the mesh.
TEST(AddCell_RejectsBadInput)
{
    Mesh m;
    std::vector<Vertex> v;
    for (int i=0; i!=4; ++i) v.push_back(m.add_vertex(Point(i,0,0)));

    // three faces cannot bound a cell
    std::vector< std::vector<Vertex> > tooFew(3, std::vector<Vertex>(3, v[0]));
    CHECK( !m.add_cell(tooFew).is_valid() );

    // a tetrahedron with one face flipped is not consistently oriented
    std::vector< std::vector<Vertex> > flipped(4);
    const int f[4][3] = { {0,2,1}, {0,1,3}, {1,2,3}, {0,2,3} }; // last one reversed
    for (int i=0; i!=4; ++i)
        for (int j=0; j!=3; ++j) flipped[i].push_back(v[f[i][j]]);
    CHECK( !m.add_cell(flipped).is_valid() );

    // nothing was added
    CHECK_EQUAL(0u, m.n_cells());
    CHECK_EQUAL(0u, m.n_faces());
    CHECK_EQUAL(0u, m.n_corners());
}


//---------------------------------------------------------------------------
// Gluing: beta3 and the shared geometry
//---------------------------------------------------------------------------

TEST(TwoTets_Glued)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    const Vertex v4 = m.add_vertex(Point(0,0,-1));

    const Cell c0 = m.add_tet(v0,v1,v2,v3);
    const Cell c1 = m.add_tet(v0,v2,v1,v4);   // shares the face (v0,v1,v2)
    CHECK( c0.is_valid() );
    CHECK( c1.is_valid() );
    checkValid(m);

    // the shared face and the three shared edges exist once, not twice
    CHECK_EQUAL(5u, m.n_vertices());
    CHECK_EQUAL(9u, m.n_edges());
    CHECK_EQUAL(7u, m.n_faces());
    CHECK_EQUAL(2u, m.n_cells());
    CHECK_EQUAL(1,  euler(m));

    // ... while the cell-local tier still counts one copy per cell
    CHECK_EQUAL(8u,  m.n_corners());
    CHECK_EQUAL(12u, m.n_edge_uses());
    CHECK_EQUAL(8u,  m.n_halffaces());

    unsigned int interior = 0;
    for (auto f : m.faces()) if (!m.is_boundary(f)) ++interior;
    CHECK_EQUAL(1u, interior);

    // beta3 is a fixed-point-free involution on the glued face
    unsigned int mated = 0;
    for (auto h : m.halfedges())
        if (m.mate(h).is_valid())
        {
            CHECK_EQUAL(h, m.mate(m.mate(h)));
            CHECK( m.mate(h) != h );
            ++mated;
        }
    CHECK_EQUAL(6u, mated);   // two triangles, three darts each

    CHECK_EQUAL(2u, m.valence(m.find_edge(v0,v1)));
    CHECK_EQUAL(1u, m.valence(m.find_edge(v0,v3)));
    CHECK_EQUAL(2u, m.valence(v0));
    CHECK_EQUAL(1u, m.valence(v3));

    // the two cells see each other through the shared face
    unsigned int nb = 0;
    for (auto n : m.cells(c0)) { CHECK_EQUAL(c1, n); ++nb; }
    CHECK_EQUAL(1u, nb);

    CHECK_CLOSE(1.0/3.0, m.volume(), 1e-12);
}

TEST(Sew_Unsew_RoundTrip)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    const Vertex v4 = m.add_vertex(Point(0,0,-1));

    m.add_tet(v0,v1,v2,v3);
    m.add_tet(v0,v2,v1,v4);

    const unsigned int nf = m.n_faces();

    // find the single interior face and take it apart
    gsVolMesh<real_t>::Halfface hf;
    for (auto f : m.faces())
        if (!m.is_boundary(f)) { hf = m.halfface(f,0); break; }
    CHECK( hf.is_valid() );

    const gsVolMesh<real_t>::Halfface other = m.opposite_halfface(hf);
    m.unsew(hf);
    checkValid(m);
    CHECK( m.is_boundary(hf) );
    CHECK( m.is_boundary(other) );
    CHECK_EQUAL(nf+1, m.n_faces());

    m.sew(hf, other);
    checkValid(m);
    CHECK( !m.is_boundary(hf) );
    CHECK_EQUAL(other, m.opposite_halfface(hf));
    // the extra face created by unsew() is marked deleted
    CHECK_EQUAL(nf, m.n_faces());
}


//---------------------------------------------------------------------------
// A block of hexahedra: the interesting valences live in its interior
//---------------------------------------------------------------------------

TEST(HexBlock_Counts)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);
    checkValid(m);

    CHECK_EQUAL(27u, m.n_vertices());
    CHECK_EQUAL(54u, m.n_edges());
    CHECK_EQUAL(36u, m.n_faces());
    CHECK_EQUAL(8u,  m.n_cells());
    CHECK_EQUAL(1,   euler(m));

    CHECK_EQUAL(64u,  m.n_corners());
    CHECK_EQUAL(96u,  m.n_edge_uses());
    CHECK_EQUAL(48u,  m.n_halffaces());
    CHECK_EQUAL(192u, m.n_halfedges());

    unsigned int bf = 0;
    for (auto f : m.faces()) if (m.is_boundary(f)) ++bf;
    CHECK_EQUAL(24u, bf);

    CHECK_CLOSE(8.0, m.volume(), 1e-12);
}

TEST(HexBlock_Valences)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);

    const int s = 3;
#   define ID(a,b,c) ids[(((a)*s+(b))*s+(c))]
    const Vertex centre = ID(1,1,1);
    const Vertex corner = ID(0,0,0);

    CHECK_EQUAL(8u, m.valence(centre));      // eight cells meet in the centre
    CHECK_EQUAL(1u, m.valence(corner));
    CHECK( !m.is_boundary(centre) );
    CHECK(  m.is_boundary(corner) );
    CHECK(  m.is_manifold(centre) );

    // an interior edge is surrounded by four cells, a corner edge by one
    const Edge interior = m.find_edge(ID(1,1,0), ID(1,1,1));
    CHECK( interior.is_valid() );
    CHECK_EQUAL(4u, m.valence(interior));
    CHECK( !m.is_boundary(interior) );
    CHECK( m.is_manifold(interior) );

    const Edge onEdge = m.find_edge(ID(0,0,0), ID(1,0,0));
    CHECK( onEdge.is_valid() );
    CHECK_EQUAL(1u, m.valence(onEdge));
    CHECK( m.is_boundary(onEdge) );
#   undef ID
}

// Each circulator has to visit exactly the entities the corresponding collect
// function returns; a ring that closes early or late shows up here.
TEST(HexBlock_CirculatorsAgree)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);

    for (auto c : m.cells())
    {
        CHECK_EQUAL(8u,  (unsigned)m.vertices(c).size());
        CHECK_EQUAL(12u, (unsigned)m.edges(c).size());
        CHECK_EQUAL(6u,  (unsigned)m.faces(c).size());
        CHECK_EQUAL(6u,  m.valence(c));
        CHECK_EQUAL(8u,  m.n_vertices(c));
        CHECK_EQUAL(12u, m.n_edges(c));
    }

    const int s = 3;
    const Vertex centre = ids[((1*s+1)*s)+1];
    CHECK_EQUAL(6u,  (unsigned)m.edges(centre).size());
    CHECK_EQUAL(12u, (unsigned)m.faces(centre).size());

    // one corner per incident cell, and one dart per (corner, incident face)
    unsigned int ncorners = 0;
    for (auto cn : m.corners(centre)) { ++ncorners; GISMO_UNUSED(cn); }
    CHECK_EQUAL(8u, ncorners);

    unsigned int ndarts = 0;
    for (auto h : m.halfedges(centre)) { ++ndarts; GISMO_UNUSED(h); }
    CHECK_EQUAL(24u, ndarts);   // eight cells, three outgoing darts each

    // the radial cycle of an interior edge closes on the cells around it
    const Edge interior = m.find_edge(ids[((1*s+1)*s)+0], centre);
    unsigned int ncells = 0;
    for (auto cc : m.cells(interior)) { ++ncells; GISMO_UNUSED(cc); }
    CHECK_EQUAL(4u, ncells);

    unsigned int nhf = 0;
    for (auto hf : m.halffaces(interior)) { ++nhf; GISMO_UNUSED(hf); }
    CHECK_EQUAL(4u, nhf);
}


//---------------------------------------------------------------------------
// Value semantics: the copy has to carry both tiers and re-bind its handles
//---------------------------------------------------------------------------

TEST(CopyIsDeepAndRebinds)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);

    Mesh cp(m);
    checkValid(cp);
    CHECK_EQUAL(m.n_vertices(),  cp.n_vertices());
    CHECK_EQUAL(m.n_edges(),     cp.n_edges());
    CHECK_EQUAL(m.n_faces(),     cp.n_faces());
    CHECK_EQUAL(m.n_cells(),     cp.n_cells());
    CHECK_EQUAL(m.n_halfedges(), cp.n_halfedges());

    // the geometry came along ...
    for (auto v : m.vertices())
        CHECK_CLOSE(0.0, (m.position(v) - cp.position(v)).norm(), 1e-15);

    // ... and writing to the copy leaves the original alone
    cp.position(Vertex(0)) = Point(42,42,42);
    cp.add_vertex(Point(7,7,7));
    CHECK_EQUAL(27u, m.n_vertices());
    CHECK_EQUAL(28u, cp.n_vertices());
    CHECK_CLOSE(0.0, m.position(Vertex(0)).norm(), 1e-15);

    Mesh as;
    as.assign(m);
    checkValid(as);
    CHECK_EQUAL(m.n_cells(), as.n_cells());
    for (auto v : m.vertices())
        CHECK_CLOSE(0.0, (m.position(v) - as.position(v)).norm(), 1e-15);
}

TEST(Properties)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 1, ids);

    gsVolMesh<real_t>::Cell_property<int> tag = m.add_cell_property<int>("C:tag", 0);
    for (auto c : m.cells()) tag[c] = 7;
    CHECK_EQUAL(7, tag[Cell(0)]);

    gsVolMesh<real_t>::Face_property<real_t> flux =
        m.add_face_property<real_t>("F:flux", 0.0);
    for (auto f : m.faces()) flux[f] = 1.5;
    CHECK_CLOSE(1.5, flux[MFace(0)], 1e-15);

    // the cell-local tier has its own namespace, so the same suffix is free
    gsVolMesh<real_t>::Halfface_property<int> hftag =
        m.add_halfface_property<int>("hf:tag", 3);
    CHECK_EQUAL(3, hftag[m.halfface(Cell(0))]);
}


//---------------------------------------------------------------------------
// Deletion and garbage collection have to renumber both tiers together
//---------------------------------------------------------------------------

TEST(DeleteCellAndCollect)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    const Vertex v4 = m.add_vertex(Point(0,0,-1));

    const Cell c0 = m.add_tet(v0,v1,v2,v3);
    m.add_tet(v0,v2,v1,v4);

    m.delete_cell(c0);
    CHECK_EQUAL(1u, m.n_cells());
    CHECK( m.garbage() );

    m.garbage_collection();
    CHECK( !m.garbage() );
    checkValid(m);

    // what is left is a lone tetrahedron
    CHECK_EQUAL(1u, m.n_cells());
    CHECK_EQUAL(4u, m.n_faces());
    CHECK_EQUAL(6u, m.n_edges());
    CHECK_EQUAL(4u, m.n_corners());
    CHECK_EQUAL(6u, m.n_edge_uses());
    CHECK_EQUAL(4u, m.n_halffaces());
    for (auto f : m.faces()) CHECK( m.is_boundary(f) );

    CHECK_CLOSE(1.0/6.0, m.volume(), 1e-12);
}

TEST(GarbageCollectionKeepsGeometry)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);

    // remove the eighth cell and compact
    Cell last;
    for (auto c : m.cells()) last = c;
    m.delete_cell(last);
    m.garbage_collection();
    checkValid(m);

    CHECK_EQUAL(7u, m.n_cells());
    CHECK_CLOSE(7.0, m.volume(), 1e-12);

    // every surviving position is still a lattice point of the block
    for (auto v : m.vertices())
    {
        const Point & p = m.position(v);
        for (index_t i = 0; i != 3; ++i)
            CHECK_CLOSE(p[i], math::round(p[i]), 1e-15);
    }
}


//---------------------------------------------------------------------------
// The boundary of the volume mesh, as a gsSurfMesh
//---------------------------------------------------------------------------

TEST(BoundaryMesh)
{
    Mesh m;
    std::vector<Vertex> ids;
    hexBlock(m, 2, ids);

    gsSurfMesh<real_t> b = m.boundary_mesh();

    CHECK_EQUAL(26u, b.n_vertices());   // 27 minus the interior one
    CHECK_EQUAL(24u, b.n_faces());
    CHECK_EQUAL(48u, b.n_edges());
    CHECK_EQUAL(2, (index_t)b.n_vertices() - (index_t)b.n_edges() + (index_t)b.n_faces());
    CHECK( b.is_quad_mesh() );

    // it is closed, and the back-map points at the right vertices
    for (auto e : b.edges()) CHECK( !b.is_boundary(e) );

    gsSurfMesh<real_t>::Vertex_property<index_t> back =
        b.get_vertex_property<index_t>("v:volvertex");
    for (auto v : b.vertices())
    {
        CHECK( back[v] >= 0 );
        CHECK_CLOSE(0.0, (b.position(v) - m.position(Vertex(back[v]))).norm(), 1e-15);
    }
}

TEST(CellMesh)
{
    Mesh m;
    const Cell c = unitHex(m);

    gsSurfMesh<real_t> b = m.cell_mesh(c);
    CHECK_EQUAL(8u,  b.n_vertices());
    CHECK_EQUAL(12u, b.n_edges());
    CHECK_EQUAL(6u,  b.n_faces());
    for (auto e : b.edges()) CHECK( !b.is_boundary(e) );
}

} // SUITE
