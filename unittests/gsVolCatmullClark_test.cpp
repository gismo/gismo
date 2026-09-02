/** @file gsVolCatmullClark_test.cpp

    @brief Tests volumetric Catmull-Clark subdivision on gsVolMesh.

    The refinement replaces every (cell, corner) pair by one cell and rewires
    every dart, so almost any mistake shows up as a broken 3-map rather than as
    a wrong number.  Every case therefore ends with is_valid_topology() and with
    a check that no cell came out inverted -- an orientation slip in either of
    the two quad families would either be rejected by add_cell() or produce a
    negative volume.

    The counts are pinned as well, because they encode the topology of the
    scheme: a corner of valence k yields a cell with 2k faces, so a trivalent
    corner yields a hexahedron.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
**/

#include "gismo_unittest.h"

#include <gsMesh2/gsSubdivisionSchemes/gsVolCatmullClark.h>

SUITE(gsVolCatmullClark_test)
{

typedef gsVolMesh<real_t>          Mesh;
typedef gsVolMesh<real_t>::Point   Point;
typedef gsVolMesh<real_t>::Vertex  Vertex;
typedef gsVolMesh<real_t>::Cell    Cell;

inline index_t euler(const Mesh & m)
{
    return (index_t)m.n_vertices() - (index_t)m.n_edges()
         + (index_t)m.n_faces()    - (index_t)m.n_cells();
}

// sound 3-map, and no cell turned inside out
inline void checkSound(const Mesh & m)
{
    std::string msg;
    const bool ok = m.is_valid_topology(&msg);
    if (!ok) gsInfo << "is_valid_topology: " << msg << "\n";
    CHECK( ok );

    for (auto c : m.cells())
        CHECK( m.volume(c) > 0 );
}

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


// One hexahedron becomes the 2x2x2 block: eight cells over the 3x3x3 lattice
// of vertex, edge, face and cell points.
TEST(Hexahedron)
{
    Mesh m;
    unitHex(m);
    checkSound(m);

    gsVolCatmullClark<real_t>::apply(m);

    CHECK_EQUAL(8u,  m.n_cells());
    CHECK_EQUAL(27u, m.n_vertices());
    CHECK_EQUAL(54u, m.n_edges());
    CHECK_EQUAL(36u, m.n_faces());
    CHECK_EQUAL(1,   euler(m));
    CHECK( m.is_hex_mesh() );
    checkSound(m);

    gsVolCatmullClark<real_t>::apply(m);
    CHECK_EQUAL(64u,  m.n_cells());
    CHECK_EQUAL(125u, m.n_vertices());
    CHECK( m.is_hex_mesh() );
    checkSound(m);
}

// Every corner of a tetrahedron is trivalent, so one step already gives four
// hexahedra.
TEST(Tetrahedron_BecomesHexahedral)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    m.add_tet(v0,v1,v2,v3);

    gsVolCatmullClark<real_t>::apply(m);

    CHECK_EQUAL(4u, m.n_cells());
    CHECK_EQUAL(1,  euler(m));
    CHECK( m.is_hex_mesh() );
    checkSound(m);

    gsVolCatmullClark<real_t>::apply(m);
    CHECK_EQUAL(32u, m.n_cells());
    CHECK( m.is_hex_mesh() );
    checkSound(m);
}

// The apex of a pyramid has four incident faces, so it yields an 8-faced cell.
// That cell has a vertex point and a cell point of valence four again, so the
// defect is not repaired by refining further -- it persists as two poles.
TEST(Pyramid_ApexIsNeverHexahedral)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(1,1,0));
    const Vertex v3 = m.add_vertex(Point(0,1,0));
    const Vertex v4 = m.add_vertex(Point(0.5,0.5,1));
    m.add_pyramid(v0,v1,v2,v3,v4);

    gsVolCatmullClark<real_t>::apply(m);

    CHECK_EQUAL(5u, m.n_cells());          // one per corner of the pyramid
    CHECK( !m.is_hex_mesh() );             // the apex cell has eight faces
    CHECK_EQUAL(1, euler(m));
    checkSound(m);

    unsigned int eight = 0;
    for (auto c : m.cells()) if (8u == m.valence(c)) ++eight;
    CHECK_EQUAL(1u, eight);

    gsVolCatmullClark<real_t>::apply(m);
    CHECK( !m.is_hex_mesh() );             // still not repaired
    checkSound(m);
}

TEST(PrismIsTrivalentEverywhere)
{
    Mesh m;
    const Vertex v0 = m.add_vertex(Point(0,0,0));
    const Vertex v1 = m.add_vertex(Point(1,0,0));
    const Vertex v2 = m.add_vertex(Point(0,1,0));
    const Vertex v3 = m.add_vertex(Point(0,0,1));
    const Vertex v4 = m.add_vertex(Point(1,0,1));
    const Vertex v5 = m.add_vertex(Point(0,1,1));
    m.add_prism(v0,v1,v2,v3,v4,v5);

    gsVolCatmullClark<real_t>::apply(m);

    CHECK_EQUAL(6u, m.n_cells());
    CHECK( m.is_hex_mesh() );
    checkSound(m);
}

// A block with interior faces exercises the gluing: the refined cells have to
// share their faces, not duplicate them.
TEST(BlockStaysGlued)
{
    Mesh m;
    const index_t n = 2, s = n+1;
    std::vector<Vertex> id(s*s*s);
    for (index_t i=0;i!=s;++i) for (index_t j=0;j!=s;++j) for (index_t k=0;k!=s;++k)
        id[(i*s+j)*s+k] = m.add_vertex(Point(i,j,k));

#   define ID(a,b,c) id[(((a)*s+(b))*s+(c))]
    for (index_t i=0;i!=n;++i) for (index_t j=0;j!=n;++j) for (index_t k=0;k!=n;++k)
        m.add_hex(ID(i,j,k),   ID(i+1,j,k),   ID(i+1,j+1,k),   ID(i,j+1,k),
                  ID(i,j,k+1), ID(i+1,j,k+1), ID(i+1,j+1,k+1), ID(i,j+1,k+1));
#   undef ID

    gsVolCatmullClark<real_t>::apply(m);

    CHECK_EQUAL(64u,  m.n_cells());
    CHECK_EQUAL(125u, m.n_vertices());
    CHECK_EQUAL(300u, m.n_edges());
    CHECK_EQUAL(240u, m.n_faces());
    CHECK_EQUAL(1,    euler(m));
    CHECK( m.is_hex_mesh() );
    checkSound(m);

    // refining the 2x2x2 block gives a 4x4x4 one, whose boundary is six 4x4
    // grids of faces; anything that failed to glue would show up as extra
    // boundary faces here
    unsigned int bnd = 0;
    for (auto f : m.faces()) if (m.is_boundary(f)) ++bnd;
    CHECK_EQUAL(96u, bnd);
}

// Until the boundary masks are implemented the interior rules are applied
// everywhere, which shrinks the boundary.  Pin that down so the day it changes
// is a deliberate one.
TEST(BoundaryStillShrinks)
{
    Mesh m;
    unitHex(m);
    const real_t v0 = m.volume();

    gsVolCatmullClark<real_t>::apply(m);
    const real_t v1 = m.volume();

    CHECK( v1 < v0 );
    checkSound(m);
}

TEST(CheckMeshRejectsEmpty)
{
    Mesh empty;
    gsVolCatmullClark<real_t> cc;
    CHECK_EQUAL((int)gsVolSubdivisionScheme<real_t>::INVALID,
                (int)cc.check_mesh(empty));

    Mesh m;
    unitHex(m);
    CHECK_EQUAL((int)gsVolSubdivisionScheme<real_t>::VALID,
                (int)cc.check_mesh(m));
}

} // SUITE
