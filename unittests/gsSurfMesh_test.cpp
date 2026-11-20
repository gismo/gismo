/** @file gsSurfMesh_test.cpp

    @brief Tests gsSurfMesh and its non-templated topology base gsSurfMeshTopology.

    The emphasis is on the two things that regressed silently in the past:

      - geometry *values*, not just element counts.  A reader that fills the
        connectivity correctly but writes garbage coordinates passes every
        count-based check, so the vertex positions are asserted explicitly.

      - value semantics across the base/derived boundary.  Copying a mesh has
        to carry the topology (gsSurfMeshTopology) *and* the point property
        (gsSurfMesh), and must re-bind the property handles to the copy;
        forgetting either half aliases or truncates the mesh silently.

      - compute_face_normal()'s triangle branch, whose cross product operands
        must both be relative to the same vertex: a normal has no magnitude to
        look wrong and no obviously invalid value, so translation invariance
        and agreement with the general-polygon branch are pinned directly.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsSurfMesh_test)
{

typedef gsSurfMesh<real_t>          Mesh;
typedef gsSurfMesh<real_t>::Point   Point;
typedef gsSurfMesh<real_t>::Normal  Normal;
typedef gsSurfMesh<real_t>::Vertex  Vertex;

// V - E + F, i.e. 2 for a closed genus-0 surface, 2-2g in general
inline index_t euler(const Mesh & m)
{ return (index_t)m.n_vertices() - (index_t)m.n_edges() + (index_t)m.n_faces(); }

inline void boundingBox(const Mesh & m, Point & lo, Point & hi)
{
    lo = m.position(Vertex(0));
    hi = lo;
    for (auto v : m.vertices())
    {
        lo = lo.cwiseMin(m.position(v));
        hi = hi.cwiseMax(m.position(v));
    }
}

inline real_t totalEdgeLength(const Mesh & m)
{
    real_t s = 0;
    for (auto e : m.edges()) s += m.edge_length(e);
    return s;
}

inline bool readMesh(Mesh & m, const std::string & relPath)
{
    const std::string p = gsFileManager::findInDataDir(relPath);
    return m.read(p);
}

/// One triangle with the given corner positions, in the given order.
Mesh triangleMesh(const Point & a, const Point & b, const Point & c)
{
    Mesh m;
    std::vector<Vertex> v;
    v.push_back(m.add_vertex(a));
    v.push_back(m.add_vertex(b));
    v.push_back(m.add_vertex(c));
    m.add_face(v);
    return m;
}


//---------------------------------------------------------------------------
// Reading: counts *and* coordinates
//---------------------------------------------------------------------------

TEST(ReadOff_Counts)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    CHECK_EQUAL(8u, m.n_vertices());
    CHECK_EQUAL(12u, m.n_edges());
    CHECK_EQUAL(6u, m.n_faces());
    CHECK( m.is_quad_mesh() );
    CHECK( !m.is_triangle_mesh() );
}

// A reader that scans coordinates through a mismatched format specifier keeps
// the connectivity intact while every position becomes garbage, so the counts
// above are not enough -- pin the actual geometry down.
TEST(ReadOff_GeometryValues)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );

    Point lo, hi;
    boundingBox(m, lo, hi);
    for (index_t i = 0; i != 3; ++i)
    {
        CHECK_CLOSE(-1.0, lo[i], 1e-12);
        CHECK_CLOSE( 1.0, hi[i], 1e-12);
    }

    // 12 edges of length 2
    CHECK_CLOSE(24.0, totalEdgeLength(m), 1e-12);

    // every vertex of this cube is a corner (+-1,+-1,+-1)
    for (auto v : m.vertices())
        for (index_t i = 0; i != 3; ++i)
            CHECK_CLOSE(1.0, math::abs(m.position(v)[i]), 1e-12);
}

// mesh.read() and gsReadFile() reach the mesh through completely different
// code paths (direct reader vs. gsFileData/gsXml); they must agree.
TEST(ReadPathsAgree)
{
    const std::string rel = "off/cube.off";

    Mesh direct;
    CHECK( readMesh(direct, rel) );

    Mesh viaXml;
    gsReadFile<real_t>(gsFileManager::findInDataDir(rel), viaXml);

    CHECK_EQUAL(direct.n_vertices(), viaXml.n_vertices());
    CHECK_EQUAL(direct.n_edges(),    viaXml.n_edges());
    CHECK_EQUAL(direct.n_faces(),    viaXml.n_faces());
    CHECK_CLOSE(totalEdgeLength(direct), totalEdgeLength(viaXml), 1e-12);

    for (auto v : direct.vertices())
        CHECK_CLOSE(0.0, (direct.position(v) - viaXml.position(v)).norm(), 1e-12);
}


//---------------------------------------------------------------------------
// Value semantics across the gsSurfMeshTopology / gsSurfMesh boundary
//---------------------------------------------------------------------------

// Copy/assign must duplicate topology *and* positions, and must re-bind the
// property handles so that the copy is independent of the original.
TEST(ValueSemantics_DeepCopy)
{
    Mesh a;
    CHECK( readMesh(a, "off/cube.off") );
    const Point p0 = a.position(Vertex(0));

    // copy constructor
    Mesh b(a);
    CHECK_EQUAL(a.n_vertices(), b.n_vertices());
    CHECK_EQUAL(a.n_edges(),    b.n_edges());
    CHECK_EQUAL(a.n_faces(),    b.n_faces());
    CHECK_CLOSE(0.0, (a.position(Vertex(0)) - b.position(Vertex(0))).norm(), 1e-12);

    // the copy owns its own point array: mutating it must not touch the source
    b.position(Vertex(0)) = Point(42, 43, 44);
    CHECK_CLOSE(0.0, (a.position(Vertex(0)) - p0).norm(), 1e-12);
    CHECK_CLOSE(42.0, b.position(Vertex(0))[0], 1e-12);

    // ... and the topology is independent too
    const size_t aFaces = a.n_faces();
    b.triangulate();
    CHECK_EQUAL(aFaces, a.n_faces());
    CHECK(b.n_faces() > aFaces);
}

TEST(ValueSemantics_AssignAndOperatorEq)
{
    Mesh a;
    CHECK( readMesh(a, "off/cube.off") );

    Mesh viaAssign;
    viaAssign.assign(a);
    CHECK_EQUAL(a.n_vertices(), viaAssign.n_vertices());
    CHECK_EQUAL(a.n_faces(),    viaAssign.n_faces());
    CHECK_CLOSE(totalEdgeLength(a), totalEdgeLength(viaAssign), 1e-12);

    Mesh viaOpEq;
    viaOpEq = a;
    CHECK_EQUAL(a.n_vertices(), viaOpEq.n_vertices());
    CHECK_EQUAL(a.n_faces(),    viaOpEq.n_faces());
    CHECK_CLOSE(totalEdgeLength(a), totalEdgeLength(viaOpEq), 1e-12);

    // both routes must reproduce the coordinates, not just the counts
    for (auto v : a.vertices())
    {
        CHECK_CLOSE(0.0, (a.position(v) - viaAssign.position(v)).norm(), 1e-12);
        CHECK_CLOSE(0.0, (a.position(v) - viaOpEq.position(v)).norm(), 1e-12);
    }
}

TEST(ValueSemantics_Move)
{
    Mesh a;
    CHECK( readMesh(a, "off/cube.off") );
    const real_t len = totalEdgeLength(a);
    const size_t nv = a.n_vertices(), nf = a.n_faces();

    Mesh moved;
    moved = give(a);
    CHECK_EQUAL(nv, moved.n_vertices());
    CHECK_EQUAL(nf, moved.n_faces());
    CHECK_CLOSE(len, totalEdgeLength(moved), 1e-12);
}


//---------------------------------------------------------------------------
// Topology invariants
//---------------------------------------------------------------------------

// V-E+F is 2-2g; these three files pin down genus 0, 1 and 2 respectively.
TEST(EulerCharacteristic)
{
    Mesh cube, torus, bitorus;
    CHECK( readMesh(cube,    "off/cube.off") );
    CHECK( readMesh(torus,   "off/octtorus.off") );
    CHECK( readMesh(bitorus, "off/bitorus.off") );

    CHECK_EQUAL( 2, euler(cube) );      // sphere
    CHECK_EQUAL( 0, euler(torus) );     // genus 1
    CHECK_EQUAL(-2, euler(bitorus) );   // genus 2
}

// quad_split() is Catmull-Clark refinement: it must preserve the topology type
// and the Euler characteristic, and grow the counts by the expected factors.
TEST(QuadSplit_Invariants)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    const size_t v0 = m.n_vertices(), e0 = m.n_edges(), f0 = m.n_faces();

    m.quad_split();

    // one new vertex per old edge and per old face
    CHECK_EQUAL(v0 + e0 + f0, m.n_vertices());
    CHECK_EQUAL(4 * f0,       m.n_faces());
    CHECK_EQUAL(2 * e0 + 4 * f0, m.n_edges());
    CHECK_EQUAL(2, euler(m));
    CHECK( m.is_quad_mesh() );

    // refinement must not move the surface outside its original bounding box
    Point lo, hi;
    boundingBox(m, lo, hi);
    for (index_t i = 0; i != 3; ++i)
    {
        CHECK(lo[i] >= -1.0 - 1e-12);
        CHECK(hi[i] <=  1.0 + 1e-12);
    }
}

TEST(Triangulate_Invariants)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    const index_t chi = euler(m);

    m.triangulate();

    CHECK( m.is_triangle_mesh() );
    CHECK_EQUAL(chi, euler(m));
    // a closed triangle mesh has 2E = 3F
    CHECK_EQUAL(3 * m.n_faces(), 2 * m.n_edges());
    for (auto f : m.faces())
        CHECK_EQUAL(3u, m.valence(f));
}

// Exercises the base class' allocation/deletion bookkeeping: after collecting
// the garbage the mesh must be dense again and still self-consistent.
TEST(DeleteAndGarbageCollection)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    const size_t f0 = m.n_faces();

    m.delete_face(*m.faces_begin());
    CHECK_EQUAL(f0 - 1, m.n_faces());

    m.garbage_collection();
    CHECK_EQUAL(f0 - 1, m.n_faces());

    // no deleted elements survive a garbage collection
    for (auto v : m.vertices()) CHECK( !m.is_deleted(v) );
    for (auto f : m.faces())    CHECK( !m.is_deleted(f) );
    for (auto e : m.edges())    CHECK( !m.is_deleted(e) );

    // removing one face of a closed cube opens a boundary
    size_t nb = 0;
    for (auto v : m.vertices()) if (m.is_boundary(v)) ++nb;
    CHECK(nb > 0);
}

// Building a mesh through the base-class face API, with geometry supplied by
// the derived class.
TEST(BuildFromScratch)
{
    Mesh m;
    Vertex a = m.add_vertex(Point(0,0,0));
    Vertex b = m.add_vertex(Point(1,0,0));
    Vertex c = m.add_vertex(Point(0,1,0));

    CHECK_EQUAL(3u, m.n_vertices());
    CHECK_CLOSE(1.0, m.position(b)[0], 1e-12);

    m.add_triangle(a, b, c);
    CHECK_EQUAL(1u, m.n_faces());
    CHECK_EQUAL(3u, m.n_edges());
    CHECK( m.is_triangle_mesh() );
    CHECK_EQUAL(3u, m.valence(*m.faces_begin()));

    // a lone triangle is all boundary
    CHECK( m.is_boundary(a) );
    CHECK( m.is_boundary(b) );
    CHECK( m.is_boundary(c) );

    // 1 + 1 + sqrt(2)
    CHECK_CLOSE(2.0 + math::sqrt(2.0), totalEdgeLength(m), 1e-12);
}


//---------------------------------------------------------------------------
// Functions recovered by hand rather than from version control
//
// These were added to the working tree but never committed, and were restored
// by transcription.  Unlike the rest of the class -- whose bodies are
// byte-identical to their previous revision and therefore carry no new risk --
// their only guarantee is the checks below, so keep them pinned.
//---------------------------------------------------------------------------

TEST(AddBatchVertices)
{
    Mesh m;
    m.add_batch_vertices(5);
    CHECK_EQUAL(5u, m.n_vertices());
    CHECK_EQUAL(0u, m.n_faces());

    // the batch is plain new vertices, so the point property must exist and be
    // addressable for each of them
    for (auto v : m.vertices())
        CHECK_CLOSE(0.0, m.position(v).norm(), 1e-12);
}

// The dual of a cube is an octahedron: vertices and faces swap, chi is kept.
TEST(DualMesh_CubeGivesOctahedron)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );

    Mesh d = m.dual_mesh();
    CHECK_EQUAL(m.n_faces(),    d.n_vertices());
    CHECK_EQUAL(m.n_vertices(), d.n_faces());
    CHECK_EQUAL(m.n_edges(),    d.n_edges());
    CHECK_EQUAL(2, euler(d));

    // the source must be left untouched by the non-mutating overload
    CHECK_EQUAL(8u, m.n_vertices());
    CHECK_EQUAL(6u, m.n_faces());

    // the in-place overload must agree with it
    Mesh inplace;
    CHECK( readMesh(inplace, "off/cube.off") );
    inplace.dual_mesh_inplace();
    CHECK_EQUAL(d.n_vertices(), inplace.n_vertices());
    CHECK_EQUAL(d.n_edges(),    inplace.n_edges());
    CHECK_EQUAL(d.n_faces(),    inplace.n_faces());

    // dual vertices are face barycentres, so they stay inside the cube
    Point lo, hi;
    boundingBox(d, lo, hi);
    for (index_t i = 0; i != 3; ++i)
    {
        CHECK(lo[i] >= -1.0 - 1e-12);
        CHECK(hi[i] <=  1.0 + 1e-12);
    }
}

TEST(AddMesh_Merge)
{
    Mesh a, b;
    CHECK( readMesh(a, "off/cube.off") );
    CHECK( readMesh(b, "off/cube.off") );
    const size_t va = a.n_vertices(), ea = a.n_edges(), fa = a.n_faces();

    gsVector<Vertex> idmap = a.add_mesh(b);

    CHECK_EQUAL((index_t)b.n_vertices(), idmap.size());
    CHECK_EQUAL(2 * va, a.n_vertices());
    CHECK_EQUAL(2 * ea, a.n_edges());
    CHECK_EQUAL(2 * fa, a.n_faces());
    // two disjoint closed surfaces
    CHECK_EQUAL(4, euler(a));

    // the returned map must send each submesh vertex to a vertex of the host
    // carrying the same position
    for (auto v : b.vertices())
    {
        const Vertex g = idmap[v.idx()];
        CHECK( g.is_valid() );
        CHECK_CLOSE(0.0, (b.position(v) - a.position(g)).norm(), 1e-12);
    }
}

// w=2 splits every edge in three, so each quad becomes 3x3.
TEST(QuadSplit_Width2)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    const size_t v0 = m.n_vertices(), e0 = m.n_edges(), f0 = m.n_faces();

    m.quad_split(2);

    CHECK_EQUAL(9 * f0, m.n_faces());                   // 3x3 per face
    CHECK_EQUAL(v0 + 2*e0 + 4*f0, m.n_vertices());      // corners + 2/edge + 4/face
    CHECK_EQUAL(2, euler(m));
    CHECK( m.is_quad_mesh() );

    Point lo, hi;
    boundingBox(m, lo, hi);
    for (index_t i = 0; i != 3; ++i)
    {
        CHECK(lo[i] >= -1.0 - 1e-12);
        CHECK(hi[i] <=  1.0 + 1e-12);
    }
}

// w=0 is a no-op and w=1 must coincide with the plain quad_split()
TEST(QuadSplit_WidthDegenerateCases)
{
    Mesh zero;
    CHECK( readMesh(zero, "off/cube.off") );
    zero.quad_split(0);
    CHECK_EQUAL(8u, zero.n_vertices());
    CHECK_EQUAL(6u, zero.n_faces());

    Mesh one, plain;
    CHECK( readMesh(one,   "off/cube.off") );
    CHECK( readMesh(plain, "off/cube.off") );
    one.quad_split(1);
    plain.quad_split();
    CHECK_EQUAL(plain.n_vertices(), one.n_vertices());
    CHECK_EQUAL(plain.n_edges(),    one.n_edges());
    CHECK_EQUAL(plain.n_faces(),    one.n_faces());
}


//---------------------------------------------------------------------------
// Normals (the triangle branch of compute_face_normal)
//---------------------------------------------------------------------------

// Uses a triangulated cube so that compute_face_normal() takes its triangle
// branch, where the cross product operands are easy to get wrong.
TEST(FaceNormals_Outward)
{
    Mesh m;
    CHECK( readMesh(m, "off/cube.off") );
    m.triangulate();
    CHECK( m.is_triangle_mesh() );

    for (auto f : m.faces())
    {
        const Point n = m.compute_face_normal(f);
        CHECK_CLOSE(1.0, n.norm(), 1e-10);

        // faces of an axis-aligned cube have axis-aligned normals
        index_t nonzero = 0;
        for (index_t i = 0; i != 3; ++i)
            if (math::abs(n[i]) > 1e-10) ++nonzero;
        CHECK_EQUAL(1, nonzero);

        // ... and they point away from the centre (the cube is origin-centred)
        CHECK( n.dot(m.face_barycenter(f)) > 0.0 );
    }
}

/// The normal of a counter-clockwise triangle must be (p1-p0) x (p2-p0),
/// normalized -- at ANY position, not only at the origin.
TEST(computeFaceNormal_matchesTheTextbookFormula)
{
    // Deliberately off-origin and not axis aligned: at the origin a sign error
    // can survive a symmetric fixture.
    const Point a(5.0, 3.0, 2.0);
    const Point b(6.0, 3.0, 2.0);
    const Point c(5.0, 4.0, 2.0);

    Mesh m = triangleMesh(a, b, c);
    const Normal n = m.compute_face_normal(*m.faces().begin());

    Normal expect = (b-a).cross(c-a);
    expect.normalize();

    CHECK_CLOSE(expect[0], n[0], 1e-12);
    CHECK_CLOSE(expect[1], n[1], 1e-12);
    CHECK_CLOSE(expect[2], n[2], 1e-12);

    // For this fixture that is +z; assert it explicitly so a future refactor
    // that flips the winding convention has to say so out loud.
    CHECK_CLOSE(0.0, n[0], 1e-12);
    CHECK_CLOSE(0.0, n[1], 1e-12);
    CHECK_CLOSE(1.0, n[2], 1e-12);
}

/// A normal is a direction: translating the mesh must not change it.
TEST(computeFaceNormal_isTranslationInvariant)
{
    const Point a(0.0, 0.0, 0.0);
    const Point b(1.0, 0.0, 0.0);
    const Point c(0.0, 1.0, 0.0);
    const Point t(5.0, 3.0, 2.0);

    Mesh m0 = triangleMesh(a,     b,     c    );
    Mesh m1 = triangleMesh(a + t, b + t, c + t);

    const Normal n0 = m0.compute_face_normal(*m0.faces().begin());
    const Normal n1 = m1.compute_face_normal(*m1.faces().begin());

    CHECK_CLOSE(n0[0], n1[0], 1e-12);
    CHECK_CLOSE(n0[1], n1[1], 1e-12);
    CHECK_CLOSE(n0[2], n1[2], 1e-12);
}

/// The triangle branch and the general-polygon branch must agree: a coplanar
/// quad and a triangle with the same winding report the same normal.
TEST(computeFaceNormal_triangleAndPolygonBranchesAgree)
{
    const Point a(5.0, 3.0, 2.0);
    const Point b(6.0, 3.0, 2.0);
    const Point c(6.0, 4.0, 2.0);
    const Point d(5.0, 4.0, 2.0);

    Mesh tri = triangleMesh(a, b, c);

    Mesh quad;
    std::vector<Vertex> qv;
    qv.push_back(quad.add_vertex(a));
    qv.push_back(quad.add_vertex(b));
    qv.push_back(quad.add_vertex(c));
    qv.push_back(quad.add_vertex(d));
    quad.add_face(qv);

    const Normal nt = tri .compute_face_normal(*tri .faces().begin());
    const Normal nq = quad.compute_face_normal(*quad.faces().begin());

    CHECK_CLOSE(nt[0], nq[0], 1e-12);
    CHECK_CLOSE(nt[1], nq[1], 1e-12);
    CHECK_CLOSE(nt[2], nq[2], 1e-12);
}

} // SUITE
