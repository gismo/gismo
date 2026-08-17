/** @file gsSurfMesh_test.cpp

    @brief Regression lock for gsSurfMesh::compute_face_normal().

    The triangle branch of compute_face_normal() used to cross the RELATIVE
    (p2-p1) with the ABSOLUTE p1, instead of with the relative (p0-p1) that the
    general-polygon branch of the same function accumulates. Two consequences,
    both silent:

      - the normal was not translation invariant: moving a mesh changed the
        normals of its triangles;
      - it was sign-inverted even at the origin, so a consistently
        counter-clockwise triangle reported an INWARD normal.

    A face normal is a quantity nothing else can sanity-check for you -- it has
    no magnitude to look wrong and no obviously invalid value -- so this pins
    both properties directly:

      1. against the textbook (p1-p0) x (p2-p0) for a triangle at a general
         position (not at the origin, where a sign error can hide behind a
         symmetric fixture);
      2. invariance under a rigid translation;
      3. agreement between the triangle branch and the polygon branch, which is
         what made the discrepancy findable in the first place -- a quad and the
         triangle cutting off one of its corners must report the same normal
         when they are coplanar.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <gsMesh2/gsSurfMesh.h>

#include <vector>

namespace {

/// One triangle with the given corner positions, in the given order.
gsSurfMesh triangleMesh(const gsSurfMesh::Point & a,
                        const gsSurfMesh::Point & b,
                        const gsSurfMesh::Point & c)
{
    gsSurfMesh m;
    std::vector<gsSurfMesh::Vertex> v;
    v.push_back(m.add_vertex(a));
    v.push_back(m.add_vertex(b));
    v.push_back(m.add_vertex(c));
    m.add_face(v);
    return m;
}

} // anonymous namespace

SUITE(gsSurfMesh_test)
{

/// The normal of a counter-clockwise triangle must be (p1-p0) x (p2-p0),
/// normalized -- at ANY position, not only at the origin.
TEST(computeFaceNormal_matchesTheTextbookFormula)
{
    // Deliberately off-origin and not axis aligned: at the origin a sign error
    // can survive a symmetric fixture.
    const gsSurfMesh::Point a(5.0, 3.0, 2.0);
    const gsSurfMesh::Point b(6.0, 3.0, 2.0);
    const gsSurfMesh::Point c(5.0, 4.0, 2.0);

    gsSurfMesh m = triangleMesh(a, b, c);
    const gsSurfMesh::Normal n = m.compute_face_normal(*m.faces().begin());

    gsSurfMesh::Normal expect = (b-a).cross(c-a);
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

/// A normal is a direction: translating the mesh must not change it. This is
/// the property the old code violated.
TEST(computeFaceNormal_isTranslationInvariant)
{
    const gsSurfMesh::Point a(0.0, 0.0, 0.0);
    const gsSurfMesh::Point b(1.0, 0.0, 0.0);
    const gsSurfMesh::Point c(0.0, 1.0, 0.0);
    const gsSurfMesh::Point t(5.0, 3.0, 2.0);

    gsSurfMesh m0 = triangleMesh(a,     b,     c    );
    gsSurfMesh m1 = triangleMesh(a + t, b + t, c + t);

    const gsSurfMesh::Normal n0 = m0.compute_face_normal(*m0.faces().begin());
    const gsSurfMesh::Normal n1 = m1.compute_face_normal(*m1.faces().begin());

    CHECK_CLOSE(n0[0], n1[0], 1e-12);
    CHECK_CLOSE(n0[1], n1[1], 1e-12);
    CHECK_CLOSE(n0[2], n1[2], 1e-12);
}

/// The triangle branch and the general-polygon branch must agree: a coplanar
/// quad and a triangle with the same winding report the same normal. The two
/// branches disagreeing is what exposed the defect.
TEST(computeFaceNormal_triangleAndPolygonBranchesAgree)
{
    const gsSurfMesh::Point a(5.0, 3.0, 2.0);
    const gsSurfMesh::Point b(6.0, 3.0, 2.0);
    const gsSurfMesh::Point c(6.0, 4.0, 2.0);
    const gsSurfMesh::Point d(5.0, 4.0, 2.0);

    gsSurfMesh tri = triangleMesh(a, b, c);

    gsSurfMesh quad;
    std::vector<gsSurfMesh::Vertex> qv;
    qv.push_back(quad.add_vertex(a));
    qv.push_back(quad.add_vertex(b));
    qv.push_back(quad.add_vertex(c));
    qv.push_back(quad.add_vertex(d));
    quad.add_face(qv);

    const gsSurfMesh::Normal nt = tri .compute_face_normal(*tri .faces().begin());
    const gsSurfMesh::Normal nq = quad.compute_face_normal(*quad.faces().begin());

    CHECK_CLOSE(nt[0], nq[0], 1e-12);
    CHECK_CLOSE(nt[1], nq[1], 1e-12);
    CHECK_CLOSE(nt[2], nq[2], 1e-12);
}

}
