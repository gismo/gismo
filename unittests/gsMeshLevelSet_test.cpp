/** @file gsMeshLevelSet_test.cpp

    @brief Correctness of gsMeshSignedDist's analytic gradient at SHARP
           mesh features.

    PURPOSE. gsMeshSignedDist has two independent notions of "inside":

      - eval_into() takes the sign from the generalized WINDING NUMBER. This
        is the trusted oracle: it is a topological quantity, robust to mesh
        orientation (the class thresholds |omega| rather than reading its
        sign) and to sharp features.
      - deriv_into() takes the sign from the closest triangle's NORMAL, as
        sgn = dot(p - closest, n). This exists because it costs ONE BVH query
        instead of the twelve a finite-difference stencil charges, and the
        winding traversal is ~45% of a query.

    Those two must agree. When they do not, phi and grad(phi) disagree about
    which side of the surface a point is on, which corrupts anything that
    consumes the interface normal -- the Nitsche terms, and Algoim's interval
    bounds and height-direction selection.

    WHEN CAN THEY DISAGREE? Let theta be the interior dihedral angle at the
    closest edge. A point whose closest point lies on that edge has
    (p - closest) inside the normal cone spanned by the two incident face
    normals n1, n2. Writing it as a*n1 + b*n2 with a,b >= 0,

        dot(p - closest, n1) = a + b*cos(phi_n),   phi_n = angle(n1, n2)

    which can only go negative when cos(phi_n) < 0, i.e. phi_n > 90 deg, i.e.

        FACE NORMAL SIGN IS CORRECT  <=>  90 deg <= theta <= 270 deg.

    So the defect is invisible on smooth meshes AND on box-like geometry:
    an L-shaped solid has only 90 and 270 degree dihedrals and sits exactly
    on the boundary of correctness. It needs genuinely SHARP features, which
    is what filedata/obj/wedge.obj (theta ~ 25 deg, convex -- fails for
    EXTERIOR points) and filedata/obj/vgroove.obj (theta ~ 330 deg, concave
    -- fails for INTERIOR points) are for. Measured on the bundled meshes:
    spot.obj has 0 such edges and homer.obj has 10, so neither is an adequate
    fixture on its own.

    WHAT IS EVIDENCE HERE. Not a reference value: the test compares the
    analytic gradient against a CENTRAL DIFFERENCE of eval_into(). Because
    eval_into's sign comes from the winding number, that finite difference is
    an independent oracle, and a sign error shows up unmissably as a ~180 deg
    angle between the two gradients. No golden numbers to refresh.

    Sample points are placed OFF the sharp edge along the two face normals
    and their bisector, at offsets large enough to be well outside the
    finite-difference step and _gradEps(), so the comparison is
    well-conditioned and does not probe the non-differentiable point itself.

    A NOTE ON WELDING, because it bounds what any of this can promise. The
    pseudonormals are accumulated per MESH vertex and per MESH edge. If a file
    carries geometrically coincident but topologically distinct vertices --
    routine in .obj and STL exported per-face -- the halfedge mesh sees two
    BOUNDARY edges rather than one shared edge, each accumulator receives only
    one incident face, and the "edge pseudonormal" collapses back to that face's
    normal. The fix then silently does nothing at exactly the features it exists
    for. gsSurfMeshBVH warns when a component has boundary edges, which is the
    detectable symptom; sliver.obj deliberately contains such a duplicated pair
    so the behaviour is pinned rather than assumed. Welding is the caller's
    responsibility.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <gsDomain/gsMeshLevelSet.h>

namespace {

/// Reads a mesh from the data dir, failing the test rather than the process
/// if it is missing.
bool loadMesh(const std::string & rel, gsSurfMesh & mesh)
{
    return gsReadSurfMesh(rel, mesh);
}

/// Central difference of phi, i.e. the gradient implied by eval_into()'s
/// winding-number sign. Independent of deriv_into() by construction.
gsVector<real_t,3> fdGradient(const gsFunction<real_t> & phi,
                              const gsVector<real_t,3> & p,
                              real_t h)
{
    gsVector<real_t,3> g;
    for (short_t c = 0; c != 3; ++c)
    {
        gsMatrix<real_t> a(3,1), b(3,1), va, vb;
        a = p; b = p;
        a(c,0) -= h;  b(c,0) += h;
        phi.eval_into(a, va);
        phi.eval_into(b, vb);
        g[c] = (vb(0,0) - va(0,0)) / (2*h);
    }
    return g;
}

/// Worst angle (degrees) between the analytic and finite-difference gradients
/// over \a pts, the number of outright SIGN disagreements (angle > 90), and
/// the number of samples actually compared.
///
/// Samples where the finite difference is not EIKONAL are skipped. phi is a
/// signed distance, so |grad phi| == 1 wherever it is differentiable; the
/// central difference reports a norm well below 1 exactly where it straddles a
/// discontinuity of the gradient -- the medial axis, or the surface itself.
/// Comparing against a finite difference there is meaningless: BOTH answers are
/// legitimate one-sided limits and they genuinely differ. Filtering on |gN| is
/// self-contained (no geometric knowledge of the mesh required) and keeps the
/// comparison to the region where a unique correct answer exists.
void compareGradients(const gsFunction<real_t>        & phi,
                      const std::vector<gsVector<real_t,3> > & pts,
                      real_t h, real_t & maxAngleDeg, index_t & nFlipped,
                      index_t & nCompared)
{
    maxAngleDeg = 0;
    nFlipped    = 0;
    nCompared   = 0;
    for (size_t i = 0; i != pts.size(); ++i)
    {
        gsMatrix<real_t> q(3,1), ga;
        q = pts[i];
        phi.deriv_into(q, ga);

        gsVector<real_t,3> gA = ga.col(0);
        gsVector<real_t,3> gN = fdGradient(phi, pts[i], h);

        const real_t la = gA.norm(), ln = gN.norm();
        if (la < 1e-12) continue;
        if (ln < 0.9 || ln > 1.1) continue;        // not eikonal here -- skip

        ++nCompared;
        const real_t c   = math::max(real_t(-1),
                           math::min(real_t(1), gA.dot(gN)/(la*ln)));
        const real_t ang = math::acos(c) * 180 / EIGEN_PI;
        maxAngleDeg = math::max(maxAngleDeg, ang);
        if (ang > 90) ++nFlipped;
    }
}

/// Sample points around a straight sharp edge running along direction \a dir
/// through \a origin, offset into the half-space spanned by the two incident
/// face normals (or their negatives, for a concave edge).
std::vector<gsVector<real_t,3> >
edgeSamples(const gsVector<real_t,3> & origin,
            const gsVector<real_t,3> & dir,
            const gsVector<real_t,3> & n1,
            const gsVector<real_t,3> & n2,
            real_t sign, index_t nAlong, index_t nAround, real_t offset)
{
    std::vector<gsVector<real_t,3> > pts;
    for (index_t i = 1; i <= nAlong; ++i)
    {
        // Stay away from the edge's own endpoints, where OTHER features win.
        const real_t t = real_t(i) / real_t(nAlong + 1);
        const gsVector<real_t,3> base = origin + dir * t;
        for (index_t j = 0; j <= nAround; ++j)
        {
            const real_t s = real_t(j) / real_t(nAround);
            gsVector<real_t,3> d = n1 * (1-s) + n2 * s;
            if (d.norm() < 1e-12) continue;
            d.normalize();
            pts.push_back(base + sign * offset * d);
        }
    }
    return pts;
}

} // anonymous namespace

SUITE(gsMeshLevelSet_test)
{

/// Control: a smooth mesh has no dihedral outside [90,270], so the analytic
/// gradient must already agree with the finite difference everywhere. This
/// guards against a "fix" that perturbs the cases that were always correct.
TEST(smoothMesh_gradientAgreesWithWindingFD)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/spot.obj", mesh))
    { CHECK(false); return; }                 // data file missing

    gsNormalizeToUnitBox(mesh, 0.8);
    gsMeshSignedDist<real_t> phi(mesh, gsUnitBox3());

    // A shell of points around the surface, at a radius where the signed
    // distance is smooth.
    std::vector<gsVector<real_t,3> > pts;
    for (index_t i = 0; i < 12; ++i)
        for (index_t j = 0; j < 12; ++j)
        {
            const real_t a = EIGEN_PI * real_t(i) / 11;
            const real_t b = 2*EIGEN_PI * real_t(j) / 12;
            gsVector<real_t,3> p(3);
            p << 0.5 + 0.45*math::sin(a)*math::cos(b),
                 0.5 + 0.45*math::sin(a)*math::sin(b),
                 0.5 + 0.45*math::cos(a);
            pts.push_back(p);
        }

    real_t  maxAng = 0;
    index_t nFlip = 0, nCmp = 0;
    compareGradients(phi, pts, 1e-5, maxAng, nFlip, nCmp);

    CHECK(nCmp > 50);          // the filter must not have emptied the sample
    CHECK_EQUAL(0, nFlip);
    CHECK(maxAng < 5.0);
}

/// Degeneracy guards. filedata/obj/sliver.obj is an ordinary unit cube except
/// for two deliberate pathologies: a corner shaved by 1e-12 (a triangle whose
/// EDGES are ~1.4e-12 long) and a split top edge (a triangle of area ~5e-15
/// with nearly collinear vertices). Together they hit both guards in
/// _computeNormals(): the zero-length-edge test that keeps acos() away from
/// 0/0, and the zero-area test that keeps a cancellation-dominated cross
/// product out of the accumulators.
///
/// This matters more than a normal being slightly off. A NaN pseudonormal
/// makes the comparison `d.dot(n) < 0` FALSE (every comparison against NaN is),
/// so the sign silently pins to +1 -- "outside" -- at every query landing on
/// that feature, with no warning and no obviously wrong value to notice.
/// Hence the explicit finiteness sweep over every feature of every triangle,
/// rather than only checking the queries a sample happens to reach.
TEST(sliverTriangles_produceFiniteNormalsAndCorrectSigns)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/sliver.obj", mesh))
    { CHECK(false); return; }

    gsSurfMeshBVH bvh;
    bvh.build(mesh);
    CHECK(bvh.numTriangles() > 0);

    // (a) No NaN or Inf anywhere in the precomputed normals.
    index_t nBadNormal = 0;
    for (index_t t = 0; t != bvh.numTriangles(); ++t)
        for (short_t f = 0; f != 7; ++f)
        {
            const gsSurfMesh::Point n = bvh.pseudoNormal(t, f);
            if (!(n.array() == n.array()).all()) { ++nBadNormal; continue; }
            if (!n.allFinite()) ++nBadNormal;
        }
    CHECK_EQUAL(0, nBadNormal);

    // (b) Signs are still right. Away from the perturbed corner and edge (both
    // within 1e-12 of the nominal cube) the exact answer is the unit cube's.
    gsMeshSignedDist<real_t> phi(mesh, gsSurfMeshBoundingBox(mesh));

    index_t nWrongSign = 0, nNonFinite = 0, nNonUnit = 0;
    for (index_t i = 0; i <= 6; ++i)
        for (index_t j = 0; j <= 6; ++j)
            for (index_t k = 0; k <= 6; ++k)
            {
                gsMatrix<real_t> q(3,1), v, g;
                q << -0.25 + 0.25*i, -0.25 + 0.25*j, -0.25 + 0.25*k;

                const bool inside = (q(0,0) > 0 && q(0,0) < 1 &&
                                     q(1,0) > 0 && q(1,0) < 1 &&
                                     q(2,0) > 0 && q(2,0) < 1);
                phi.eval_into (q, v);
                phi.deriv_into(q, g);

                if (!math::isfinite(v(0,0)) || !g.col(0).allFinite())
                { ++nNonFinite; continue; }

                if ((v(0,0) < 0) != inside) ++nWrongSign;

                // phi is a signed distance, so |grad phi| == 1 off the medial
                // axis. A pinned/garbage normal shows up here even when the
                // sign happens to come out right.
                const real_t gl = g.col(0).norm();
                if (gl < 0.99 || gl > 1.01) ++nNonUnit;
            }

    CHECK_EQUAL(0, nNonFinite);
    CHECK_EQUAL(0, nWrongSign);
    CHECK_EQUAL(0, nNonUnit);
}

/// Coincident-but-distinct vertices must not split the pseudonormal.
///
/// sliver.obj stores the corner (1,1,1) under TWO vertex indices. Without
/// welding, each index accumulates only the faces incident to it, so the
/// pseudonormal at that corner is partial -- in the worst case a single face
/// normal, which has two zero components. With welding it sees all three cube
/// faces meeting there, so every component is strictly positive (the corner
/// bisector points into the (+,+,+) octant).
///
/// That component-sign test is the assertion rather than an exact vector,
/// because the corner also carries a 1e-12 shave: the exact bisector depends on
/// the shaved geometry, but "all three components positive" separates a welded
/// three-face average from any single face normal cleanly.
TEST(duplicatedVertices_areWeldedForNormalAccumulation)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/sliver.obj", mesh))
    { CHECK(false); return; }

    gsMeshSignedDist<real_t> phi(mesh, gsSurfMeshBoundingBox(mesh));

    gsMatrix<real_t> q(3,1), g;
    q << 1, 1, 1;                       // exactly the duplicated corner
    phi.deriv_into(q, g);

    CHECK(g.col(0).allFinite());
    CHECK_CLOSE(1.0, g.col(0).norm(), 1e-9);
    CHECK(g(0,0) > 0.1);
    CHECK(g(1,0) > 0.1);
    CHECK(g(2,0) > 0.1);
}

/// Locks that the gradient ON a sharp edge is the angle-weighted BISECTOR and
/// not one of the two face normals.
///
/// This is the one assertion that distinguishes the pseudonormal machinery
/// from the old face-normal code: everywhere off the surface the gradient is
/// (p - closest)/dist and the normal only supplies a sign, so a wrong normal
/// hides. Exactly ON the surface the dist <= _gradEps() branch returns the
/// normal itself, which is observable.
///
/// wedge.obj is used rather than sliver.obj because it has no duplicated
/// vertices: the two coincident-but-distinct corner indices in sliver.obj
/// split the per-vertex accumulator between them, so the pseudonormal there is
/// partial by construction (see the note on welding in this file's header).
TEST(sharpEdge_gradientOnSurfaceIsTheBisector)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/wedge.obj", mesh))
    { CHECK(false); return; }

    gsMeshSignedDist<real_t> phi(mesh, gsSurfMeshBoundingBox(mesh));

    // A point exactly ON the apex edge (x=0, z=0), midway along it. The two
    // slanted faces have outward normals at +-(90+12.5) deg from +x in the
    // xz-plane, so their normalized sum -- the edge pseudonormal -- is exactly
    // -x. A face normal would instead be (-sin12.5, 0, -+cos12.5), i.e. ~77.5
    // degrees away.
    gsMatrix<real_t> q(3,1), g;
    q << 0, 0.5, 0;
    phi.deriv_into(q, g);

    CHECK(g.col(0).allFinite());
    CHECK_CLOSE(-1.0, g(0,0), 1e-9);
    CHECK_CLOSE( 0.0, g(1,0), 1e-9);
    CHECK_CLOSE( 0.0, g(2,0), 1e-9);
}

/// Locks WHERE the winding-number fallback fires. This is the cost model: the
/// fallback costs a winding traversal, so a mesh reporting zero sharp features
/// runs the fast path exclusively and pays nothing for the fix. If spot.obj
/// ever reports non-zero here, the FCM drivers just got slower and the reason
/// is a mesh change, not a code change.
TEST(sharpFeatureCount_isZeroOnSmoothMeshAndPositiveOnSharpOnes)
{
    gsSurfMesh smooth, wedge, groove;
    if (!loadMesh("obj/spot.obj",    smooth) ||
        !loadMesh("obj/wedge.obj",   wedge)  ||
        !loadMesh("obj/vgroove.obj", groove))
    { CHECK(false); return; }

    gsSurfMeshBVH bSmooth, bWedge, bGroove;
    bSmooth.build(smooth);
    bWedge .build(wedge);
    bGroove.build(groove);

    // spot.obj has ZERO sharp EDGES (every dihedral is inside [90,270]) but
    // two sharp VERTICES: 1845 and 2921, both of degree 6, where a pair of
    // incident face normals is 91.2 degrees apart -- marginally over. Vertex
    // sharpness is the stricter condition, because the normal cone at a vertex
    // is spanned by ALL incident faces, not just an adjacent pair; a faceted
    // cone can step around in safe increments and still have opposite facets
    // pointing more than 90 degrees apart.
    //
    // 2 vertices x degree 6 = 12 (triangle, feature) slots, out of
    // 5856 x 6 = 35136. So the winding fallback fires on ~0.03% of features
    // here: the cost of the fix on this mesh is nil.
    CHECK_EQUAL(12, bSmooth.numSharpFeatures());

    // The two fixtures exist precisely to be outside that window.
    CHECK(bWedge .numSharpFeatures() > 0);
    CHECK(bGroove.numSharpFeatures() > 0);
}

/// Sharp CONVEX edge (wedge.obj, dihedral ~25 deg). Exterior points whose
/// closest point is the apex edge are the failure case: the face normal of
/// whichever incident triangle the BVH happened to return can be more than
/// 90 deg away from (p - closest), inverting the sign.
TEST(sharpConvexEdge_gradientSignIsCorrect)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/wedge.obj", mesh))
    { CHECK(false); return; }

    gsMeshSignedDist<real_t> phi(mesh, gsSurfMeshBoundingBox(mesh));

    // Apex edge runs along +y at (x,z) = (0,0); see the mesh's own header.
    // The two slanted faces make +-12.5 deg with the xy-plane, so their
    // outward normals are rotated +-(90+12.5) deg from +x in the xz-plane.
    const real_t half = 12.5 * EIGEN_PI / 180;
    gsVector<real_t,3> origin(3), dir(3), n1(3), n2(3);
    origin << 0, 0, 0;
    dir    << 0, 1, 0;
    n1     << -math::sin(half), 0, -math::cos(half);
    n2     << -math::sin(half), 0,  math::cos(half);

    // OUTSIDE the wedge is -x of the apex; offset along the normal cone.
    std::vector<gsVector<real_t,3> > pts =
        edgeSamples(origin, dir, n1, n2, +1, 8, 16, 0.02);

    real_t  maxAng = 0;
    index_t nFlip = 0, nCmp = 0;
    compareGradients(phi, pts, 1e-6, maxAng, nFlip, nCmp);

    CHECK(nCmp > 50);          // the filter must not have emptied the sample
    CHECK_EQUAL(0, nFlip);
    CHECK(maxAng < 5.0);
}

/// Sharp CONCAVE edge (vgroove.obj, dihedral ~330 deg). The mirror case:
/// here it is INTERIOR points whose closest point is the groove bottom, and
/// the face normal can report them as outside. A rectangular slot (270 deg)
/// would not exercise this at all.
TEST(sharpConcaveEdge_gradientSignIsCorrect)
{
    gsSurfMesh mesh;
    if (!loadMesh("obj/vgroove.obj", mesh))
    { CHECK(false); return; }

    gsMeshSignedDist<real_t> phi(mesh, gsSurfMeshBoundingBox(mesh));

    // Groove bottom edge runs along +y at (x,z) = (0.5, 0.5). The two groove
    // walls make +-15 deg with the vertical, so their outward normals (which
    // point INTO the groove, i.e. away from the solid) are rotated
    // correspondingly about y.
    const real_t half = 15.0 * EIGEN_PI / 180;
    gsVector<real_t,3> origin(3), dir(3), n1(3), n2(3);
    origin << 0.5, 0, 0.5;
    dir    << 0,   1, 0;
    n1     <<  math::cos(half), 0, math::sin(half);
    n2     << -math::cos(half), 0, math::sin(half);

    // INSIDE the solid is below the groove bottom: offset along -(normals).
    std::vector<gsVector<real_t,3> > pts =
        edgeSamples(origin, dir, n1, n2, -1, 8, 16, 0.02);

    real_t  maxAng = 0;
    index_t nFlip = 0, nCmp = 0;
    compareGradients(phi, pts, 1e-6, maxAng, nFlip, nCmp);

    CHECK(nCmp > 50);          // the filter must not have emptied the sample
    CHECK_EQUAL(0, nFlip);
    CHECK(maxAng < 5.0);
}

}
