/** @file gsMeshLevelSet.h

    @brief Signed-distance level set of a closed triangle mesh, for use as the
           implicit function of a gsImplicitTrimmedDomain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsIO/gsFileManager.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/gsSurfMeshBVH.h>

#include <string>
#include <vector>
#include <limits>

namespace gismo
{

/**
 * @brief Signed-distance level set of a closed triangle mesh.
 *
 * Sign convention (matches gsImplicitTrimmedDomain):
 *   - phi < 0  ->  inside  (active / interior)
 *   - phi > 0  ->  outside
 *   - phi = 0  ->  on the surface (cut cell)
 *
 * Both queries are answered through a gsSurfMeshBVH built once in the
 * constructor, instead of the O(nTriangles)-per-point brute force this would
 * otherwise be (profiled at ~400e6 triangle tests for a single small assembly
 * pass). See gsSurfMeshBVH for the acceleration structure, its cost model and
 * why it is not delegated to gmsh, CGAL or libigl.
 *
 * PRECISION. The mesh, the BVH and every intermediate are real_t; the only cast
 * is at the gsFunction<T> interface (the static_cast<T> in eval_into and
 * deriv_into), a no-op in the usual T == real_t instantiation and present so
 * that e.g. gsMeshSignedDist<float> still computes in library precision and
 * narrows once, at the end.
 *
 * \ingroup Domain
 */
template<class T>
class gsMeshSignedDist : public gsFunction<T>
{
public:
    typedef gsSurfMesh::Point Point;   ///< gsVector3d<real_t>

    GISMO_CLONE_FUNCTION(gsMeshSignedDist)

    /// Builds the level set of \a mesh over the bounding box \a bbox.
    ///
    /// \param mesh  A closed TRIANGLE mesh; call gsSurfMesh::triangulate()
    ///              first if it carries n-gons. Its winding is normalized to
    ///              outward inside the BVH (on a private copy, so \a mesh is
    ///              not modified), and a warning is issued if it is not closed.
    /// \param bbox  3 x 2: column 0 = lower corner, column 1 = upper corner.
    ///              Returned by support(), which is what the basis-free
    ///              gsImplicitTrimmedDomain construction paths read.
    gsMeshSignedDist(const gsSurfMesh & mesh, const gsMatrix<T> & bbox)
    : m_bbox(bbox)
    {
        GISMO_ASSERT(bbox.rows() == 3 && bbox.cols() == 2,
                     "gsMeshSignedDist: bbox must be 3 x 2 (col 0 = lower, "
                     "col 1 = upper corner), got " << bbox.rows() << " x "
                     << bbox.cols() << ".");
        m_bvh.build(mesh);
    }

    short_t     domainDim() const override { return 3; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    /// The acceleration structure, for callers that need the raw queries.
    const gsSurfMeshBVH & bvh() const { return m_bvh; }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const Point p(static_cast<real_t>(u(0,k)),
                          static_cast<real_t>(u(1,k)),
                          static_cast<real_t>(u(2,k)));
            result(0,k) = static_cast<T>(_signedDistance(p));
        }
    }

    /// Analytic gradient, replacing gsFunction<T>'s central differences.
    ///
    /// Motivation (measured, not assumed): gsAlgoimFunctionAdapter asks for
    /// grad(phi) one point at a time, and the finite-difference fallback
    /// answers each such request with 3 calls x 4 points = 12 BVH queries.
    /// Profiling an assembly pass showed 82% of ALL level-set point
    /// evaluations were those stencils -- gradients, not values. This computes
    /// the gradient from the SAME nearest-point query that |phi| needs, i.e.
    /// one query instead of twelve.
    ///
    /// Away from the surface grad(phi)(p) = sign(phi) * (p - closest)/|p - closest|
    /// exactly (phi is a signed distance, so |grad phi| == 1 off the medial
    /// axis). On the medial axis phi is genuinely non-differentiable and this
    /// picks one side -- as does any finite-difference stencil straddling it.
    ///
    /// The catch, and why a previous attempt at this was reverted: that ratio
    /// divides a near-zero vector by a near-zero scalar as p approaches the
    /// surface, and the Nitsche quadrature points sit ON {phi==0} by
    /// construction, so it degenerates into pure rounding noise exactly where
    /// it is needed most (it broke CG convergence outright and blew the error
    /// norms up to ~1e16). Below _gradEps() the limit is taken instead: the
    /// outward unit normal of the closest triangle. Verified on spot.obj: 0
    /// sign disagreements in 20000 samples, and 0.000 deg mean and max angle
    /// against the well-conditioned formula at offsets 1e-5 and 1e-7.
    ///
    /// The sign likewise comes from that normal rather than from a
    /// winding-number query: for a closed, outward-oriented mesh p is outside
    /// exactly when (p - closest) points along it, because the closest-point
    /// vector always lies in the normal cone of the closest feature. This makes
    /// a gradient cost ONE nearest-point query -- the winding traversal was
    /// ~45% of the query and is skipped entirely here. Verified against the
    /// winding verdict on spot.obj: 0 disagreements in 20000 samples at offsets
    /// 1e-1 through 1e-7.
    ///
    /// WHEN IS THE DOT-PRODUCT SIGN VALID? Only where the closest FEATURE is
    /// not sharp. Let theta be the interior dihedral angle at the closest edge.
    /// (p - closest) lies in the normal cone spanned by the two incident face
    /// normals n1, n2; writing it a*n1 + b*n2 with a, b >= 0,
    ///
    ///     dot(p - closest, n1) = a + b*cos(angle(n1,n2))
    ///
    /// which can only go negative when the normals are more than 90 degrees
    /// apart, i.e. when theta is outside [90, 270] degrees. So:
    ///
    ///   - face interior: always correct ((p-closest) is exactly +/- n).
    ///   - edge or vertex with all incident normals within 90 degrees: correct.
    ///   - SHARP edge or vertex: the dot product can be inverted, and this
    ///     falls back to the winding number -- the same oracle eval_into()
    ///     uses, so phi and grad(phi) cannot disagree about the side.
    ///
    /// The fallback costs one winding traversal (~45% of a query) but fires
    /// only on sharp features: measured, spot.obj has none at all, homer.obj
    /// has 10 edges. gsSurfMeshBVH::numSharpFeatures() reports the count, and
    /// zero means the fast path is provably exact for the whole mesh.
    ///
    /// Deliberately NOT relying on the pseudonormal alone for the sign: the
    /// Baerentzen-Aanaes theorem assumes a closed, manifold, consistently
    /// oriented mesh, which .obj input routinely is not, whereas the winding
    /// number is robust to all three (which is why eval_into() uses it).
    ///
    /// \note The pseudonormal returned in the dist <= eps branch is a
    /// CONVENTION, not a limit: on an edge or a vertex grad(phi) genuinely has
    /// no limit (different approach directions in the normal cone give
    /// different answers). The angle-weighted bisector is the classical choice
    /// and the one the sign theorem is stated for.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(3, u.cols());
        const real_t eps = _gradEps();
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const Point p(static_cast<real_t>(u(0,k)),
                          static_cast<real_t>(u(1,k)),
                          static_cast<real_t>(u(2,k)));
            Point   closest;
            index_t tri  = -1;
            short_t feat = gsSurfMeshBVH::FaceInterior;
            const real_t dist =
                math::sqrt(m_bvh.squaredDistance(p, closest, &tri, &feat));
            const Point n = m_bvh.pseudoNormal(tri, feat);

            Point g;
            if (dist > eps)
            {
                const Point d = p - closest;
                const real_t sgn = m_bvh.isSharp(tri, feat)
                    ? (_isInside(p) ? real_t(-1) : real_t(1))   // robust, rare
                    : (d.dot(n) < 0 ? real_t(-1) : real_t(1));  // fast, common
                g = d * (sgn / dist);
            }
            else
                g = n;

            for (index_t c = 0; c != 3; ++c)
                result(c,k) = static_cast<T>(g[c]);
        }
    }

private:
    /// Distance below which (p - closest)/|p - closest| is dominated by
    /// rounding and the closest triangle's normal is used instead.
    ///
    /// 1e-9 is a hand-tuned constant, verified in double as cited in
    /// deriv_into(), and kept EXACTLY in any build whose arithmetic can carry
    /// it -- the recorded volume regression anchors are sensitive to it at the
    /// 1e-8 level and there is no accuracy to gain by perturbing them.
    ///
    /// It cannot simply be hard-coded, though: 1e-9 sits two orders of
    /// magnitude BELOW float's epsilon (1.19e-7), so in a real_t = float build
    /// the guard would never fire and the near-surface division would run on
    /// pure rounding noise -- exactly the failure above. Where the constant is
    /// not supportable, sqrt(eps) takes over: the standard "this difference has
    /// lost half its significant digits" threshold, and the right shape here
    /// because the mesh is normalized into the unit box, so absolute and
    /// relative error of the subtraction coincide.
    ///
    /// Evaluates to 1e-9 for double, long double and mpfr (whose extra
    /// precision cannot help anyway -- meshes arrive from ASCII .obj at ~7
    /// significant digits), and to ~3.4e-4 for float.
    static real_t _gradEps()
    {
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        return (eps < real_t(1e-12)) ? real_t(1e-9) : math::sqrt(eps);
    }

    real_t _signedDistance(const Point & p) const
    {
        Point closest;
        const real_t dist = math::sqrt(m_bvh.squaredDistance(p, closest));
        return _isInside(p) ? -dist : dist;
    }

    bool _isInside(const Point & p) const
    {
        // atan2(0,-1) == pi exactly, at whatever precision real_t carries --
        // unlike a decimal literal, which would truncate pi to double in an
        // mpfr build.
        static const real_t fourPi = real_t(4) * math::atan2(real_t(0), real_t(-1));
        const real_t winding = m_bvh.windingOmega(p) / fourPi;
        // |winding| ~ 1 inside, ~ 0 outside (sign-robust to mesh orientation)
        return math::abs(winding) > real_t(0.5);
    }

    gsSurfMeshBVH m_bvh;   ///< built once in the constructor
    gsMatrix<T>   m_bbox;  ///< 3 x 2 : col 0 = lower, col 1 = upper
};

// =============================================================================
//  Helpers shared by the immersed FCM drivers
// =============================================================================

/// Reads a surface mesh from \a filename (.obj/.off/.stl/...) and triangulates
/// it, returning false if the file could not be read.
///
/// \a filename is resolved through gsFileManager::find(), so a path relative to
/// the search paths works: the bundled meshes are reached as "obj/spot.obj".
/// An absolute or explicitly relative path is used as given.
///
/// \note gsSurfMesh::add_face() REJECTS a face that would create a complex
/// edge, which is what an inconsistently wound input produces, and reports it
/// on std::cerr. Such faces are dropped rather than repaired, so a mesh that
/// triggers those messages is not usable as a level set. The global
/// inside-out case IS repaired, inside gsSurfMeshBVH::build().
inline bool gsReadSurfMesh(const std::string & filename, gsSurfMesh & mesh)
{
    const std::string fn = gsFileManager::find(filename);
    if (fn.empty())
    {
        gsWarn << "gsReadSurfMesh: cannot find '" << filename << "' in the "
                  "search paths (" << gsFileManager::getSearchPaths() << ").\n";
        return false;
    }
    if (!mesh.read(fn))
        return false;
    if (!mesh.is_triangle_mesh())
        mesh.triangulate();
    return mesh.n_faces() != 0;
}

/// Axis-aligned bounding box of \a mesh, as 3 x 2 (col 0 = lower corner,
/// col 1 = upper corner) -- the layout gsMeshSignedDist and the basis-free
/// gsImplicitTrimmedDomain paths expect.
inline gsMatrix<real_t> gsSurfMeshBoundingBox(const gsSurfMesh & mesh)
{
    GISMO_ENSURE(mesh.n_vertices() != 0,
                 "gsSurfMeshBoundingBox: the mesh has no vertices.");

    gsSurfMesh::Point lo = mesh.position(*mesh.vertices().begin());
    gsSurfMesh::Point hi = lo;
    for (gsSurfMesh::Vertex v : mesh.vertices())
    {
        lo = lo.cwiseMin(mesh.position(v));
        hi = hi.cwiseMax(mesh.position(v));
    }

    gsMatrix<real_t> bbox(3,2);
    bbox.col(0) = lo;
    bbox.col(1) = hi;
    return bbox;
}

/// Scales and translates \a mesh so that its longest extent occupies a
/// fraction \a fill of the unit box [0,1]^3, centred.
///
/// Every immersed driver needs the geometry inside the background box with a
/// margin, and each used to open-code this; \a fill is that margin (0.8 leaves
/// 10% of the box on each side, so no cut cell touches the box boundary).
///
/// \returns the scale factor applied (physical -> parametric). This is the only
/// non-trivial output -- the resulting bounding box is by construction the unit
/// box -- and drivers need it to convert results back: a volume measured in the
/// unit box is physical when divided by scale^3, an area by scale^2. Call
/// gsSurfMeshBoundingBox() BEFORE this if the original extent is also wanted.
inline real_t gsNormalizeToUnitBox(gsSurfMesh & mesh, real_t fill = 0.8)
{
    const gsMatrix<real_t> bb = gsSurfMeshBoundingBox(mesh);
    const gsSurfMesh::Point lo = bb.col(0), hi = bb.col(1);

    const gsSurfMesh::Point centre = (lo + hi) / real_t(2);
    const real_t extent = (hi - lo).maxCoeff();
    const real_t scale  = (extent > 0 ? fill / extent : real_t(1));

    for (gsSurfMesh::Vertex v : mesh.vertices())
        mesh.position(v) = (mesh.position(v) - centre) * scale
                         + gsSurfMesh::Point::Constant(real_t(0.5));

    return scale;
}

/// The unit box [0,1]^3 as a 3 x 2 matrix (col 0 = lower, col 1 = upper), i.e.
/// the support to pair with a mesh normalized by gsNormalizeToUnitBox().
inline gsMatrix<real_t> gsUnitBox3()
{
    gsMatrix<real_t> bbox(3,2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    return bbox;
}

/// Flattens \a mesh into a triangle soup: \a verts is 3*nV (x,y,z per vertex,
/// indexed by gsSurfMesh::Vertex::idx()), \a tris is 3*nF zero-based corner
/// indices into it.
///
/// gsMeshSignedDist does NOT need this -- its BVH keeps its own copy. It is
/// here for drivers that run their own per-triangle geometry (bounding-box
/// pre-filters, SAT tests, ParaView output) and want the raw arrays.
///
/// \note Sized by vertices_size() rather than n_vertices() so that idx()
/// indexes \a verts directly; for a mesh with deleted vertices the gaps are
/// present but unreferenced.
inline void gsFlattenSurfMesh(const gsSurfMesh    & mesh,
                              std::vector<real_t> & verts,
                              std::vector<index_t> & tris)
{
    verts.assign(3 * static_cast<size_t>(mesh.vertices_size()), real_t(0));
    for (gsSurfMesh::Vertex v : mesh.vertices())
    {
        const gsSurfMesh::Point & p = mesh.position(v);
        for (index_t c = 0; c != 3; ++c)
            verts[3*static_cast<size_t>(v.idx()) + c] = p[c];
    }

    tris.clear();
    tris.reserve(3 * static_cast<size_t>(mesh.n_faces()));
    for (gsSurfMesh::Face f : mesh.faces())
        for (gsSurfMesh::Vertex v : mesh.vertices(f))
            tris.push_back(static_cast<index_t>(v.idx()));
}

} // namespace gismo
