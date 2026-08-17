/** @file gsSurfMeshBVH.h

    @brief A static bounding-volume hierarchy (BVH) over the triangles of a
           gsSurfMesh, answering the two queries a mesh signed-distance level
           set needs:
             - nearest point on the surface to a query point (for |phi|),
             - generalized winding number (for sign(phi)).

    Without this structure both queries are O(nTriangles) per point: every
    evaluation loops over the whole triangle list. For a mesh of a few thousand
    triangles queried at every quadrature/octree-classification point of an FCM
    assembly, this dominates the runtime by orders of magnitude (profiled:
    ~400e6 triangle tests for a single small assembly pass). See
    gsDomain/gsMeshLevelSet.h for the consumer.

    Nearest point: a standard AABB-tree branch-and-bound search, pruning
    subtrees whose box lower-bound distance already exceeds the current best.

    Winding number: the "fast winding number" method (Barill, Dickson,
    Schmidt, Levin, Jacobson, 2018). Each internal node stores the
    area-weighted centroid and the summed area-weighted normal ("dipole
    moment") of the triangles in its subtree, plus a radius bounding how far
    those triangles are scattered from the centroid. For a query point far
    enough from a node (dist > beta * radius), the node's contribution to the
    winding number is approximated in O(1) via the dipole formula instead of
    summing every triangle; near-field triangles (which matter most, since the
    level set is queried mostly near/inside cut cells close to the surface)
    still get exact per-triangle solid angles.

    WHY THIS EXISTS RATHER THAN A THIRD-PARTY DEPENDENCY. The query is
    "exact distance and inside/outside against a raw triangle soup", which is a
    geometry-processing primitive, not a meshing one:
      - gmsh has no such API. Its getElementByCoordinates() is CONTAINMENT
        backed by an octree, and on a miss it falls back to a brute-force sweep
        over every element while mutating a process-global tolerance -- unusable
        from the OpenMP cut-cell loops. Its Distance FIELD replaces surfaces by
        a SAMPLED point cloud (nanoflann kd-tree over points), so it answers a
        different, approximate question. Its real BVH (contrib/hxt) and its
        point-to-surface projector (SurfaceProjector, libOL) are internal and
        not exposed through the public API.
      - CGAL's AABB_tree + Side_of_triangle_mesh match exactly, but are
        GPL-3.0-or-later, incompatible with G+Smo's MPL 2.0.
      - libigl (MPL 2.0) does match: igl::AABB::squared_distance and
        igl::FastWindingNumberForSoups implement the same two queries with the
        same method. Adopting it remains a live option; see the note on
        windingOmega() for where its implementation is more advanced than this
        one (a second-order far-field expansion, where this uses first order).
        Its igl::pseudonormal_test is what the sharp-feature handling here was
        modelled on -- though this file derives the closest FEATURE directly
        from the Ericson branch rather than reconstructing it from barycentric
        coordinates, which removes libigl's tolerance heuristics for slivers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsSurfMesh.h>

#include <algorithm>
#include <array>
#include <map>
#include <limits>
#include <vector>

namespace gismo
{

/// A static AABB-tree over the triangles of a gsSurfMesh, answering
/// nearest-point and generalized-winding-number queries in O(log n) amortized
/// instead of the brute-force O(n).
///
/// The tree COPIES the triangle geometry it needs at build() time, so it stays
/// valid independently of the source mesh -- required because its owner
/// (gsMeshSignedDist) is copy-constructed by GISMO_CLONE_FUNCTION, and
/// convenient because the copy is where the outward-orientation normalization
/// is applied without touching the caller's mesh.
class gsSurfMeshBVH
{
public:
    typedef gsSurfMesh::Point Point;   ///< gsVector3d<real_t>

    gsSurfMeshBVH() : m_root(-1) {}

    /// Builds the tree from the triangles of \a mesh.
    ///
    /// \a mesh must be a triangle mesh; call gsSurfMesh::triangulate() first if
    /// it is not. The triangles are copied and, per connected component,
    /// normalized to OUTWARD winding (see _extractOriented) -- everything
    /// downstream may therefore assume triangleNormal() points out of the solid.
    void build(const gsSurfMesh & mesh)
    {
        m_nodes.clear();
        m_orderedTris.clear();
        m_orderedNrm.clear();
        m_root = -1;

        GISMO_ENSURE(mesh.is_triangle_mesh(),
                     "gsSurfMeshBVH: expects a triangle mesh; call "
                     "gsSurfMesh::triangulate() before building the tree.");

        std::vector<Point>   tri;    // 3 points per triangle, outward-wound
        std::vector<index_t> triV;   // matching corner vertex indices
        _extractOriented(mesh, tri, triV);

        const size_t nt = tri.size() / 3;
        if (nt == 0) return;

        std::vector<TriNormals> tn;
        _computeNormals(tri, triV, mesh.vertices_size(), tn);

        m_nodes.reserve(2 * nt);
        m_orderedTris.reserve(3 * nt);
        m_orderedNrm.reserve(nt);

        std::vector<TriInfo> infos(nt);
        for (size_t t = 0; t != nt; ++t)
        {
            const Point & a = tri[3*t+0];
            const Point & b = tri[3*t+1];
            const Point & c = tri[3*t+2];
            TriInfo & ti = infos[t];
            ti.a = a; ti.b = b; ti.c = c;
            // direction = normal, magnitude = area
            ti.normalArea = (b-a).cross(c-a) * real_t(0.5);
            ti.centroid   = (a+b+c) / real_t(3);
            ti.bmin = a.cwiseMin(b).cwiseMin(c);
            ti.bmax = a.cwiseMax(b).cwiseMax(c);
            ti.nrm  = tn[t];
        }

        std::vector<index_t> order(nt);
        for (size_t t = 0; t != nt; ++t) order[t] = static_cast<index_t>(t);

        m_root = _buildRecursive(infos, order, 0, static_cast<index_t>(nt));
    }

    bool empty() const { return m_root < 0; }

    /// Number of triangles in the tree.
    index_t numTriangles() const
    { return static_cast<index_t>(m_orderedTris.size() / 3); }

    /// Nearest point on the mesh to \a p; returns the SQUARED distance and
    /// writes the closest point into \a closest.
    ///
    /// \a nearestTri, if non-null, receives the index of the winning triangle
    /// in this tree's internal (leaf-contiguous) ordering -- feed it to
    /// triangleNormal(), not to any index of the source mesh. It is left at -1
    /// for an empty tree.
    /// \a feature, if non-null, receives which feature of that triangle the
    /// closest point lies on (see Feature); pass it to pseudoNormal()/isSharp().
    real_t squaredDistance(const Point & p, Point & closest,
                           index_t * nearestTri = nullptr,
                           short_t * feature    = nullptr) const
    {
        real_t  best = std::numeric_limits<real_t>::max();
        index_t tri  = -1;
        short_t feat = FaceInterior;
        closest = p;
        if (!empty())
            _queryNearest(m_root, p, best, closest, tri, feat);
        if (nearestTri) *nearestTri = tri;
        if (feature)    *feature    = feat;
        return best;
    }

    /// Outward pseudonormal of \a feature on triangle \a i: the angle-weighted
    /// vertex normal, the edge normal, or the face normal, as appropriate.
    ///
    /// On a vertex or an edge this is a CONVENTION, not a limit: grad(phi) has
    /// no limit there (approaching along different directions of the normal
    /// cone gives different answers). It is the classical choice, and the one
    /// for which Baerentzen & Aanaes proved the sign test correct.
    Point pseudoNormal(index_t i, short_t feature) const
    {
        if (i < 0 || i >= static_cast<index_t>(m_orderedNrm.size()))
            return Point::Zero();
        if (feature >= VertexA && feature <= VertexC)
            return m_orderedNrm[i].vn[feature];
        if (feature >= EdgeAB && feature <= EdgeBC)
            return m_orderedNrm[i].en[feature - EdgeAB];
        return triangleNormal(i);
    }

    /// True when a plain dot-product-with-the-normal test is NOT a safe sign
    /// source at \a feature of triangle \a i, i.e. when the incident face
    /// normals span more than 90 degrees (interior dihedral outside
    /// [90, 270] degrees). Callers should fall back to a winding-number query
    /// there; see gsMeshSignedDist::deriv_into().
    ///
    /// Always false for a face interior, where (p - closest) is exactly
    /// +/- the face normal.
    bool isSharp(index_t i, short_t feature) const
    {
        if (i < 0 || i >= static_cast<index_t>(m_orderedNrm.size())) return false;
        if (feature < VertexA || feature > EdgeBC) return false;
        return m_orderedNrm[i].sharp[feature] != 0;
    }

    /// Number of (triangle, feature) pairs flagged sharp. Zero means the fast
    /// sign path is provably correct for every query against this mesh.
    index_t numSharpFeatures() const
    {
        index_t n = 0;
        for (size_t t = 0; t != m_orderedNrm.size(); ++t)
            for (index_t k = 0; k != 6; ++k)
                if (m_orderedNrm[t].sharp[k]) ++n;
        return n;
    }

    /// Unit normal of triangle \a i in this tree's internal ordering (as
    /// reported by squaredDistance()), from its vertex winding order.
    ///
    /// build() normalizes the winding per component, so for a closed mesh this
    /// points out of the solid -- which is what lets it stand in for grad(phi)
    /// on {phi==0}.
    ///
    /// \warning This is the FACE normal, which is only the right answer when
    /// the closest point lies in the INTERIOR of a face. On a shared edge or
    /// vertex grad(phi) has no limit at all, and the conventional choice is the
    /// angle-weighted pseudonormal -- use pseudoNormal(i, feature), which
    /// reduces to this for a face interior. See gsMeshLevelSet.h::deriv_into().
    Point triangleNormal(index_t i) const
    {
        if (i < 0 || 3*i + 2 >= static_cast<index_t>(m_orderedTris.size()))
            return Point::Zero();
        const Point & a = m_orderedTris[3*i+0];
        const Point & b = m_orderedTris[3*i+1];
        const Point & c = m_orderedTris[3*i+2];
        const Point n = (b-a).cross(c-a);
        const real_t L = n.norm();
        return (L > 0 ? Point(n/L) : Point::Zero());
    }

    /// 4*pi*(generalized winding number) at \a p, i.e. what a brute-force sum
    /// of the per-triangle signed solid angles would give.
    ///
    /// \a beta is the Barnes-Hut acceptance parameter: a subtree is replaced by
    /// its dipole approximation once dist(p, node) > beta * node.radius.
    ///
    /// Only the SIGN of the winding number is consumed downstream
    /// (gsMeshSignedDist thresholds |omega/4pi| at 0.5), which tolerates far
    /// more far-field error than the winding value itself, so beta can be
    /// pushed well below the conservative 2.0 this started at. Measured on
    /// spot.obj (5856 triangles, 200k query points jittered within one h=1/16
    /// cut cell of the surface -- the distribution that dominates an FCM
    /// assembly), against a brute-force solid-angle sum as ground truth:
    ///
    ///   beta   cost/point   inside/outside verdicts wrong
    ///   2.0    4.26 us      0 / 200000
    ///   1.5    2.21 us      0 / 200000
    ///   1.2    1.43 us      0 / 200000   <-- default
    ///   1.0    0.88 us      7 / 200000   <-- accuracy cliff starts here
    ///   0.8    0.33 us   24989 / 200000
    ///
    /// 1.2 is therefore 3x cheaper than 2.0 at identical classification, with
    /// the cliff a clear factor away. Since the nearest-point query costs
    /// ~1.36 us/point, this halves the cost of a whole signed-distance call.
    /// Do not lower it without re-running that comparison: below 1.0 the
    /// far-field error crosses the 0.5 threshold and silently inverts
    /// inside/outside classification near the surface. (Timings were taken in
    /// a real_t = double build; beta is a GEOMETRIC criterion, so the verdict
    /// counts carry over to other precisions and only the cost column scales.)
    ///
    /// \note libigl's fast_winding_number defaults to beta = 2.0 but with a
    /// SECOND-order (quadrupole) expansion, where this uses first order
    /// (dipole) only. Higher order buys far-field accuracy at equal beta --
    /// which is worth having if the winding VALUE is ever needed, but buys
    /// nothing for the sign test above, where beta = 1.2 already misclassifies
    /// nothing. Raise the order before raising beta if that ever changes.
    real_t windingOmega(const Point & p, real_t beta = real_t(1.2)) const
    {
        return empty() ? real_t(0) : _queryWinding(m_root, p, beta);
    }

    /// Which feature of a triangle the closest point landed on. The numbering
    /// is the storage order used by pseudoNormal()/isSharp().
    enum Feature
    {
        VertexA = 0, VertexB = 1, VertexC = 2,   ///< corners 0, 1, 2
        EdgeAB  = 3, EdgeAC  = 4, EdgeBC  = 5,   ///< corners (0,1), (0,2), (1,2)
        FaceInterior = 6
    };

    /// Closest point on triangle (a,b,c) to \a p (Ericson, Real-Time Collision
    /// Detection). Exposed because callers occasionally need it standalone.
    ///
    /// \a feature, if non-null, receives which feature the result lies on. The
    /// routine already branches on exactly that, so this costs nothing and --
    /// unlike reconstructing it afterwards from barycentric coordinates, as
    /// libigl's pseudonormal_test does -- needs no tolerance and no special
    /// case for sliver triangles. Note also that the three edge branches divide
    /// by |ab|^2, |ac|^2 and |bc|^2 respectively (e.g. d1-d3 == |ab|^2), so they
    /// are division-safe for any triangle with non-zero edges.
    ///
    /// The branch guards are inclusive (<=, >=), so a point exactly on an edge
    /// is reported as that edge rather than as the face interior. Two triangles
    /// sharing an edge therefore agree on the feature class of a point on it.
    static Point closestPointOnTriangle(const Point & p, const Point & a,
                                        const Point & b, const Point & c,
                                        short_t * feature = nullptr)
    {
        // Keep the assignments adjacent to their returns: the source order is
        // A, B, AB, C, AC, BC, interior -- NOT A, B, C, AB, AC, BC.
        const Point ab = b-a, ac = c-a, ap = p-a;
        const real_t d1 = ab.dot(ap), d2 = ac.dot(ap);
        if (d1 <= 0 && d2 <= 0)
        { if (feature) *feature = VertexA; return a; }

        const Point bp = p-b;
        const real_t d3 = ab.dot(bp), d4 = ac.dot(bp);
        if (d3 >= 0 && d4 <= d3)
        { if (feature) *feature = VertexB; return b; }

        const real_t vc = d1*d4 - d3*d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0)
        { if (feature) *feature = EdgeAB; return a + ab * (d1 / (d1-d3)); }

        const Point cp = p-c;
        const real_t d5 = ab.dot(cp), d6 = ac.dot(cp);
        if (d6 >= 0 && d5 <= d6)
        { if (feature) *feature = VertexC; return c; }

        const real_t vb = d5*d2 - d1*d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0)
        { if (feature) *feature = EdgeAC; return a + ac * (d2 / (d2-d6)); }

        const real_t va = d3*d6 - d5*d4;
        if (va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0)
        { if (feature) *feature = EdgeBC;
          return b + (c-b) * ((d4-d3) / ((d4-d3) + (d5-d6))); }

        if (feature) *feature = FaceInterior;
        const real_t denom = real_t(1) / (va+vb+vc);
        return a + ab*(vb*denom) + ac*(vc*denom);
    }

    /// Signed solid angle subtended by triangle (a,b,c) at \a p (Van Oosterom
    /// & Strackee). Summed over a closed mesh this equals 4*pi*(winding number).
    static real_t signedSolidAngle(const Point & p, const Point & a,
                                   const Point & b, const Point & c)
    {
        const Point va = a-p, vb = b-p, vc = c-p;
        const real_t la = va.norm(), lb = vb.norm(), lc = vc.norm();
        const real_t numer = va.dot(vb.cross(vc));
        const real_t denom = la*lb*lc
                           + va.dot(vb)*lc
                           + vb.dot(vc)*la
                           + vc.dot(va)*lb;
        return real_t(2) * math::atan2(numer, denom);
    }

private:
    /// Angle-weighted vertex and uniform edge pseudonormals, plus a per-feature
    /// flag saying whether the plain face normal would be an unsafe sign
    /// source there.
    ///
    /// WHY. `deriv_into` takes the sign of phi from `dot(p - closest, n)`. For
    /// a point whose closest point lies on an EDGE, (p - closest) lies in the
    /// normal cone spanned by the two incident face normals n1, n2; writing it
    /// a*n1 + b*n2 with a,b >= 0,
    ///
    ///     dot(p - closest, n1) = a + b*cos(angle(n1,n2))
    ///
    /// which can only go negative when the two normals are more than 90 degrees
    /// apart, i.e. when the interior dihedral angle theta is outside
    /// [90, 270] degrees. So the face normal is a PROVABLY correct sign source
    /// on any feature that is not "sharp" in that sense, and only there.
    ///
    /// The angle-weighted pseudonormal (Baerentzen & Aanaes) is the classical
    /// repair, and it is what this stores -- but its correctness theorem
    /// assumes a closed, manifold, consistently oriented mesh, which .obj files
    /// routinely are not. So sharpness is recorded too, and the caller is
    /// expected to fall back to a winding-number query on sharp features
    /// (see gsMeshSignedDist::deriv_into). The pseudonormal is then used for
    /// the gradient DIRECTION, where phi is non-differentiable anyway and the
    /// bisector is a convention rather than a limit.
    ///
    /// A face interior never needs either: there (p - closest) is exactly
    /// +/- the face normal, so the dot product cannot have the wrong sign.
    struct TriNormals
    {
        Point vn[3];        ///< angle-weighted pseudonormal at corners 0,1,2
        Point en[3];        ///< edge pseudonormal for (0,1), (0,2), (1,2)
        char  sharp[6];     ///< per feature 0..5: face normal unsafe as a sign source
    };

    /// True when two unit normals are more than 90 degrees apart, i.e. when a
    /// point in the cone they span can have a negative dot product with one of
    /// them. Uses a small negative tolerance so that the exactly-90-degree case
    /// (a cube edge, an L-shaped solid) stays on the SAFE side: there the worst
    /// dot product is exactly zero, and `dot < 0` is false, so the fast path is
    /// still correct.
    static bool _normalsSpanTooWide(const Point & n1, const Point & n2)
    {
        return n1.dot(n2) < -1e-12;
    }

    struct TriInfo
    {
        Point a, b, c;
        Point normalArea;  // 0.5*cross(b-a,c-a): direction = normal, magnitude = area
        Point centroid;
        Point bmin, bmax;
        TriNormals nrm;    // pseudonormals + sharpness, see TriNormals
    };

    struct Node
    {
        Point   bmin, bmax;         // AABB of the subtree (nearest-point query)
        Point   centroid;           // area-weighted centroid (winding dipole)
        Point   normalArea;         // sum of triangle (area * normal) (winding dipole)
        real_t  totalArea;          // sum of triangle areas, for merging centroids
        real_t  radius;             // max |vertex - centroid| over the subtree
        index_t left, right;        // child node indices, -1 if leaf
        index_t triStart, triCount; // valid range into m_orderedTris if leaf
    };

    static const index_t LEAF_SIZE = 8;

    /// Fills \a tn (one entry per triangle) with the angle-weighted vertex and
    /// uniform edge pseudonormals and the per-feature sharpness flags.
    ///
    /// Runs AFTER _extractOriented(), on the already-outward-wound \a tri and
    /// its matching corner indices \a triV, so every face normal here is the
    /// outward one and no flip bookkeeping is needed.
    ///
    /// Degeneracy is handled by omission rather than by tolerance: a triangle
    /// with a zero-length cross product contributes nothing to any accumulator,
    /// and a corner whose two incident edges have zero length contributes no
    /// angle. That matters because a single NaN reaching a pseudonormal would
    /// make `dot(d, n) < 0` false and silently pin the sign to "outside" at
    /// every query landing on that feature.
    static void _computeNormals(const std::vector<Point>   & tri,
                                const std::vector<index_t> & triV,
                                size_t nVert,
                                std::vector<TriNormals>    & tn)
    {
        const size_t nt = tri.size() / 3;
        tn.assign(nt, TriNormals());

        // WELD coincident vertex INDICES (positions and topology untouched).
        //
        // The accumulators below are keyed on mesh vertices and mesh edges, so
        // a file that stores geometrically coincident but topologically
        // distinct vertices -- routine in .obj/STL exported per face -- splits
        // each accumulator between the duplicates. The shared edge is then seen
        // as two boundary edges, each receiving a single incident face, and the
        // "edge pseudonormal" silently collapses back to that face's normal:
        // the whole mechanism no-ops at exactly the features it exists for.
        //
        // Matching is EXACT, deliberately. A tolerance would be the thing that
        // eats legitimate thin geometry -- the 1e-12 corner shave in
        // sliver.obj, or a genuinely sharp fin -- and per-face exporters emit
        // bit-identical coordinates for shared corners, so exact matching
        // catches the real case with no false positives. (Comparison-based
        // ordering also treats -0.0 and +0.0 as one key, which is what we want.)
        //
        // This is a RELABELLING for accumulation only: `tri`, the BVH tree and
        // the reported triangle indices are unaffected.
        std::vector<index_t> canon(nVert);
        for (size_t i = 0; i != nVert; ++i) canon[i] = static_cast<index_t>(i);
        {
            typedef std::array<real_t,3> Key3;
            std::map<Key3, index_t> seen;
            for (size_t k = 0; k != triV.size(); ++k)
            {
                const Key3 key = {{ tri[k][0], tri[k][1], tri[k][2] }};
                const std::pair<typename std::map<Key3,index_t>::iterator,bool>
                    ins = seen.insert(std::make_pair(key, triV[k]));
                canon[triV[k]] = ins.first->second;
            }
        }

        // Per-triangle outward face normal (zero for a degenerate triangle).
        std::vector<Point> fn(nt, Point::Zero());
        for (size_t t = 0; t != nt; ++t)
        {
            const Point n = (tri[3*t+1]-tri[3*t+0]).cross(tri[3*t+2]-tri[3*t+0]);
            const real_t L = n.norm();
            if (L > 0) fn[t] = n / L;
        }

        // Accumulate angle-weighted normals per MESH vertex, and the incident
        // face list needed for the vertex sharpness test.
        std::vector<Point>                vAcc(nVert, Point::Zero());
        std::vector<std::vector<index_t> > vFaces(nVert);
        for (size_t t = 0; t != nt; ++t)
        {
            if (fn[t].isZero()) continue;
            for (index_t k = 0; k != 3; ++k)
            {
                const Point & P = tri[3*t + k];
                const Point   e1 = tri[3*t + (k+1)%3] - P;
                const Point   e2 = tri[3*t + (k+2)%3] - P;
                const real_t  l1 = e1.norm(), l2 = e2.norm();
                if (l1 <= 0 || l2 <= 0) continue;          // degenerate corner
                const real_t cs = math::max(real_t(-1),
                                  math::min(real_t(1), e1.dot(e2)/(l1*l2)));
                const index_t v = canon[triV[3*t + k]];
                vAcc  [v] += fn[t] * math::acos(cs);       // ANGLE weighting
                vFaces[v].push_back(static_cast<index_t>(t));
            }
        }

        // Accumulate per EDGE, keyed on the unordered corner-vertex pair. Two
        // triangles sharing an edge hit the same key and therefore end up with
        // a bitwise identical edge pseudonormal.
        typedef std::pair<index_t,index_t> Key;
        std::map<Key, std::pair<Point,index_t> > eAcc;   // key -> (sum n, count)
        std::map<Key, bool>                      eSharp;
        // Ericson feature order: edge 0 = corners (0,1), 1 = (0,2), 2 = (1,2).
        static const index_t EC[3][2] = {{0,1},{0,2},{1,2}};
        for (size_t t = 0; t != nt; ++t)
        {
            if (fn[t].isZero()) continue;
            for (index_t e = 0; e != 3; ++e)
            {
                const index_t u = canon[triV[3*t + EC[e][0]]],
                              w = canon[triV[3*t + EC[e][1]]];
                if (u == w) continue;      // collapsed by welding: not an edge
                const Key key(math::min(u,w), math::max(u,w));
                std::pair<Point,index_t> & acc = eAcc[key];
                if (acc.second == 0) acc.first = Point::Zero();
                // Sharp as soon as ANY previously seen incident face normal is
                // more than 90 degrees from this one.
                if (acc.second != 0 && _normalsSpanTooWide(acc.first, fn[t]))
                    eSharp[key] = true;
                acc.first += fn[t];
                ++acc.second;
            }
        }

        // Scatter back per triangle.
        for (size_t t = 0; t != nt; ++t)
        {
            TriNormals & r = tn[t];
            for (index_t k = 0; k != 6; ++k) r.sharp[k] = 0;

            for (index_t k = 0; k != 3; ++k)
            {
                const index_t v = canon[triV[3*t + k]];
                const real_t  L = vAcc[v].norm();
                r.vn[k] = (L > 0 ? Point(vAcc[v]/L) : fn[t]);
                // Nothing usable accumulated here (every incident triangle was
                // degenerate, or they cancelled exactly). A zero normal is as
                // dangerous as a NaN one: `d.dot(n) < 0` is then false and the
                // sign pins to "outside". Route the feature to the winding
                // fallback instead of trusting it.
                if (r.vn[k].isZero()) { r.sharp[k] = 1; continue; }

                // A vertex is unsafe if ANY pair of its incident face normals
                // spans more than 90 degrees -- adjacent-pair checks are not
                // enough (a faceted cone can step round in safe increments and
                // still have opposite facets pointing apart).
                const std::vector<index_t> & inc = vFaces[v];
                bool sharp = false;
                for (size_t i = 0; i < inc.size() && !sharp; ++i)
                    for (size_t j = i+1; j < inc.size() && !sharp; ++j)
                        sharp = _normalsSpanTooWide(fn[inc[i]], fn[inc[j]]);
                r.sharp[k] = sharp ? 1 : 0;
            }

            for (index_t e = 0; e != 3; ++e)
            {
                const index_t u = canon[triV[3*t + EC[e][0]]],
                              w = canon[triV[3*t + EC[e][1]]];
                if (u == w) { r.en[e] = fn[t]; r.sharp[3+e] = 1; continue; }
                const Key key(math::min(u,w), math::max(u,w));
                const std::pair<Point,index_t> & acc = eAcc[key];
                const real_t L = acc.first.norm();
                r.en[e] = (L > 0 ? Point(acc.first/L) : fn[t]);
                // Same reasoning as for vertices above: an unusable edge
                // normal must fall back to the winding number, not be trusted.
                r.sharp[3+e] = (eSharp.count(key) || r.en[e].isZero()) ? 1 : 0;
            }
        }
    }

    /// Copies the mesh triangles into \a tri (3 points per triangle) with
    /// OUTWARD winding per connected component.
    ///
    /// Only the GLOBAL inversion has to be repaired here. Local inconsistency
    /// -- neighbouring faces disagreeing -- cannot occur: a halfedge mesh is
    /// structurally incapable of representing it, and gsSurfMesh::add_face()
    /// rejects such a face outright ("complex edge") rather than storing it.
    /// What remains is a component that is consistently wound but inside-out,
    /// which the signed volume (1/6)*sum dot(a, b x c) detects: it is
    /// translation-invariant for a closed surface and negative exactly when the
    /// component is inverted.
    ///
    /// The check is skipped for a component touching a boundary edge, since an
    /// open surface encloses no well-defined volume; that case is warned about,
    /// because a signed distance is only meaningful against a closed surface.
    ///
    /// O(nFaces) via one BFS over face adjacency.
    static void _extractOriented(const gsSurfMesh & mesh, std::vector<Point> & tri,
                                 std::vector<index_t> & triV)
    {
        const index_t nf = static_cast<index_t>(mesh.n_faces());
        tri.clear();  triV.clear();
        tri.reserve(3 * static_cast<size_t>(nf));
        triV.reserve(3 * static_cast<size_t>(nf));
        if (nf == 0) return;

        // Triangles in face order, plus a face-index map for the BFS below.
        // The CORNER VERTEX INDICES are tracked alongside the positions and
        // swapped with them by the flip below, so downstream code never has to
        // ask whether a component was inverted -- triV always describes the
        // winding actually stored in tri. It also means the per-edge lookup in
        // _computeNormals() is keyed on real vertex pairs rather than on
        // halfedge slot arithmetic, which is where the local edge ordering
        // (halfedge(f) gives edges CA, AB, BC, not AB, AC, BC) would otherwise
        // silently mismatch the Ericson feature numbering.
        std::vector<index_t> fidx;
        fidx.reserve(static_cast<size_t>(nf));
        for (gsSurfMesh::Face f : mesh.faces())
        {
            for (gsSurfMesh::Vertex v : mesh.vertices(f))
            { tri.push_back(mesh.position(v)); triV.push_back(v.idx()); }
            fidx.push_back(f.idx());
        }

        // Component labelling by BFS over halfedge adjacency, and per-component
        // signed volume. Faces are addressed by their POSITION in tri (0..nf),
        // so map mesh face indices onto that.
        std::vector<index_t> slot(mesh.faces_size(), -1);
        for (size_t i = 0; i != fidx.size(); ++i)
            slot[fidx[i]] = static_cast<index_t>(i);

        std::vector<char>    visited(nf, 0);
        std::vector<index_t> comp, queue;
        index_t nOpen = 0;

        for (index_t seed = 0; seed != nf; ++seed)
        {
            if (visited[seed]) continue;
            comp.clear(); queue.clear();
            queue.push_back(seed); visited[seed] = 1;
            bool closed = true;

            for (size_t qi = 0; qi != queue.size(); ++qi)
            {
                const index_t s = queue[qi];
                comp.push_back(s);
                const gsSurfMesh::Face f(fidx[s]);
                for (gsSurfMesh::Halfedge h : mesh.halfedges(f))
                {
                    const gsSurfMesh::Halfedge o = mesh.opposite_halfedge(h);
                    const gsSurfMesh::Face     g = mesh.face(o);
                    if (!g.is_valid()) { closed = false; continue; } // boundary edge
                    const index_t t = slot[g.idx()];
                    if (t < 0 || visited[t]) continue;
                    visited[t] = 1;
                    queue.push_back(t);
                }
            }

            if (!closed) { ++nOpen; continue; } // no enclosed volume to sign

            real_t vol6 = 0;
            for (size_t i = 0; i != comp.size(); ++i)
            {
                const index_t t = comp[i];
                vol6 += tri[3*t+0].dot(tri[3*t+1].cross(tri[3*t+2]));
            }
            if (vol6 < 0)
                for (size_t i = 0; i != comp.size(); ++i)
                {
                    std::swap(tri [3*comp[i]+1], tri [3*comp[i]+2]);
                    std::swap(triV[3*comp[i]+1], triV[3*comp[i]+2]);
                }
        }

        if (nOpen)
            gsWarn << "gsSurfMeshBVH: " << nOpen << " connected component(s) of "
                      "the mesh have boundary edges, i.e. the surface is not "
                      "closed. Outward orientation could not be verified there, "
                      "and a signed distance against an open surface is not "
                      "well defined -- expect a wrong sign near those parts.\n";
    }

    /// Builds the subtree over infos[order[begin..end)], returns its node index.
    index_t _buildRecursive(std::vector<TriInfo> & infos, std::vector<index_t> & order,
                            index_t begin, index_t end)
    {
        const index_t nodeIdx = static_cast<index_t>(m_nodes.size());
        m_nodes.emplace_back();

        // Bounding box + area-weighted moments of this range.
        Point  bmin = Point::Constant( std::numeric_limits<real_t>::max());
        Point  bmax = Point::Constant(-std::numeric_limits<real_t>::max());
        Point  normalArea      = Point::Zero();
        Point  areaCentroidSum = Point::Zero();
        real_t totalArea = 0;
        for (index_t i = begin; i != end; ++i)
        {
            const TriInfo & ti = infos[order[i]];
            bmin = bmin.cwiseMin(ti.bmin);
            bmax = bmax.cwiseMax(ti.bmax);
            const real_t area = ti.normalArea.norm();
            normalArea      += ti.normalArea;
            areaCentroidSum += ti.centroid * area;
            totalArea       += area;
        }
        // Degenerate (zero-area) ranges: fall back to the unweighted centroid
        // so windingOmega() never divides by zero.
        Point centroid;
        if (totalArea > 0)
            centroid = areaCentroidSum / totalArea;
        else
        {
            centroid = Point::Zero();
            for (index_t i = begin; i != end; ++i)
                centroid += infos[order[i]].centroid;
            centroid /= real_t(end - begin);
        }

        real_t radius = 0;
        for (index_t i = begin; i != end; ++i)
        {
            const TriInfo & ti = infos[order[i]];
            radius = math::max(radius, (ti.a - centroid).norm());
            radius = math::max(radius, (ti.b - centroid).norm());
            radius = math::max(radius, (ti.c - centroid).norm());
        }

        const index_t count = end - begin;
        if (count <= LEAF_SIZE)
        {
            Node & n = m_nodes[nodeIdx]; // re-fetch: emplace_back may have reallocated
            n.bmin = bmin; n.bmax = bmax;
            n.centroid = centroid; n.normalArea = normalArea; n.totalArea = totalArea;
            n.radius = radius;
            n.left = n.right = -1;
            n.triStart = static_cast<index_t>(m_orderedTris.size() / 3);
            n.triCount = count;
            for (index_t i = begin; i != end; ++i)
            {
                const TriInfo & ti = infos[order[i]];
                m_orderedTris.push_back(ti.a);
                m_orderedTris.push_back(ti.b);
                m_orderedTris.push_back(ti.c);
                // MUST stay in lockstep with m_orderedTris: _buildRecursive
                // permutes `order`, so triangle k of the tree is not triangle k
                // of the mesh and a missing push here silently mis-indexes
                // every pseudonormal.
                m_orderedNrm.push_back(ti.nrm);
            }
            return nodeIdx;
        }

        // Split on the longest axis of the box, at the median centroid.
        short_t axis = 0;
        (bmax - bmin).maxCoeff(&axis);

        const index_t mid = begin + count / 2;
        std::nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
            [&](index_t i1, index_t i2)
            { return infos[i1].centroid[axis] < infos[i2].centroid[axis]; });

        const index_t left  = _buildRecursive(infos, order, begin, mid);
        const index_t right = _buildRecursive(infos, order, mid, end);

        Node & n = m_nodes[nodeIdx]; // re-fetch: recursive calls reallocated m_nodes
        n.bmin = bmin; n.bmax = bmax;
        n.centroid = centroid; n.normalArea = normalArea; n.totalArea = totalArea;
        n.radius = radius;
        n.left = left; n.right = right;
        n.triStart = n.triCount = 0;
        return nodeIdx;
    }

    static real_t _squaredDistPointAABB(const Point & p,
                                        const Point & bmin, const Point & bmax)
    {
        // Per component: 0 inside the slab, the overshoot outside it.
        const Point d = (bmin - p).cwiseMax(Point::Zero())
                      + (p - bmax).cwiseMax(Point::Zero());
        return d.squaredNorm();
    }

    void _queryNearest(index_t nodeIdx, const Point & p, real_t & best,
                       Point & closest, index_t & bestTri, short_t & bestFeat) const
    {
        const Node & n = m_nodes[nodeIdx];
        if (_squaredDistPointAABB(p, n.bmin, n.bmax) >= best) return; // pruned

        if (n.left < 0) // leaf
        {
            for (index_t i = n.triStart; i != n.triStart + n.triCount; ++i)
            {
                short_t feat = FaceInterior;
                const Point q = closestPointOnTriangle(p, m_orderedTris[3*i+0],
                                                          m_orderedTris[3*i+1],
                                                          m_orderedTris[3*i+2],
                                                          &feat);
                const real_t d2 = (p-q).squaredNorm();
                if (d2 < best)
                { best = d2; closest = q; bestTri = i; bestFeat = feat; }
            }
            return;
        }

        // Visit the nearer child first so the sibling is more likely pruned.
        const Node & L = m_nodes[n.left];
        const Node & R = m_nodes[n.right];
        const real_t dL = _squaredDistPointAABB(p, L.bmin, L.bmax);
        const real_t dR = _squaredDistPointAABB(p, R.bmin, R.bmax);
        if (dL <= dR)
        {
            if (dL < best) _queryNearest(n.left,  p, best, closest, bestTri, bestFeat);
            if (dR < best) _queryNearest(n.right, p, best, closest, bestTri, bestFeat);
        }
        else
        {
            if (dR < best) _queryNearest(n.right, p, best, closest, bestTri, bestFeat);
            if (dL < best) _queryNearest(n.left,  p, best, closest, bestTri, bestFeat);
        }
    }

    real_t _queryWinding(index_t nodeIdx, const Point & p, real_t beta) const
    {
        const Node & n = m_nodes[nodeIdx];
        if (n.left < 0) // leaf: always exact (near field matters most here)
        {
            real_t omega = 0;
            for (index_t i = n.triStart; i != n.triStart + n.triCount; ++i)
                omega += signedSolidAngle(p, m_orderedTris[3*i+0],
                                             m_orderedTris[3*i+1],
                                             m_orderedTris[3*i+2]);
            return omega;
        }

        const Point  r    = n.centroid - p;
        const real_t dist = r.norm();
        if (dist > beta * n.radius)
        {
            // Far-field dipole approximation: Omega ~ normalArea . (centroid-p) / dist^3
            // (sign/formula verified numerically against signedSolidAngle()).
            return n.normalArea.dot(r) / (dist*dist*dist);
        }

        return _queryWinding(n.left, p, beta) + _queryWinding(n.right, p, beta);
    }

    std::vector<Node>       m_nodes;
    std::vector<Point>      m_orderedTris; // 3 per triangle, leaves index contiguous ranges
    std::vector<TriNormals> m_orderedNrm;  // 1 per triangle, same ordering as above
    index_t                 m_root;
};

} // namespace gismo
