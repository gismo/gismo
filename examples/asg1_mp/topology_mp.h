/** @file topology_mp.h

    @brief Header-only multi-patch topology utilities for the AS-G1
           multi-patch basis construction (Phase A).

    Enumerates the three building blocks of the Kapl-Sangalli-Takacs
    decomposition  Global = Interior + Edge + Vertex  for an arbitrary
    planar multi-patch geometry:

      * `EdgeInfo`   -- one per topological edge (interface or boundary)
      * `VertexInfo` -- one per topological vertex (collection of the
                        patch-corners that meet at one physical point)

    Public API (in namespace `gismo::asg1mp`):

      * `enumerateEdges(mp)`     -> std::vector<EdgeInfo>
      * `enumerateVertices(mp)`  -> std::vector<VertexInfo>
      * helper `tangentDirection(side)`, `cornerOfSideEnd(...)`

    The enumeration mirrors the pattern used by
    gsUnstructuredSplines::gsApproxC1Spline (which walks
    `mp.vertices()` and inspects each corner's containing sides), but
    keeps everything in light POD structs so the rest of the
    construction is topology-driven and does not depend on any external
    optional module.

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include <vector>
#include <gismo.h>

namespace gismo {
namespace asg1mp {

// ====================================================================
// Edge (interface or boundary) descriptor
// ====================================================================

/// One topological edge of the multi-patch domain.
///
/// For an *interior* edge (`isBoundary == false`) the edge is shared
/// by two patches and `ifc` holds the `boundaryInterface`.  For a
/// *boundary* edge (`isBoundary == true`) only `side1` is meaningful.
template <class T>
struct EdgeInfo
{
    index_t            id        = -1;     ///< index in enumeration
    bool               isBoundary = false; ///< true: domain boundary edge
    boundaryInterface  ifc;                ///< valid iff !isBoundary

    patchSide          side1;              ///< first  patch side
    patchSide          side2;              ///< second patch side (interior only)

    /// Tangential directions run opposite across the interface
    /// (`!ifc.dirOrientation(side1, tangentDir(side1))`).
    bool               flipped   = false;

    /// 8 gluing-data numbers in patch-1's tangential parametrisation:
    ///   [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1].
    /// For boundary edges this stays (1,1,0,0,1,1,0,0).
    gsVector<T>        gluingData;

    bool               gluingOk  = false;  ///< gluing data solve succeeded

    EdgeInfo()
    : side1(0, boundary::west), side2(0, boundary::west)
    {
        gluingData.setZero(8);
        gluingData(0) = gluingData(1) = T(1);     // alpha_1 = 1
        gluingData(4) = gluingData(5) = T(1);     // alpha_2 = 1
    }
};

// ====================================================================
// Vertex descriptor
// ====================================================================

/// One topological vertex of the multi-patch domain: the set of
/// patch-corners that are glued to a single physical point.
template <class T>
struct VertexInfo
{
    index_t                  id        = -1;   ///< index in enumeration
    gsVector<T>              physPoint;         ///< x^(i) in R^2
    bool                     isBoundary = false;///< lies on domain boundary
    index_t                  valence   = 0;     ///< number of incident patches

    /// All patch-corners glued at this vertex (as returned by
    /// `gsBoxTopology::getCornerList`).
    std::vector<patchCorner> corners;

    /// For every corner, the two patch-sides of that patch meeting at
    /// the corner (lexicographic side order from `getContainingSides`).
    std::vector<std::pair<patchSide, patchSide> > cornerSides;

    /// For every corner, whether each of the two containing sides is an
    /// interface (true) or a domain boundary (false).
    std::vector<std::pair<bool, bool> >           cornerSideIsInterface;
};

// ====================================================================
// Small helpers
// ====================================================================

/// Tangential parametric direction of a 2-D box side
/// (the direction that varies along the side).
inline short_t tangentDirection(const boxSide& s)
{
    return 1 - s.direction();
}

// ====================================================================
// Edge enumeration
// ====================================================================

/// Enumerate every interface and boundary edge of `mp`.
/// `mp.computeTopology()` is assumed to have been called already.
template <class T>
std::vector<EdgeInfo<T> > enumerateEdges(const gsMultiPatch<T>& mp)
{
    std::vector<EdgeInfo<T> > edges;
    edges.reserve(mp.nInterfaces() + mp.nBoundary());

    // Interior edges (interfaces)
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        EdgeInfo<T> e;
        e.id         = static_cast<index_t>(edges.size());
        e.isBoundary = false;
        e.ifc        = *it;
        e.side1      = it->first();
        e.side2      = it->second();
        const short_t tdir1 = tangentDirection(e.side1);
        e.flipped    = !it->dirOrientation(e.side1, tdir1);
        edges.push_back(e);
    }

    // Boundary edges
    const std::vector<patchSide>& bnd = mp.boundaries();
    for (size_t i = 0; i < bnd.size(); ++i)
    {
        EdgeInfo<T> e;
        e.id         = static_cast<index_t>(edges.size());
        e.isBoundary = true;
        e.side1      = bnd[i];
        edges.push_back(e);
    }
    return edges;
}

// ====================================================================
// Vertex enumeration
// ====================================================================

/// Enumerate every topological vertex of `mp`.
/// `mp.computeTopology()` is assumed to have been called already.
template <class T>
std::vector<VertexInfo<T> > enumerateVertices(const gsMultiPatch<T>& mp)
{
    const short_t dim = 2;
    std::vector<VertexInfo<T> > result;

    // `vertices()` returns, for each topological vertex, the list of all
    // equivalent patch-corners (already equivalenced via getCornerList).
    std::vector<std::vector<patchCorner> > allVerts = mp.vertices();
    result.reserve(allVerts.size());

    for (size_t v = 0; v < allVerts.size(); ++v)
    {
        VertexInfo<T> vi;
        vi.id      = static_cast<index_t>(v);
        vi.corners = allVerts[v];
        vi.valence = static_cast<index_t>(vi.corners.size());

        // Physical point from the first corner.
        const patchCorner& pc0 = vi.corners.front();
        gsVector<bool> parB = pc0.parameters(dim);   // dim entries in {0,1}
        gsMatrix<T> par(dim, 1);
        for (short_t d = 0; d < dim; ++d)
            par(d, 0) = parB(d) ? T(1) : T(0);
        gsMatrix<T> phys;
        mp.patch(pc0.patch).eval_into(par, phys);
        vi.physPoint = phys.col(0);

        // For each corner, record its two containing sides and whether
        // each is an interface.  A vertex is a boundary vertex iff any
        // incident side is a domain boundary.
        bool boundaryVertex = false;
        vi.cornerSides.reserve(vi.corners.size());
        vi.cornerSideIsInterface.reserve(vi.corners.size());
        for (size_t k = 0; k < vi.corners.size(); ++k)
        {
            std::vector<patchSide> sides;
            vi.corners[k].getContainingSides(dim, sides);
            GISMO_ASSERT(sides.size() == 2,
                         "A 2-D corner must have exactly two sides.");
            vi.cornerSides.push_back(std::make_pair(sides[0], sides[1]));

            const bool isI0 = mp.isInterface(sides[0]);
            const bool isI1 = mp.isInterface(sides[1]);
            vi.cornerSideIsInterface.push_back(std::make_pair(isI0, isI1));
            if (!isI0 || !isI1) boundaryVertex = true;
        }
        vi.isBoundary = boundaryVertex;

        result.push_back(vi);
    }
    return result;
}

} // namespace asg1mp
} // namespace gismo
