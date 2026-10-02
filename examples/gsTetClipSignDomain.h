/** @file gsTetClipSignDomain.h

    @brief A gsTrimmedDomain<3,real_t> whose kd-tree leaves are exactly the
    level-0 cells of a uniform gsTetClip::Grid3, with each leaf's sign
    copied -- never resampled -- from an already-computed
    gsTetClip::CellStatus classification (gsTetMeshClip.h /
    gsImmersedLookupRule.h). Also provides the background-space helpers
    (backgroundBasis, identityBoxGeometry) that this class and the
    driver's FE space / geometry map share.

    Theory / design: gsTrimmedDomain<d,T>::init(gsTensorBSplineBasis<d,T>,
    samples) drives ITS OWN classification by sampling this->sign() at
    Lobatto points inside every leaf -- appropriate when "in/out" is a
    genuine function of the point, evaluated for the first time. Here it is
    not: the exact clip status of every cell was already computed once by
    the tet-mesh clipper, and sampling it again would at best waste that
    work and at worst (a status decided by a KahanSum'd clipped volume
    against a tolerance, gsTetMeshClip.h) disagree with the clipper's own
    verdict through a completely different, coarser test. sign() therefore
    does no classification of its own: it locates the ONE background cell a
    query's centroid falls in (colOf, gsTetMeshClip.h) and returns that
    cell's already-known status, converted once to the trimmed-domain sign
    convention. Its only substantive job is enforcing -- via a single-cell
    GISMO_ENSURE -- the one invariant this reuse depends on: every call must
    carry points from ONE cell, never a query that would need two different
    stored answers averaged into one.

    Sign convention (gsTrimmedDomain.h:141-150 vs gsTetMeshClip.h's
    CellStatus): the two conventions are OPPOSITE on Full/Empty and agree on
    Cut. Converted once, by an explicit switch, in the constructor -- not by
    negation, so that an out-of-range status value is caught by the switch's
    default branch instead of silently producing a plausible-looking wrong
    sign.

    | gsTetClip::CellStatus | value | trimmed-domain sign | selected by    |
    |------------------------|-------|----------------------|----------------|
    | Full                   |  +1   | -1 (interior)         | InteriorSign   |
    | Cut                    |   0   |  0 (cut)              | BoundarySign   |
    | Empty                  |  -1   | +1 (exterior)         | ExteriorSign   |

    Audit of every sign()/inDomain()/onBoundary() call site in src/ matching
    `.sign(`, `->sign(`, `::sign(`, unqualified `sign(`, `inDomain(`,
    `onBoundary(`, excluding the look-alikes gsCutTreeData::sign(), the
    domain ITERATOR's own sign(), Eigen's .array().sign() and the
    GhostFace/SkeletonFace::sign(size_t) functors):

      | call site                                   | points passed          | multi-leaf possible? |
      |----------------------------------------------|-------------------------|------------------------|
      | gsTrimmedDomain.h:563 _classifyLeaf           | Lobatto samples of ONE leaf (:561-562) | No -- its only non-forbidden caller is _classifyTree's Phase B (:620), which runs AFTER Phase A has already split every leaf to a single level-0 cell (:578-607); no sign() call happens before that split |
      | gsTrimmedDomain.h:266,276 inDomain/onBoundary | arbitrary               | Would be pointwise by contract, but zero callers exist anywhere in src/ |
      | ghost/skeleton faces (:316-363 -> _elementSignGrid :445-504) | -- | Read only the CACHED leaf sign, never call sign() |
      | dof elimination (driver code, e.g. poisson2_ghost_penalty_example.cpp:268-293) | -- | Walks begin<InteriorSign>/begin<BoundarySign>; cached signs again |
      | gsAssembler/gsQuadrature.h cut-cell/Algoim factories | -- | Never call sign(); dynamic_cast the domain to gsImplicitTrimmedDomain<1/2/3> and warn-once-then-fall-back-to-plain-Gauss on any other domain type (including this one) if no custom quadrature factory is installed |

    Conclusion: in src/, sign() is only ever fed the samples of a single
    leaf, so the per-CELL (not pointwise) semantics implemented below are
    safe for every existing call path. inDomain()/onBoundary() on this class
    therefore return a per-cell constant, not a pointwise verdict -- callers
    that need a true pointwise test must not use them on this class.

    Only the tensor-basis init() overload is valid here. The two ADAPTIVE
    overloads (HTB, size-based) are forbidden outright: they drive
    _classifyTreeAdaptive(), whose split decision consumes the sign itself,
    so a leaf CAN be queried while still spanning several level-0 cells --
    exactly the multi-cell query the centroid-lookup trick cannot answer
    correctly. The bbox+numCells overload is not adaptive (it also runs the
    single-cell-before-any-sign()-call _classifyTree Phase A,
    gsTrimmedDomain.h:825), so it is not forbidden for that reason; it is
    simply not used, because its breaks come from the UNIFORM
    `gsKnotVector(lo,hi,numCells-1,2)` constructor (gsTrimmedDomain.h:
    816-817), recomputing every interior break instead of reproducing the
    grid lines bitwise -- the same reason gridKnots() below exists rather
    than reusing that constructor. All four init() overloads are protected
    (gsTrimmedDomain.h:429-979), so the constructor below -- which calls
    Base::init(basis, samples) and nothing else -- is the only entry point;
    outside code cannot call init() on this class at all.

    Complexity: construction makes n^3 sign() calls (one per leaf; Phase A
    subdivision is a constant number of tree operations per leaf and costs
    no level-set evaluation), each call passed the samples^3 Lobatto points
    of its own leaf and costing O(n + samples^3): colOf is an O(n) linear
    scan of the grid lines, called three times per call (once per axis) on
    the CENTROID only -- not once per one of the up to samples^3 points in
    \a u -- and the constant array-fill of the returned sign vector is
    O(samples^3). Overall O(n^3 * (n + samples^3)) for an n^3 grid. A single
    sign() call thereafter, e.g. from the guard-test checks, is likewise
    O(m) for m input points (three colOf scans plus an O(m) fill).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "gsTetMeshClip.h"
#include <gsDomain/gsTrimmedDomain.h>

#include <limits>
#include <vector>

namespace gsTetClip
{
using namespace gismo;

/// Builds the knot container {L[0] x(p+1), L[1..n-1], L[n] x(p+1)} bitwise
/// from the grid-line array \a L (size n+1, e.g. gsTetClip::gridLines'
/// output) and wraps it in an open gsKnotVector of degree \a p, via the
/// knotContainer constructor `gsKnotVector(knotContainer, degree)`
/// (gsKnotVector.h:208), which stores the values verbatim.
///
/// Deliberately NOT the uniform constructor
/// `gsKnotVector(first,last,interior,mult_ends)` (gsKnotVector.h:505-510):
/// that RECOMPUTES every interior knot as
/// `first + i*((last-first)/(interior+1))` (gsKnotVector.hpp:675,686-688),
/// which equals L[i] = x0 + i*h bitwise only when (last-first)/n happens to
/// round to exactly h -- true on the checks' own box (h = 2/n, a power of
/// two) but not in general. gsTetClip's own cell lookups (colOf,
/// TetClipSignDomain::sign() below) compare against L itself, so the basis'
/// breaks must literally BE L, not a recomputation that merely agrees with
/// it to rounding.
inline gsKnotVector<real_t> gridKnots(const std::vector<real_t> & L, short_t p)
{
    GISMO_ASSERT(L.size() >= 2, "gridKnots: need at least 2 grid lines, got " << L.size());
    std::vector<real_t> kc;
    kc.reserve(L.size() + 2*(size_t)p);
    for (short_t m = 0; m <= p; ++m) kc.push_back(L.front());
    for (size_t i = 1; i + 1 < L.size(); ++i) kc.push_back(L[i]);
    for (short_t m = 0; m <= p; ++m) kc.push_back(L.back());
    return gsKnotVector<real_t>(kc, p);
}

/// The degree-\a p tensor B-spline basis whose breakpoints are bitwise the
/// grid lines of \a g (X, Y, Z from gridLines()): n elements per direction,
/// C^{p-1} across the interior lines, `knots(j).numElements()==g.n`. Used
/// both as the discretisation basis (via TetClipSignDomain's constructor)
/// and, wrapped in a single-basis gsMultiBasis, as the driver's FE space.
inline gsTensorBSplineBasis<3,real_t> backgroundBasis(const Grid3 & g, short_t p)
{
    std::vector<real_t> X, Y, Z;
    gridLines(g, X, Y, Z);
    return gsTensorBSplineBasis<3,real_t>(gridKnots(X,p), gridKnots(Y,p), gridKnots(Z,p));
}

/// An exactly affine identity map on the box of \a g: a degree-1 basis with
/// NO interior knots per direction (one element covering the whole box,
/// independent of g.n) and control points equal to its own anchors, so
/// physical == parametric everywhere on the box, not merely at grid lines --
/// a degree-1 open-knot B-spline with control points equal to the domain
/// endpoints reproduces the identity function exactly at any parametric
/// value via its linear-precision property, not only at knots/elements.
/// This is deliberately NOT tied to any discretisation degree p, and its
/// own element count need not match backgroundBasis(g,p)'s: the assembler
/// evaluates this geometry map at quadrature points supplied by the FE
/// space's OWN elements, and spline evaluation does not require the two
/// bases to share element boundaries, only the same parameter range (both
/// helpers build that range from the same Grid3 g).
///
/// gsImmersedLookupRule.h's lookup rules (ClipStreamer's own volRule()/
/// bdrRule() and everything built on it) key and produce PHYSICAL
/// coordinates directly; those only coincide with the parametric points the
/// assembler's quadrature-and-geometry pipeline expects when the driver's
/// OWN background geometry map -- built from this function -- is the
/// identity, as it is here. Since J = I everywhere on the box, a physical
/// normal equals the corresponding parametric one, so a boundary rule's
/// stored physical normals (gsImmersedLookupRule.h's BdrNormalField) need
/// no further transformation. Deliberately NOT
/// gsNurbsCreator<>::BSplineCube(r,x,y,z): its knots are on [0,1]^3
/// (gsNurbsCreator.hpp:759-774), so e.g. BSplineCube(2,-1,-1,-1) maps
/// [0,1]^3 -> [-1,1]^3 with Jacobian 2*I, not the identity.
inline gsMultiPatch<real_t> identityBoxGeometry(const Grid3 & g)
{
    std::vector<real_t> X, Y, Z;
    gridLines(g, X, Y, Z);
    const std::vector<real_t> Lx = { X.front(), X.back() };
    const std::vector<real_t> Ly = { Y.front(), Y.back() };
    const std::vector<real_t> Lz = { Z.front(), Z.back() };
    gsTensorBSplineBasis<3,real_t> gb(gridKnots(Lx,1), gridKnots(Ly,1), gridKnots(Lz,1));
    gsMatrix<real_t> coefs = gb.anchors().transpose(); // anchors(): d x N; makeGeometry wants N x d.
    return gsMultiPatch<real_t>(*gb.makeGeometry(give(coefs)));
}

/// gsTrimmedDomain<3,real_t> whose kd-tree has exactly one leaf per level-0
/// cell of a uniform gsTetClip::Grid3, with leaf signs copied from an
/// already-computed gsTetClip::CellStatus classification rather than
/// resampled (see the file header for the full design rationale, the sign
/// table, the src/ call-site audit and the complexity bound).
class TetClipSignDomain : public gismo::gsTrimmedDomain<3,real_t>
{
    typedef gismo::gsTrimmedDomain<3,real_t> Base;

public:
    /// \param grid    the uniform background grid; every one of its
    ///                 grid.n^3 cells becomes exactly one kd-tree leaf.
    /// \param status   one gsTetClip::CellStatus value per cell id, size
    ///                 grid.n^3, Grid3's own id = i + n*(j + n*k) convention.
    /// \param deg      the discretisation degree p of the space that will be
    ///                 assembled on this domain. Sets degree() via
    ///                 tbasis.maxDegree() inside Base::init
    ///                 (gsTrimmedDomain.h:832), which gsQuadrature uses to
    ///                 size rules; it must equal the degree of the space
    ///                 actually assembled here (gsTrimmedDomain.h:200-223).
    /// \param samples  Lobatto samples per direction Base::init uses to
    ///                 probe each single-cell leaf. The sample COUNT never
    ///                 changes the classification result here (every sample
    ///                 of one cell returns the same stored constant, see
    ///                 sign() below); >= 2 is enforced defensively because
    ///                 whether gsLobattoRule accepts fewer was not verified.
    TetClipSignDomain(const Grid3 & grid, const std::vector<int> & status,
                       short_t deg, index_t samples = 5)
    : m_grid(grid)
    {
        GISMO_ENSURE(grid.n >= 1, "TetClipSignDomain: grid.n must be >= 1, got " << grid.n);
        GISMO_ENSURE(grid.h > 0, "TetClipSignDomain: grid.h must be positive, got " << grid.h);
        const size_t n3 = (size_t)grid.n*(size_t)grid.n*(size_t)grid.n;
        GISMO_ENSURE(status.size() == n3, "TetClipSignDomain: status.size()=" << status.size()
                    << " does not match grid.n^3=" << n3);
        GISMO_ENSURE(samples >= 2, "TetClipSignDomain: samples must be >= 2, got " << samples);

        gridLines(grid, m_X, m_Y, m_Z);

        // All member state sign() reads must be filled BEFORE Base::init:
        // init() calls this->sign() during its own classification pass.
        m_sign.resize(n3);
        for (size_t id = 0; id != n3; ++id)
        {
            switch (status[id])
            {
            case Full:  m_sign[id] = -1; break;
            case Cut:   m_sign[id] =  0; break;
            case Empty: m_sign[id] =  1; break;
            default:
                GISMO_ENSURE(false, "TetClipSignDomain: status[" << id << "]=" << status[id]
                            << " is not a gsTetClip::CellStatus value.");
            }
        }

        m_tol = 1e-10*grid.h + 64*std::numeric_limits<real_t>::epsilon()*boxCornerAbsMax(grid);

        // The basis is a temporary: Base::init copies its breaks
        // (gsTrimmedDomain.h:836-837) and the domain holds no basis
        // reference afterwards (gsTrimmedDomain.h:44-49).
        Base::init(backgroundBasis(grid, deg), samples);

        GISMO_ENSURE(1 == this->numLevels(), "TetClipSignDomain: Base::init produced "
                    << this->numLevels() << " kd-tree levels; this class supports a "
                    "single-level background grid only.");
        GISMO_ENSURE(n3 == this->numElements<AnySign>(), "TetClipSignDomain: Base::init "
                    "produced " << this->numElements<AnySign>() << " elements, expected "
                    "grid.n^3=" << n3 << ".");
    }

    /// Returns the CELL status (converted to the trimmed-domain sign, see
    /// the file header table) for every column of \a u, as a constant --
    /// not a pointwise evaluation of anything. The cell is located from the
    /// column centroid via colOf() (gsTetMeshClip.h), and every column's
    /// bounding box is enforced (within \a m_tol) to lie in that ONE cell:
    /// this is the invariant every call in src/ already satisfies (see the
    /// file header's audit table), and the only overload of
    /// gsTrimmedDomain::init() this class permits (the tensor-basis one)
    /// guarantees it structurally, since its Phase A splits every leaf to a
    /// single cell before any sign() call. Thread-safe: reads only members
    /// that are immutable after construction (grid lines, the sign table,
    /// m_tol), no caches, no mutable state -- required because
    /// gsTrimmedDomain classifies leaves under `#pragma omp parallel`.
    /// inDomain()/onBoundary() (gsTrimmedDomain's non-virtual wrappers
    /// around sign()) are therefore per-cell for this class too, not
    /// pointwise.
    gsVector<short_t> sign(const gsMatrix<real_t> & u) override
    {
        GISMO_ENSURE(3 == u.rows() && u.cols() > 0, "TetClipSignDomain::sign: expects 3D "
                    "points, got u.rows()=" << u.rows() << ", u.cols()=" << u.cols());
        const gsVector<real_t> lo = u.rowwise().minCoeff();
        const gsVector<real_t> hi = u.rowwise().maxCoeff();
        const gsVector<real_t> c  = u.rowwise().mean();
        const index_t i = colOf(c[0], m_X), j = colOf(c[1], m_Y), k = colOf(c[2], m_Z);
        GISMO_ENSURE(lo[0] >= m_X[i]-m_tol && hi[0] <= m_X[i+1]+m_tol &&
                     lo[1] >= m_Y[j]-m_tol && hi[1] <= m_Y[j+1]+m_tol &&
                     lo[2] >= m_Z[k]-m_tol && hi[2] <= m_Z[k+1]+m_tol,
                     "TetClipSignDomain::sign: the " << u.cols() << " points span more than one "
                     "background cell (centroid cell " << i << "," << j << "," << k << "); this "
                     "domain is valid only under the single-cell init(tbasis, samples).");
        const index_t n = m_grid.n;
        gsVector<short_t> s(u.cols());
        s.setConstant(m_sign[(size_t)i + (size_t)n*((size_t)j + (size_t)n*(size_t)k)]);
        return s;
    }

    /// 3x2, row j = direction j, col 0/1 = lower/upper corner (the
    /// gsTrimmedDomain convention, :743-745), taken from the grid extremes.
    gsMatrix<real_t> boundingBox() const override
    {
        gsMatrix<real_t> bb(3,2);
        bb(0,0) = m_X.front(); bb(0,1) = m_X.back();
        bb(1,0) = m_Y.front(); bb(1,1) = m_Y.back();
        bb(2,0) = m_Z.front(); bb(2,1) = m_Z.back();
        return bb;
    }

    const Grid3 & grid() const { return m_grid; }
    short_t cellSign(index_t id) const { return m_sign[id]; }

private:
    Grid3 m_grid;
    std::vector<real_t> m_X, m_Y, m_Z;
    std::vector<short_t> m_sign;
    real_t m_tol;
};

} // namespace gsTetClip
