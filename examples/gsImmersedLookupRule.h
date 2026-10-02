/** @file gsImmersedLookupRule.h

    @brief Streamed, cell-major tet-clip quadrature and the lookup
    quadrature-rule/normal-field/RAII-scope machinery built on top of it.

    gsTetMeshClip.h's tetClipQuadrature() computes and STORES a CellRule3 /
    CellBdrRule3 for every cell of the background grid at once (a Grid3 of
    n^3 cells holds two full arrays of that size, each cell's nodes/weights
    materialized whether or not the assembler ever visits that cell). This
    file instead streams the SAME clip on demand, one cell at a time, from a
    compact per-cell bucket of overlapping (tet, cell) / (triangle, cell)
    pairs (ClipStreamer) -- its volRule()/bdrRule() output is bitwise
    identical to tetClipQuadrature()'s per-cell CellRule3/CellBdrRule3 (see
    the driver's `--study check`, which verifies this directly), but no
    per-cell rule is ever held in memory beyond the call that produced it.

    On top of the streamer sit:
    - VolLookupRule / BdrLookupRule, gsQuadRule<real_t> subclasses keyed on
      the element-box midpoint, dispatching to a tensor-Gauss rule (Full
      cells), an empty rule (Empty cells) or a VolCellSource/BdrCellSource
      (Cut cells / every cell respectively). The source is either a
      ClipStreamer (an on-the-fly clip rule) or a VolCellTable /
      BdrCellTable (a table-backed rule, holding precomputed,
      e.g. compressed, per-cell rules) -- both classes serve either source unchanged;
    - BdrNormalField, a gsFunction<real_t> returning the stored physical
      unit normals for exactly the nodes a BdrLookupRule produced, keyed by
      centroid-plus-tolerance candidate search and a bitwise node match;
    - QuadratureScope and its VolumeQuadratureScope/BoundaryQuadratureScope
      specializations, RAII scopes that install a quadrature factory on a
      gsExprAssembler or gsExprEvaluator and restore both the factory state
      and the "quDim" option on exit.

    All physical quantities (nodes, weights, normals) are in the SAME
    physical coordinates as gsTetMeshClip.h's TetMesh/Grid3; the background
    geometry map used with these rules must be the identity (see
    makeBackground() in immersed_tetmesh_poisson_example.cpp).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "gsTetMeshClip.h"

#include <cstring>
#include <limits>
#include <vector>

namespace gsTetClip
{
using namespace gismo;

//----------------------------------------------------------------------------
// CellIndex: the grid, its precomputed lines, and the per-cell status,
// shared read-only by every lookup rule / source / field below.
//----------------------------------------------------------------------------

/// Shared, read-only description of the background Cartesian grid that
/// keys every lookup rule in this file: the Grid3 itself, its precomputed
/// grid-line positions X/Y/Z (the SAME arrays every colOf() call here
/// compares against -- the ULP trap documented at gsTetMeshClip.h:287-295
/// applies verbatim), and the per-cell CellStatus (Empty/Cut/Full), filled
/// once by ClipStreamer's status pass or copied in by the caller from a
/// TetClipQuadrature. Units: physical coordinates. Thread-safety: safe to
/// read concurrently once built; never mutated after construction
/// completes (ClipStreamer's constructor is the only writer).
struct CellIndex
{
    Grid3 grid;
    std::vector<real_t> X, Y, Z;
    std::vector<int> status; // CellStatus per cell id, size n^3

    /// Flat cell id, Grid3's own convention: i + n*(j + n*k).
    size_t id(index_t i, index_t j, index_t k) const
    { return (size_t)i + (size_t)grid.n*((size_t)j + (size_t)grid.n*(size_t)k); }

    /// The cell containing physical point (x,y,z), by the SAME colOf()
    /// lookup tetClipQuadrature uses against these SAME grid lines.
    size_t cellOfPoint(real_t x, real_t y, real_t z) const
    { return id(colOf(x,X), colOf(y,Y), colOf(z,Z)); }

    /// Inverse of id(): recovers (i,j,k) from a flat cell id.
    void ijk(size_t cellId, index_t & i, index_t & j, index_t & k) const
    {
        const size_t n = (size_t)grid.n;
        i = (index_t)(cellId % n);
        j = (index_t)((cellId / n) % n);
        k = (index_t)(cellId / (n*n));
    }
};

/// Buckets every tet/boundary-triangle index of \a M into the background
/// cells its bounding box overlaps, using the SAME bbox and colOf() ranges
/// as tetClipQuadrature's own volume/boundary loops (gsTetMeshClip.h, the
/// blocks preceding its `k=k0..k1,j=j0..j1,i=i0..i1` loops), so that a
/// ClipStreamer built on these buckets visits exactly the same (tet, cell)
/// / (triangle, cell) pairs tetClipQuadrature does. \a tets and \a tris are
/// resized to n^3 and assigned fresh (any prior content is discarded).
/// Because the outer loop over t (resp. bt) runs in increasing order and a
/// bucket is only ever appended to, every bucket ends up sorted in
/// increasing tet/triangle id -- the same ORDER tetClipQuadrature's own
/// loops visit them in, which the streamed rules' bitwise equality to
/// tetClipQuadrature depends on.
///
/// Complexity: O((#tets + #tris) * (n + overlapped cells)), n = cells per
/// direction (colOf is an O(n) linear scan of the grid lines).
inline void buildCellBuckets(const TetMesh & M, const CellIndex & idx,
                             std::vector<std::vector<index_t> > & tets,
                             std::vector<std::vector<index_t> > & tris)
{
    const index_t n = idx.grid.n;
    const size_t nCells = (size_t)n*(size_t)n*(size_t)n;
    tets.assign(nCells, std::vector<index_t>());
    tris.assign(nCells, std::vector<index_t>());

    for (size_t t = 0; t != M.tet.size(); ++t)
    {
        const std::array<index_t,4> & tt = M.tet[t];
        real_t xmin=M.P[tt[0]][0], xmax=xmin, ymin=M.P[tt[0]][1], ymax=ymin, zmin=M.P[tt[0]][2], zmax=zmin;
        for (int v = 1; v != 4; ++v)
        {
            const Vec3 & p = M.P[tt[v]];
            xmin = math::min(xmin,p[0]); xmax = math::max(xmax,p[0]);
            ymin = math::min(ymin,p[1]); ymax = math::max(ymax,p[1]);
            zmin = math::min(zmin,p[2]); zmax = math::max(zmax,p[2]);
        }
        const index_t i0 = colOf(xmin,idx.X), i1 = colOf(xmax,idx.X);
        const index_t j0 = colOf(ymin,idx.Y), j1 = colOf(ymax,idx.Y);
        const index_t k0 = colOf(zmin,idx.Z), k1 = colOf(zmax,idx.Z);
        for (index_t k = k0; k <= k1; ++k)
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
            tets[idx.id(i,j,k)].push_back((index_t)t);
    }

    for (size_t bt = 0; bt != M.bdrTri.size(); ++bt)
    {
        const std::array<index_t,4> & bTri = M.bdrTri[bt];
        const Vec3 & A = M.P[bTri[0]], & B = M.P[bTri[1]], & C = M.P[bTri[2]];
        const real_t xmin = math::min(A[0], math::min(B[0],C[0])), xmax = math::max(A[0], math::max(B[0],C[0]));
        const real_t ymin = math::min(A[1], math::min(B[1],C[1])), ymax = math::max(A[1], math::max(B[1],C[1]));
        const real_t zmin = math::min(A[2], math::min(B[2],C[2])), zmax = math::max(A[2], math::max(B[2],C[2]));
        const index_t i0 = colOf(xmin,idx.X), i1 = colOf(xmax,idx.X);
        const index_t j0 = colOf(ymin,idx.Y), j1 = colOf(ymax,idx.Y);
        const index_t k0 = colOf(zmin,idx.Z), k1 = colOf(zmax,idx.Z);
        for (index_t k = k0; k <= k1; ++k)
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
            tris[idx.id(i,j,k)].push_back((index_t)bt);
    }
}

//----------------------------------------------------------------------------
// Per-manifold sources of a cell's clip rule: either the on-the-fly
// streamer (ClipStreamer) or a stored table (VolCellTable/BdrCellTable,
// holding precomputed, e.g. compressed, per-cell rules).
//----------------------------------------------------------------------------

/// Source of a background cell's volume quadrature rule.
class VolCellSource
{
public:
    virtual ~VolCellSource() {}

    /// Volume rule of cell \a id: nodes 3 x m (physical coordinates),
    /// weights m (physical, already includes the sub-tet Jacobian). \a m
    /// may be 0 (an Empty cell, or a Cut cell whose only pieces were
    /// dropped by the exact-zero test). Must be safe to call concurrently.
    virtual void volRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const = 0;
};

/// Source of a background cell's boundary quadrature rule.
class BdrCellSource
{
public:
    virtual ~BdrCellSource() {}

    /// Boundary rule of cell \a id: nodes 3 x m, PHYSICAL surface weights
    /// m, PHYSICAL outward unit normals 3 x m (constant per owning
    /// boundary triangle, repeated per node). \a m may be 0. Must be safe
    /// to call concurrently.
    virtual void bdrRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights,
                         gsMatrix<real_t> & normals) const = 0;
};

/// Streamed, cell-major exact tet-clip quadrature: reproduces
/// tetClipQuadrature()'s per-cell CellRule3/CellBdrRule3 bitwise, one cell
/// at a time, from a compact bucket of overlapping (tet, cell) /
/// (triangle, cell) pairs -- it never materializes an n^3 array of
/// per-cell rules the way tetClipQuadrature() does. Only the buckets and
/// the per-cell status are held as members; volRule()/bdrRule() clip their
/// cell's bucket afresh on every call. Units, tensor shapes and the
/// degree-p collapsed-Gauss rules are exactly gsTetMeshClip.h's (see its
/// file header): volDeg = bdrDeg = 6p, 392 volume nodes per surviving
/// sub-tet and 49 boundary nodes per fan triangle at p=2.
///
/// Complexity: one volRule(id,...) call costs
/// O(Sum over cell id's tet bucket of (surviving pieces x nodes per
/// piece)); one bdrRule(id,...) call is the boundary analogue over the
/// triangle bucket. The constructor's status pass costs one FULL clip (it
/// must call volRule on every non-empty cell to classify it), i.e. the
/// same asymptotic cost as tetClipQuadrature()'s volume loop, run once
/// under `#pragma omp parallel for schedule(dynamic)` (each iteration
/// writes only its own status[id], so there is no race; the three cell
/// counts are tallied serially afterwards).
class ClipStreamer : public VolCellSource, public BdrCellSource
{
public:
    /// Builds the cell index, the five collapsed-Gauss Rule1Ds (exactly as
    /// tetClipQuadrature's own setup, gsTetMeshClip.h), the tet/triangle
    /// buckets, and runs the status pass. \a mesh must outlive this object
    /// (volRule()/bdrRule() dereference it on every call).
    ClipStreamer(memory::shared_ptr<const TetMesh> mesh, const Grid3 & grid, index_t p)
    : m_mesh(give(mesh)), m_index(memory::make_shared(new CellIndex())),
      m_numCut(0), m_numFull(0), m_numEmpty(0)
    {
        m_index->grid = grid;
        gridLines(grid, m_index->X, m_index->Y, m_index->Z);
        const size_t nCells = (size_t)grid.n*(size_t)grid.n*(size_t)grid.n;
        m_index->status.assign(nCells, (int)Empty);

        const index_t volDeg = volDegree(p);
        m_Ru = gauss01(ceilHalf(volDeg+3));
        m_Rv = gauss01(ceilHalf(volDeg+2));
        m_Rw = gauss01(ceilHalf(volDeg+1));
        const index_t bdrDeg = bdrDegree(p);
        m_Rbu = gauss01(ceilHalf(bdrDeg+2));
        m_Rbv = gauss01(ceilHalf(bdrDeg+1));

        buildCellBuckets(*m_mesh, *m_index, m_tetBuckets, m_triBuckets);

        // Status block, copied verbatim from tetClipQuadrature
        // (gsTetMeshClip.h): the classification tolerance depends only on
        // the grid, so these are computed once, outside the parallel loop.
        const real_t L = boxCornerAbsMax(grid);
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        const real_t hh = grid.h*grid.h;
        const real_t tolAbs = 16*eps*L*hh;
        const real_t box = grid.h*hh;

        const long nCellsL = (long)nCells;
        #pragma omp parallel for schedule(dynamic)
        for (long idL = 0; idL < nCellsL; ++idL)
        {
            const size_t id = (size_t)idL;
            if (m_tetBuckets[id].empty()) continue; // already Empty
            gsMatrix<real_t> nodes; gsVector<real_t> weights;
            volRule(id, nodes, weights);
            if (0 == weights.size()) continue; // already Empty
            CellRule3 tmp; tmp.weights = weights;
            const real_t vol = cellVolume(tmp);
            m_index->status[id] = (vol >= box-(1e-14*box+tolAbs) && vol <= box*(1.0+1e-14)+tolAbs)
                                 ? (int)Full : (int)Cut;
        }

        for (size_t id = 0; id != nCells; ++id)
        {
            if      (Full  == m_index->status[id]) ++m_numFull;
            else if (Cut   == m_index->status[id]) ++m_numCut;
            else                                   ++m_numEmpty;
        }
    }

    /// Volume loop body of tetClipQuadrature (gsTetMeshClip.h), restricted
    /// to cell \a id's own box and its own tet bucket, packed exactly as
    /// tetClipQuadrature packs a CellRule3. const, thread-safe.
    void volRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        index_t i,j,k; m_index->ijk(id, i,j,k);
        const real_t xL=m_index->X[i], xR=m_index->X[i+1];
        const real_t yB=m_index->Y[j], yT=m_index->Y[j+1];
        const real_t zN=m_index->Z[k], zF=m_index->Z[k+1];

        std::vector<real_t> vx,vy,vz,vw;
        for (index_t t : m_tetBuckets[id])
        {
            const std::array<index_t,4> & tt = m_mesh->tet[t];
            const Tet4 T = { m_mesh->P[tt[0]], m_mesh->P[tt[1]], m_mesh->P[tt[2]], m_mesh->P[tt[3]] };
            const std::vector<Tet4> pieces = clipTetToCell(T, xL,xR, yB,yT, zN,zF);
            for (const Tet4 & piece : pieces)
            {
                const real_t Dp = det3(sub3(piece[1],piece[0]), sub3(piece[2],piece[0]), sub3(piece[3],piece[0]));
                if (0.0 == Dp) continue;
                addCollapsedTetNodes(piece[0],piece[1],piece[2],piece[3], math::abs(Dp),
                                     m_Ru,m_Rv,m_Rw, vx,vy,vz,vw);
            }
        }
        nodes.resize(3, (index_t)vx.size());
        weights.resize((index_t)vx.size());
        for (size_t kk = 0; kk != vx.size(); ++kk)
        {
            nodes(0,(index_t)kk) = vx[kk]; nodes(1,(index_t)kk) = vy[kk]; nodes(2,(index_t)kk) = vz[kk];
            weights[(index_t)kk] = vw[kk];
        }
    }

    /// Boundary loop body of tetClipQuadrature, restricted to cell \a id's
    /// own box and its own triangle bucket, packed exactly as
    /// tetClipQuadrature packs a CellBdrRule3. const, thread-safe.
    void bdrRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights,
                gsMatrix<real_t> & normals) const override
    {
        index_t i,j,k; m_index->ijk(id, i,j,k);
        const real_t xL=m_index->X[i], xR=m_index->X[i+1];
        const real_t yB=m_index->Y[j], yT=m_index->Y[j+1];
        const real_t zN=m_index->Z[k], zF=m_index->Z[k+1];

        std::vector<real_t> bx,by,bz,bw,bnx,bny,bnz;
        for (index_t bt : m_triBuckets[id])
        {
            const std::array<index_t,4> & bTri = m_mesh->bdrTri[bt];
            const Vec3 & A = m_mesh->P[bTri[0]], & B = m_mesh->P[bTri[1]], & C = m_mesh->P[bTri[2]];
            const Vec3 crAB = cross3(sub3(B,A), sub3(C,A));
            const real_t triMag = math::sqrt(dot3(crAB,crAB));
            const Vec3 normal = scale3(1.0/triMag, crAB);

            const Poly3 poly = clipTriangleToBox3(A,B,C, xL,xR, yB,yT, zN,zF);
            if (poly.size() < 3) continue;
            for (size_t f = 1; f+1 < poly.size(); ++f)
            {
                const Vec3 & Q0 = poly[0], & Q1 = poly[f], & Q2 = poly[f+1];
                const Vec3 cr = cross3(sub3(Q1,Q0), sub3(Q2,Q0));
                if (0.0==cr[0] && 0.0==cr[1] && 0.0==cr[2]) continue;
                const real_t mag = math::sqrt(dot3(cr,cr));
                addCollapsedTriNodes(Q0,Q1,Q2, mag, m_Rbu,m_Rbv, normal, bx,by,bz,bw, bnx,bny,bnz);
            }
        }
        nodes.resize(3, (index_t)bx.size());
        weights.resize((index_t)bx.size());
        normals.resize(3, (index_t)bx.size());
        for (size_t kk = 0; kk != bx.size(); ++kk)
        {
            nodes(0,(index_t)kk) = bx[kk]; nodes(1,(index_t)kk) = by[kk]; nodes(2,(index_t)kk) = bz[kk];
            weights[(index_t)kk] = bw[kk];
            normals(0,(index_t)kk) = bnx[kk]; normals(1,(index_t)kk) = bny[kk]; normals(2,(index_t)kk) = bnz[kk];
        }
    }

    memory::shared_ptr<const CellIndex> index() const { return m_index; }
    const std::vector<std::vector<index_t> > & tetBuckets() const { return m_tetBuckets; }
    const std::vector<std::vector<index_t> > & triBuckets() const { return m_triBuckets; }
    index_t numCut()   const { return m_numCut; }
    index_t numFull()  const { return m_numFull; }
    index_t numEmpty() const { return m_numEmpty; }

private:
    memory::shared_ptr<const TetMesh> m_mesh;
    memory::shared_ptr<CellIndex> m_index;
    Rule1D m_Ru, m_Rv, m_Rw, m_Rbu, m_Rbv;
    std::vector<std::vector<index_t> > m_tetBuckets, m_triBuckets;
    index_t m_numCut, m_numFull, m_numEmpty;
};

/// Table-backed VolCellSource: volRule(id,...) copies a stored CellRule3
/// verbatim. \a cell has size n^3 (one entry per background cell; an
/// Empty or never-filled entry carries 0 columns). Holds precomputed
/// (e.g. compressed) per-cell quadrature. The driver's `--study check`
/// fills it by moving tetClipQuadrature()'s Cut-cell rules in, and
/// verifies that a VolLookupRule over the table is bitwise identical to
/// one over a ClipStreamer.
struct VolCellTable : public VolCellSource
{
    std::vector<CellRule3> cell;
    void volRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    { nodes = cell[id].nodes; weights = cell[id].weights; }
};

/// Table-backed BdrCellSource, analogous to VolCellTable.
struct BdrCellTable : public BdrCellSource
{
    std::vector<CellBdrRule3> cell;
    void bdrRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights,
                gsMatrix<real_t> & normals) const override
    { nodes = cell[id].nodes; weights = cell[id].weights; normals = cell[id].normals; }
};

//----------------------------------------------------------------------------
// Lookup quadrature rules, keyed on the element-box midpoint.
//----------------------------------------------------------------------------

/// gsQuadRule<real_t> for the background grid's volume integral, keyed on
/// the element-box midpoint. Dispatches to a tensor-Gauss rule of FIXED
/// order (p+1,p+1,p+1) on Full cells, an empty rule on Empty cells, and
/// \a src's volRule() on Cut cells. \a src is either a ClipStreamer (an
/// on-the-fly clip rule) or a VolCellTable (a table-backed rule); this
/// ONE class serves both, since VolCellSource hides the difference.
///
/// Units: physical coordinates (the background geometry map must be the
/// identity, see makeBackground() in the driver). Tensor shapes: nodes
/// 3 x m, weights m. Thread-safety: mapTo() is const and touches only
/// read-only members, so concurrent calls (as from
/// gsExprAssembler::assemble's `#pragma omp parallel`) are safe as long as
/// \a src's own volRule() is (ClipStreamer's and VolCellTable's are).
///
/// Keying: the element box's midpoint lies h/2 from every grid line, so
/// ULP-level disagreement between the background basis's knots and
/// CellIndex's X/Y/Z cannot move it into a neighbouring cell; the
/// GISMO_ENSURE that |lower-X[i]| and |upper-X[i+1]| (per axis) are
/// <= 1e-8*h rejects a background basis whose elements are not the
/// Grid3 cells. BdrLookupRule uses the same keying.
///
/// Full-cell Gauss order (decision): \a nGaussFull is PINNED by the
/// caller as (p+1,p+1,p+1), not derived from the QuadratureFactory's
/// \a degrees argument (the factory ignores it), so the assembler and
/// evaluator paths use IDENTICAL Full
/// rules by construction. p+1 Gauss points per direction are exact on
/// Q_{2p+1} superset Q_{2p}, which covers the stiffness/mass integrands of
/// a degree-p background tensor space on the identity map; this matches
/// the assembler's own default quA=1, quB=1 (gsExprAssembler.h:1049-1050).
class VolLookupRule : public gsQuadRule<real_t>
{
public:
    VolLookupRule(memory::shared_ptr<const CellIndex> idx,
                  memory::shared_ptr<const VolCellSource> src,
                  const gsVector<index_t> & nGaussFull)
    : m_idx(give(idx)), m_src(give(src)), m_full(nGaussFull)
    { }

    using gsQuadRule<real_t>::mapTo;

    void mapTo(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
              gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        GISMO_ENSURE(3 == lower.size(), "VolLookupRule::mapTo: expects a 3D element box, got "
                    "lower.size()=" << lower.size());
        const gsVector<real_t> mid = 0.5*(lower+upper);
        const size_t id = m_idx->cellOfPoint(mid[0], mid[1], mid[2]);

        index_t i,j,k; m_idx->ijk(id, i,j,k);
        const real_t h = m_idx->grid.h;
        GISMO_ENSURE(math::abs(lower[0]-m_idx->X[i])   <= 1e-8*h && math::abs(upper[0]-m_idx->X[i+1]) <= 1e-8*h &&
                    math::abs(lower[1]-m_idx->Y[j])   <= 1e-8*h && math::abs(upper[1]-m_idx->Y[j+1]) <= 1e-8*h &&
                    math::abs(lower[2]-m_idx->Z[k])   <= 1e-8*h && math::abs(upper[2]-m_idx->Z[k+1]) <= 1e-8*h,
                    "VolLookupRule::mapTo: the element box [" << lower.transpose() << "] x ["
                    << upper.transpose() << "] does not match the Grid3 lattice at cell (" << i
                    << "," << j << "," << k << "); the background basis must share the lookup "
                    "grid's knot lines.");

        switch (m_idx->status[id])
        {
        case Full:
            m_full.mapTo(lower, upper, nodes, weights);
            break;
        case Empty:
            nodes.resize(3,0); weights.resize(0);
            break;
        default: // Cut
            m_src->volRule(id, nodes, weights);
            break;
        }
    }

private:
    memory::shared_ptr<const CellIndex> m_idx;
    memory::shared_ptr<const VolCellSource> m_src;
    gsGaussRule<real_t> m_full;
};

/// gsQuadRule<real_t> for the background grid's boundary integral, keyed
/// exactly as VolLookupRule. Calls \a src's bdrRule() for EVERY cell
/// status (the boundary rule does not depend on the volume classification
/// -- a Full cell can still carry a zero-column boundary rule; a Cut cell
/// always carries the clipped one). The weights returned are PHYSICAL
/// surface weights, so a Nitsche expression built on this rule carries NO
/// separate measure factor; the outward unit normals are NOT returned by
/// mapTo() (gsQuadRule's interface has no slot for them) -- they come from
/// BdrNormalField, registered on the same \a src.
class BdrLookupRule : public gsQuadRule<real_t>
{
public:
    BdrLookupRule(memory::shared_ptr<const CellIndex> idx,
                 memory::shared_ptr<const BdrCellSource> src)
    : m_idx(give(idx)), m_src(give(src))
    { }

    using gsQuadRule<real_t>::mapTo;

    void mapTo(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
              gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        GISMO_ENSURE(3 == lower.size(), "BdrLookupRule::mapTo: expects a 3D element box, got "
                    "lower.size()=" << lower.size());
        const gsVector<real_t> mid = 0.5*(lower+upper);
        const size_t id = m_idx->cellOfPoint(mid[0], mid[1], mid[2]);

        index_t i,j,k; m_idx->ijk(id, i,j,k);
        const real_t h = m_idx->grid.h;
        GISMO_ENSURE(math::abs(lower[0]-m_idx->X[i])   <= 1e-8*h && math::abs(upper[0]-m_idx->X[i+1]) <= 1e-8*h &&
                    math::abs(lower[1]-m_idx->Y[j])   <= 1e-8*h && math::abs(upper[1]-m_idx->Y[j+1]) <= 1e-8*h &&
                    math::abs(lower[2]-m_idx->Z[k])   <= 1e-8*h && math::abs(upper[2]-m_idx->Z[k+1]) <= 1e-8*h,
                    "BdrLookupRule::mapTo: the element box [" << lower.transpose() << "] x ["
                    << upper.transpose() << "] does not match the Grid3 lattice at cell (" << i
                    << "," << j << "," << k << "); the background basis must share the lookup "
                    "grid's knot lines.");

        gsMatrix<real_t> normals; // discarded: normals come from BdrNormalField
        m_src->bdrRule(id, nodes, weights, normals);
    }

private:
    memory::shared_ptr<const CellIndex> m_idx;
    memory::shared_ptr<const BdrCellSource> m_src;
};

//----------------------------------------------------------------------------
// Boundary normal field.
//----------------------------------------------------------------------------

/// gsFunction<real_t> returning the stored outward unit normals for
/// exactly the node set a BdrLookupRule (on the SAME \a src) produced for
/// one cell, and throwing on any other point set. domainDim()=3,
/// targetDim()=3; no deriv_into is provided (a per-cell CONSTANT-per-
/// triangle normal has no meaningful derivative here).
///
/// (i) Registration: use `A.getCoeff(nf)` / `ev.getVariable(nf)` with NO
/// geometry argument (gsExprAssembler.h:346-347, gsExprEvaluator.h:174-
/// 175). Such a variable is evaluated at exactly the active rule's own
/// nodes (gsExprHelper.h precompute -> compute(m_points, ...)), which is
/// what eval_into()'s bitwise node match below requires. The COMPOSED
/// `getCoeff(nf, G)` evaluates at G(u) instead (gsExprHelper.h, the
/// mutMap-composed compute() branch), which is not the rule's own node
/// set -- it WILL throw here.
/// (ii) In streamer mode (\a src a ClipStreamer) every evaluation re-clips
/// the candidate cell's triangle bucket to find and verify the matching
/// node set, about 2x the cost of one boundary clip; a BdrCellTable
/// avoids this by returning the stored rule directly.
/// (iii) Thread-safety: eval_into() is const and touches only \a m_idx and
/// \a m_src (both read-only); no mutable state (no "last cell served"
/// cache) is kept, so concurrent evaluation is safe.
class BdrNormalField : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(BdrNormalField)

    BdrNormalField(memory::shared_ptr<const CellIndex> idx,
                  memory::shared_ptr<const BdrCellSource> src)
    : m_idx(give(idx)), m_src(give(src))
    { }

    short_t domainDim() const override { return 3; }
    short_t targetDim() const override { return 3; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        if (0 == u.cols()) { result.resize(3,0); return; }
        GISMO_ENSURE(3 == u.rows(), "BdrNormalField::eval_into: expects 3D points, got "
                    "u.rows()=" << u.rows());

        // All nodes of one cell's boundary rule lie in that cell's closed
        // box up to a few ULPs of interpolation roundoff, so their
        // centroid does too -- but it can sit within ULPs of a grid line
        // (a sliver piece touching a cell face), so a bare colOf(centroid)
        // is not robust: widen by tol on every axis and check every
        // resulting candidate cell (at most 8) by an EXACT bitwise node
        // match, never by proximity.
        const gsVector<real_t> c = u.rowwise().mean();
        const real_t tol = 64*std::numeric_limits<real_t>::epsilon()
                          * math::max(boxCornerAbsMax(m_idx->grid), m_idx->grid.h);
        const index_t i0 = colOf(c[0]-tol, m_idx->X), i1 = colOf(c[0]+tol, m_idx->X);
        const index_t j0 = colOf(c[1]-tol, m_idx->Y), j1 = colOf(c[1]+tol, m_idx->Y);
        const index_t k0 = colOf(c[2]-tol, m_idx->Z), k1 = colOf(c[2]+tol, m_idx->Z);

        for (index_t k = k0; k <= k1; ++k)
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
        {
            const size_t id = m_idx->id(i,j,k);
            gsMatrix<real_t> n; gsVector<real_t> w; gsMatrix<real_t> nrm;
            m_src->bdrRule(id, n, w, nrm);
            if (n.cols() == u.cols() && n.cols() > 0 &&
                0 == std::memcmp(n.data(), u.data(), (size_t)(3*u.cols())*sizeof(real_t)))
            { result = nrm; return; }
        }

        GISMO_ERROR("BdrNormalField: the " << u.cols() << " evaluation points (centroid "
                   << c.transpose() << ") are not the nodes of any boundary cell rule; the "
                   "field is only valid on points produced by a BdrLookupRule on the same "
                   "source");
    }

private:
    memory::shared_ptr<const CellIndex> m_idx;
    memory::shared_ptr<const BdrCellSource> m_src;
};

//----------------------------------------------------------------------------
// QuadratureFactory builders. gsExprAssembler<real_t>::QuadratureFactory
// and gsExprEvaluator<real_t>::QuadratureFactory are the same std::function
// type (gsExprAssembler.h:91-97, gsExprEvaluator.h:73-79), so one factory
// value serves both.
//----------------------------------------------------------------------------

/// Builds a QuadratureFactory returning a fresh VolLookupRule per call,
/// backed by \a src and keyed by \a idx, with the FIXED Full-cell Gauss
/// order \a nGaussFull (see VolLookupRule's own doxygen). Rejects any call
/// made for a face/boundary integral (\a fixedDirection != -1): a volume
/// rule handed a face would silently integrate the wrong manifold
/// (gsExprAssembler.h:1714-1720).
///
/// Trap: gsExprAssembler::assemble() calls this factory (via
/// makeQuadratureRule) inside `#pragma omp parallel`
/// (gsExprAssembler.h:1431-1471); a GISMO_ENSURE firing there terminates
/// the process rather than propagating as a catchable exception.
/// assembleBdr() and the ghost/skeleton face loop are serial, so throws
/// from those call sites DO propagate normally.
inline gsExprAssembler<real_t>::QuadratureFactory
makeVolLookupFactory(memory::shared_ptr<const CellIndex> idx,
                     memory::shared_ptr<const VolCellSource> src,
                     const gsVector<index_t> & nGaussFull)
{
    return [idx, src, nGaussFull]
           (const gsDomain<real_t> &, const gsBasis<real_t> *, const gsOptionList &,
            index_t, short_t fixedDirection, const gsVector<short_t> &) -> gsQuadRule<real_t>::uPtr
    {
        GISMO_ENSURE(-1 == fixedDirection, "makeVolLookupFactory: a volume quadrature factory "
                    "was handed fixedDirection=" << fixedDirection << " (a face/boundary call); "
                    "a volume rule on a face integrates the wrong manifold.");
        return gsQuadRule<real_t>::uPtr(new VolLookupRule(idx, src, nGaussFull));
    };
}

/// Builds a QuadratureFactory returning a fresh BdrLookupRule per call.
/// Accepts ONLY \a fixedDirection == 0: an immersed
/// `assembleBdr(bContainer(1,patchSide(0,boundary::none)), ...)` call
/// passes `direction() = (m_index-1)/2 = 0` for `boundary::none = 0`
/// (gsExprAssembler.h:1565, gsBoundary.h:60,113), which this ENSURE also
/// catches for a y/z-direction ghost/skeleton face -- but an
/// x-DIRECTION ghost face is indistinguishable from the immersed boundary
/// by fixedDirection alone, so only the enclosing BoundaryQuadratureScope
/// (which must never span a ghost/skeleton face loop) protects against
/// that case; this factory cannot.
inline gsExprAssembler<real_t>::QuadratureFactory
makeBdrLookupFactory(memory::shared_ptr<const CellIndex> idx,
                     memory::shared_ptr<const BdrCellSource> src)
{
    return [idx, src]
           (const gsDomain<real_t> &, const gsBasis<real_t> *, const gsOptionList &,
            index_t, short_t fixedDirection, const gsVector<short_t> &) -> gsQuadRule<real_t>::uPtr
    {
        GISMO_ENSURE(0 == fixedDirection, "makeBdrLookupFactory: a boundary quadrature factory "
                    "was handed fixedDirection=" << fixedDirection << "; expected 0 (immersed "
                    "assembleBdr(patchSide(0,boundary::none))).");
        return gsQuadRule<real_t>::uPtr(new BdrLookupRule(idx, src));
    };
}

//----------------------------------------------------------------------------
// RAII quadrature-factory scopes.
//----------------------------------------------------------------------------

/// Tests whether an int option named \a label is present in \a o, using
/// only PUBLIC gsOptionList API: gsOptionList::exists() (gsOptionList.h)
/// is private, and getInt() throws (via GISMO_ENSURE, which also prints to
/// stderr) on an absent label -- unusable here since an absent "quDim" is
/// the NORMAL case on every scope entry, not an error to report. askInt()
/// is public and silent for an absent label (it only ever warns, and only
/// under GISMO_WITH_XDEBUG, if the label exists as a DIFFERENT option
/// type), so presence is tested with two sentinel defaults: a stored value
/// cannot equal both index_t extremes at once, so both askInt() calls
/// returning their own default means no such int option is stored.
inline bool hasIntOption(const gsOptionList & o, const std::string & label)
{
    const index_t lo = std::numeric_limits<index_t>::min(), hi = std::numeric_limits<index_t>::max();
    return !(o.askInt(label, lo) == lo && o.askInt(label, hi) == hi);
}

/// RAII scope that installs a quadrature factory on a gsExprAssembler or
/// gsExprEvaluator (\a ExprObj) and, on destruction, restores BOTH the
/// installed-factory state (cleared) and the "quDim" option (to whatever
/// it was before the scope -- present with its old value, or absent).
///
/// gsExprAssembler/gsExprEvaluator expose NO factory getter (only
/// setQuadratureFactory/clearQuadratureFactory/hasCustomQuadrature), so a
/// pre-existing factory cannot be saved: entering a scope while one is
/// already installed (a nested scope, or a leaked factory from a scope
/// whose destructor never ran) is a programming error the constructor
/// ENSUREs against, rather than one this class could silently repair.
///
/// A scope must NEVER span a ghost/skeleton face loop:
/// gsExprAssembler::assembleGhost/assembleSkeleton use the installed
/// factory, if any, for EVERY parametric direction of a face
/// (gsExprAssembler.h:1721-1725), not only the manifold the scope itself
/// was built for. A gsExprEvaluator built from a gsExprAssembler does NOT
/// copy the assembler's installed factory (gsExprEvaluator.h's converting
/// constructor copies only exprData and options), so it needs its own
/// scope even while sharing the assembler's exprData.
template<class ExprObj>
class QuadratureScope
{
public:
    QuadratureScope(ExprObj & obj, typename ExprObj::QuadratureFactory factory, index_t quDim)
    : m_obj(obj), m_hadQuDim(hasIntOption(obj.options(), "quDim")),
      m_oldQuDim(m_hadQuDim ? obj.options().getInt("quDim") : -1)
    {
        GISMO_ENSURE(!m_obj.hasCustomQuadrature(), "QuadratureScope: a factory is already "
                    "installed (nested scope or leaked factory); gsExprAssembler/"
                    "gsExprEvaluator expose no factory getter, so it could not be restored.");
        GISMO_ENSURE(static_cast<bool>(factory), "QuadratureScope: empty factory.");
        m_obj.options().addInt("quDim", "Quadrature manifold: -1 volume, >=0 surface", quDim);
        m_obj.setQuadratureFactory(give(factory));
    }

    ~QuadratureScope()
    {
        m_obj.clearQuadratureFactory();
        if (m_hadQuDim) m_obj.options().setInt("quDim", m_oldQuDim);
        else            m_obj.options().remove("quDim");
    }

    QuadratureScope(const QuadratureScope &) = delete;
    QuadratureScope & operator=(const QuadratureScope &) = delete;

private:
    ExprObj & m_obj;
    const bool m_hadQuDim;
    const index_t m_oldQuDim;
};

/// QuadratureScope with quDim = -1 (volume quadrature).
template<class ExprObj>
class VolumeQuadratureScope : public QuadratureScope<ExprObj>
{
public:
    VolumeQuadratureScope(ExprObj & obj, typename ExprObj::QuadratureFactory factory)
    : QuadratureScope<ExprObj>(obj, give(factory), -1)
    { }
};

/// QuadratureScope with quDim = 2 (surface quadrature; any quDim >= 0
/// selects surface quadrature in gsQuadrature's own askInt("quDim",-1) --
/// 2 follows the existing convention of the 2D immersed drivers).
template<class ExprObj>
class BoundaryQuadratureScope : public QuadratureScope<ExprObj>
{
public:
    BoundaryQuadratureScope(ExprObj & obj, typename ExprObj::QuadratureFactory factory)
    : QuadratureScope<ExprObj>(obj, give(factory), 2)
    { }
};

} // namespace gsTetClip
