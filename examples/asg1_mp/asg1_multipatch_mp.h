/** @file asg1_multipatch_mp.h

    @brief Header-only AS-G1 conforming multi-patch basis assembler.

    Generalises the two-patch v4 construction
    (`examples/asg1_v4/as_g1_two_patch_basis_v4.h`) to an arbitrary
    planar multi-patch geometry following the Kapl-Sangalli-Takacs
    decomposition

        Global  =  Patch-Interior  +  Edge  ( +  Vertex ).

    Output: one global-to-patch sparse matrix `G[p]` per patch,
        c_p = G[p] * globalCoefficientVector,
    with the column layout

        [ interior(p0) | ... | interior(pP-1)
        | edge0_trace | edge0_deriv | edge1_trace | edge1_deriv | ... ]
        ( | vertex0(6) | ... )         <-- added in the vertex phase

    Interior columns are identity blocks; each interior edge contributes
    shared trace and d-derivative columns written into BOTH adjacent
    patches with the same `flipped`-aware sign/index logic as v4.

    This first iteration implements the Interior + Edge spaces (no
    edge-end truncation, no vertex space yet), which is exact for any
    multi-patch geometry whose interfaces do not meet at a shared
    interior corner (e.g. a 1xN patch strip).  The vertex space and
    edge-end truncation are layered on top in later phases.

    Public API (in namespace `gismo::asg1mp`):

      * struct `AsG1MpOptions<T>`
      * struct `AsG1Multipatch<T>`
      * `buildAsG1Multipatch(mp, opts)`
      * `smoothnessCheckMP(B, nCheck)`

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include <gismo.h>

#include "topology_mp.h"
#include "gluing_data_mp.h"
#include "../asg1_v4/argyris_embedding_v4.h"

#include <gsMSplines/gsMappedBasis.h>

namespace gismo {
namespace asg1mp {

// ====================================================================
// gsVertexBasisMP -- analytic AS-G1 vertex shape function (one bfID).
//
// Ported from gsUnstructuredSplines (gsApproxC1Utils.h, gsVertexBasis,
// P. Weinmueller & A. Farahat).  Evaluates, in the *standard* frame
// (vertex at (0,0), incident edges along +u/+v), the closed-form
// Kapl-Sangalli-Takacs vertex basis function for a prescribed 6-jet
// column `bfID` of `Phi`.  We feed it *exact* linear alpha/beta from
// the (reparametrized) geometry, so for piecewise-bilinear patches the
// result is an exact degree-p spline recovered exactly by
// interpolation at anchors (no L2 fit -> exact and fast).
// ====================================================================
template <class T>
class gsVertexBasisMP : public gismo::gsFunction<T>
{
protected:
    const gsGeometry<T> &           m_geo;
    gsBasis<T> &                    m_basis;
    std::vector<gsBSpline<T> >      m_alpha;
    std::vector<gsBSpline<T> >      m_beta;
    std::vector<gsBSplineBasis<T> > m_basis_plus;
    std::vector<gsBSplineBasis<T> > m_basis_minus;
    const gsMatrix<T>               m_Phi;
    const std::vector<bool>         m_kindOfEdge;
    const index_t                   m_bfID;

public:
    gsVertexBasisMP(const gsGeometry<T> & geo, gsBasis<T> & basis,
                    std::vector<gsBSpline<T> > alpha, std::vector<gsBSpline<T> > beta,
                    std::vector<gsBSplineBasis<T> > basis_plus,
                    std::vector<gsBSplineBasis<T> > basis_minus,
                    const gsMatrix<T> Phi, const std::vector<bool> kindOfEdge,
                    const index_t bfID)
        : m_geo(geo), m_basis(basis), m_alpha(alpha), m_beta(beta),
          m_basis_plus(basis_plus), m_basis_minus(basis_minus),
          m_Phi(Phi), m_kindOfEdge(kindOfEdge), m_bfID(bfID),
          _vertexBasis_piece(nullptr) {}

    ~gsVertexBasisMP() { delete _vertexBasis_piece; }
    GISMO_CLONE_FUNCTION(gsVertexBasisMP)
    short_t domainDim() const { return 2; }
    short_t targetDim() const { return 1; }

    mutable gsVertexBasisMP<T> * _vertexBasis_piece;
    const gsFunction<T> & piece(const index_t k) const
    {
        GISMO_UNUSED(k);
        _vertexBasis_piece = new gsVertexBasisMP(*this);
        return *_vertexBasis_piece;
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(targetDim(), u.cols());
        result.setZero();

        gsMatrix<T> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<T> > c_0, c_1;
        std::vector<gsMatrix<T> > c_0_plus, c_1_plus, c_2_plus;
        std::vector<gsMatrix<T> > c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
        std::vector<gsMatrix<T> > c_0_minus, c_1_minus;
        for (index_t i = 0; i < 2; i++)
        {
            gsMatrix<T> b_0, b_1;
            gsMatrix<T> b_0_plus, b_1_plus, b_2_plus;
            gsMatrix<T> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
            gsMatrix<T> b_0_minus, b_1_minus;

            m_basis.component(i).evalSingle_into(0, u.row(i), b_0);
            m_basis.component(i).evalSingle_into(1, u.row(i), b_1);

            m_basis_plus[i].evalSingle_into(0, u.row(i), b_0_plus);
            m_basis_plus[i].evalSingle_into(1, u.row(i), b_1_plus);
            m_basis_plus[i].evalSingle_into(2, u.row(i), b_2_plus);

            m_basis_plus[i].derivSingle_into(0, u.row(i), b_0_plus_deriv);
            m_basis_plus[i].derivSingle_into(1, u.row(i), b_1_plus_deriv);
            m_basis_plus[i].derivSingle_into(2, u.row(i), b_2_plus_deriv);

            m_basis_minus[i].evalSingle_into(0, u.row(i), b_0_minus);
            m_basis_minus[i].evalSingle_into(1, u.row(i), b_1_minus);

            gsMatrix<T> b_1_0, b_1_minus_0;
            m_basis.component(i).derivSingle_into(1, zero.row(i), b_1_0);
            m_basis_minus[i].derivSingle_into(1, zero.row(i), b_1_minus_0);

            T factor_b_1 = 1.0/b_1_0(0,0);
            c_0.push_back(b_0 + b_1);
            c_1.push_back(factor_b_1 * b_1);

            T factor_b_1_minus = 1.0/b_1_minus_0(0,0);
            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(factor_b_1_minus * b_1_minus);

            gsMatrix<T> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
            m_basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

            T factor_c_1_plus = 1.0/der_b_1_plus_0(0,0);
            T factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
            T factor_c_2_plus = 1.0/der2_b_2_plus_0(0,0);

            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
            c_2_plus.push_back(factor_c_2_plus * b_2_plus);

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
            c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
        }

        std::vector<gsMatrix<T> > alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;
        gsMatrix<T> temp_mat;
        for (index_t i = 0; i < 2; ++i)
        {
            if (m_kindOfEdge[i])
            {
                m_alpha[i].eval_into(u.row(i), temp_mat);   alpha.push_back(temp_mat);
                m_alpha[i].eval_into(zero.row(0), temp_mat); alpha_0.push_back(temp_mat);
                m_alpha[i].deriv_into(zero.row(0), temp_mat); alpha_deriv.push_back(temp_mat);
                m_beta[i].eval_into(u.row(i), temp_mat);    beta.push_back(temp_mat);
                m_beta[i].eval_into(zero.row(0), temp_mat); beta_0.push_back(temp_mat);
                m_beta[i].deriv_into(zero.row(0), temp_mat); beta_deriv.push_back(temp_mat);
            }
            else
            {
                temp_mat.setOnes(1, u.cols());  alpha.push_back(temp_mat);
                temp_mat.setOnes(1, 1);         alpha_0.push_back(temp_mat);
                temp_mat.setZero(1, 1);         alpha_deriv.push_back(temp_mat);
                temp_mat.setZero(1, u.cols());  beta.push_back(temp_mat);
                temp_mat.setZero(1, 1);         beta_0.push_back(temp_mat);
                temp_mat.setZero(1, 1);         beta_deriv.push_back(temp_mat);
            }
        }

        gsMatrix<T> geo_jac  = m_geo.jacobian(zero);
        gsMatrix<T> geo_der2 = m_geo.deriv2(zero);

        gsMatrix<T> dd_ik_plus, dd_ik_minus, dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1.0/(alpha_0[0](0,0)) * (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0));
        dd_ik_plus  =  1.0/(alpha_0[1](0,0)) * (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1));

        gsMatrix<T> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
        geo_deriv2_12.row(0) = geo_der2.row(2); geo_deriv2_12.row(1) = geo_der2.row(5);
        geo_deriv2_11.row(0) = geo_der2.row(0); geo_deriv2_11.row(1) = geo_der2.row(3);
        geo_deriv2_22.row(0) = geo_der2.row(1); geo_deriv2_22.row(1) = geo_der2.row(4);

        gsMatrix<T> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<T> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1.0/(alpha_squared_u(0,0)) *
            ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) + beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
             (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0)) * alpha_deriv[0](0,0));
        dd_ik_plus_deriv = 1.0/(alpha_squared_v(0,0)) *
            ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) + beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
             (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1)) * alpha_deriv[1](0,0));

        std::vector<gsMatrix<T> > d_ik;
        d_ik.push_back(m_Phi.col(0));
        d_ik.push_back(m_Phi.block(1,0,2,6).transpose() * geo_jac.col(0));
        d_ik.push_back(m_Phi.block(1,0,2,6).transpose() * geo_jac.col(1));
        d_ik.push_back((geo_jac(0,0) * m_Phi.col(3) + geo_jac(1,0) * m_Phi.col(4)) * geo_jac(0,1) +
                       (geo_jac(0,0) * m_Phi.col(4) + geo_jac(1,0) * m_Phi.col(5)) * geo_jac(1,1) +
                       m_Phi.block(0,1,6,1) * geo_der2.row(2) +
                       m_Phi.block(0,2,6,1) * geo_der2.row(5));

        std::vector<gsMatrix<T> > d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(m_Phi.col(0));
        d_ilik_minus.push_back(m_Phi.block(1,0,2,6).transpose() * geo_jac.col(0));
        d_ilik_minus.push_back((geo_jac(0,0) * m_Phi.col(3) + geo_jac(1,0) * m_Phi.col(4)) * geo_jac(0,0) +
                               (geo_jac(0,0) * m_Phi.col(4) + geo_jac(1,0) * m_Phi.col(5)) * geo_jac(1,0) +
                               m_Phi.block(0,1,6,1) * geo_der2.row(0) +
                               m_Phi.block(0,2,6,1) * geo_der2.row(3));
        d_ilik_minus.push_back(m_Phi.block(1,0,2,6).transpose() * dd_ik_minus);
        d_ilik_minus.push_back((geo_jac(0,0) * m_Phi.col(3) + geo_jac(1,0) * m_Phi.col(4)) * dd_ik_minus(0,0) +
                               (geo_jac(0,0) * m_Phi.col(4) + geo_jac(1,0) * m_Phi.col(5)) * dd_ik_minus(1,0) +
                               m_Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                               m_Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

        d_ilik_plus.push_back(m_Phi.col(0));
        d_ilik_plus.push_back(m_Phi.block(1,0,2,6).transpose() * geo_jac.col(1));
        d_ilik_plus.push_back((geo_jac(0,1) * m_Phi.col(3) + geo_jac(1,1) * m_Phi.col(4)) * geo_jac(0,1) +
                              (geo_jac(0,1) * m_Phi.col(4) + geo_jac(1,1) * m_Phi.col(5)) * geo_jac(1,1) +
                              m_Phi.block(0,1,6,1) * geo_der2.row(1) +
                              m_Phi.block(0,2,6,1) * geo_der2.row(4));
        d_ilik_plus.push_back(m_Phi.block(1,0,2,6).transpose() * dd_ik_plus);
        d_ilik_plus.push_back((geo_jac(0,1) * m_Phi.col(3) + geo_jac(1,1) * m_Phi.col(4)) * dd_ik_plus(0,0) +
                              (geo_jac(0,1) * m_Phi.col(4) + geo_jac(1,1) * m_Phi.col(5)) * dd_ik_plus(1,0) +
                              m_Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                              m_Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));

        result = d_ilik_minus.at(0)(m_bfID,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                     beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                 d_ilik_minus.at(1)(m_bfID,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                     beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                 d_ilik_minus.at(2)(m_bfID,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                     beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                 d_ilik_minus.at(3)(m_bfID,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                 d_ilik_minus.at(4)(m_bfID,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1)));

        result += d_ilik_plus.at(0)(m_bfID,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                      beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                  d_ilik_plus.at(1)(m_bfID,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                      beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                  d_ilik_plus.at(2)(m_bfID,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                      beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                  d_ilik_plus.at(3)(m_bfID,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                  d_ilik_plus.at(4)(m_bfID,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0)));

        result -= d_ik.at(0)(m_bfID,0) * c_0.at(0).cwiseProduct(c_0.at(1)) +
                  d_ik.at(2)(m_bfID,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                  d_ik.at(1)(m_bfID,0) * c_1.at(0).cwiseProduct(c_0.at(1)) +
                  d_ik.at(3)(m_bfID,0) * c_1.at(0).cwiseProduct(c_1.at(1));
    }
};

// ====================================================================
// Options
// ====================================================================

template <class T>
struct AsG1MpOptions
{
    index_t      refinements    = 0;
    short_t      minDegree      = 3;
    index_t      numGaussPerSpan = 0;
    T            eps            = T(1e-8);
    // Gluing data is computed on the actual geometry (Direct) -- the same
    // approach as the original two-patch code.  This is exact for
    // piecewise-bilinear patches and the correct (approximate) AS-G1 data
    // for curved patches.  BilinearTemplate is kept only as an option for
    // non-AS-G1 input that needs a vertex-compatible template; the vertex
    // space here does NOT use the stored gluing data, so Direct is the
    // right default.
    GluingSource gluingSource   = GluingSource::Direct;
    bool         truncate       = true;  ///< Phase D edge-end truncation
    bool         vertex         = true;  ///< Phase E interior-vertex space
    index_t      vtxWindow      = 4;     ///< corner block size per direction
    index_t      vtxSpanCover   = 2;     ///< tangent knot spans to constrain
    index_t      vtxSamples     = 24;    ///< samples per incident edge
    T            vtxNullTol     = T(1e-7); ///< rel. SVD null-space threshold
    bool         verbose        = true;
};

// Number of near-vertex edge functions removed at an interior-vertex end.
static const index_t ASG1_TRACE_DROP = 3;   ///< smoother (trace) cols / end
static const index_t ASG1_DERIV_DROP = 2;   ///< lower-degree (deriv) cols / end

// ====================================================================
// Per-edge embedding bundle
// ====================================================================

template <class T>
struct EdgeEmbedding
{
    index_t            edgeId = -1;
    gsSparseMatrix<T>  E1, E2;          ///< per-side embeddings
    index_t            nInt1 = 0, nInt2 = 0;
    index_t            nLD = 0, nSm = 0; ///< full shared deriv / trace counts
    // Phase D truncation: number of cols dropped at the low / high
    // tangent end of side1 (an end is dropped iff its vertex is interior).
    index_t            traceLo = 0, traceHi = 0;
    index_t            derivLo = 0, derivHi = 0;
    index_t            nSmKept = 0, nLDKept = 0;
    index_t            gOffTrace = 0;    ///< global column of first kept trace dof
    index_t            gOffDeriv = 0;    ///< global column of first kept deriv dof
    bool               flipped = false;
};

// ====================================================================
// Output bundle
// ====================================================================

/// Six (in general) closed-support vertex basis functions for one
/// interior vertex, expressed over a set of local patch DOFs.
template <class T>
struct VertexFuncs
{
    index_t                                 vertexId = -1;
    std::vector<std::pair<index_t,index_t> > localDofs; ///< (patch, tensorDof)
    gsMatrix<T>                             W;          ///< nLocal x k
    index_t                                 gColBase = 0;
};

template <class T>
struct AsG1Multipatch
{
    gsMultiPatch<T>                 mp;       ///< refined / elevated copy
    std::vector<EdgeInfo<T> >       edges;
    std::vector<VertexInfo<T> >     verts;
    std::vector<EdgeEmbedding<T> >  edgeEmb;  ///< one per interior edge

    std::vector<gsSparseMatrix<T> > G;        ///< global-to-patch, size nPatches
    index_t                         nGlobal = 0;

    // global column offsets
    std::vector<index_t>            interiorOff;   ///< nPatches+1 (prefix sums)
    index_t                         edgeBlockOff = 0;
    index_t                         vertexBlockOff = 0;

    // per-patch interior DOF index lists (into the tensor basis)
    std::vector<std::vector<index_t> > interiorDofs;

    // Phase E: one entry per interior vertex that received a vertex block
    std::vector<VertexFuncs<T> >    vertexFuncs;
};

// ====================================================================
// Helpers
// ====================================================================

/// Tensor B-spline basis of patch p (throws if not a tensor B-spline).
template <class T>
inline const gsTensorBSplineBasis<2,T>& asTB(const gsMultiPatch<T>& mp, index_t p)
{
    return dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(p).basis());
}

/// 2-D corner index (1..4) at the low/high tangent end of a side:
///   1=(0,0) 2=(1,0) 3=(0,1) 4=(1,1).
/// `tdEnd` is 0 for the low tangent end, 1 for the high end.
inline index_t cornerAtSideEnd(const patchSide& s, int tdEnd)
{
    const short_t d = s.direction();          // normal direction
    const bool    p = s.parameter();          // fixed coordinate (0/1)
    index_t upar, vpar;
    if (d == 0) { upar = p ? 1 : 0; vpar = tdEnd; }   // west/east
    else        { vpar = p ? 1 : 0; upar = tdEnd; }   // south/north
    return 1 + upar + 2 * vpar;
}

/// Is the vertex at the given (low/high) tangent end of side `s` an
/// interior (extraordinary) vertex?  Looks the corner up in `verts`.
template <class T>
bool sideEndIsInterior(
    const patchSide& s, int tdEnd,
    const std::vector<VertexInfo<T> >& verts)
{
    const index_t cIdx = cornerAtSideEnd(s, tdEnd);
    for (const VertexInfo<T>& v : verts)
        for (const patchCorner& c : v.corners)
            if (c.patch == s.patch && c.m_index == cIdx)
                return !v.isBoundary;
    return false;   // unknown -> treat as boundary (keep cols)
}

// ====================================================================
// Phase E -- interior vertex space (local null-space construction)
// ====================================================================

/// Tensor DOF indices of the `m x m` corner block at corner `cIdx`
/// (1..4) of a 2-D tensor B-spline basis.
template <class T>
void cornerBlockDofs(const gsTensorBSplineBasis<2,T>& tb, index_t cIdx,
                     index_t m, std::vector<index_t>& dofs)
{
    const index_t n0 = tb.size(0), n1 = tb.size(1);
    const bool pu = ((cIdx - 1) & 1);
    const bool pv = (((cIdx - 1) >> 1) & 1);
    const index_t mu = std::min(m, n0), mv = std::min(m, n1);
    dofs.clear();
    for (index_t dj = 0; dj < mv; ++dj)
        for (index_t di = 0; di < mu; ++di)
        {
            const index_t i = pu ? (n0 - 1 - di) : di;
            const index_t j = pv ? (n1 - 1 - dj) : dj;
            dofs.push_back(i + j * n0);
        }
}

/// Corner index (1..4) at which patch `p` touches vertex `vi`, or -1.
template <class T>
index_t cornerIdxForPatch(const VertexInfo<T>& vi, index_t p)
{
    for (const patchCorner& c : vi.corners)
        if (c.patch == p) return c.m_index;
    return -1;
}

/// Build the vertex basis functions for one interior vertex.
///
/// The functions are sought as combinations of the *dropped* (by Phase D
/// truncation) near-vertex edge columns of every incident interface.
/// Each dropped column already extends along its own edge with the
/// correct gluing-data Hermite shape and is G1 across that edge; we look
/// for the combinations that are *also* G1 across all the other incident
/// edges (and single valued at the corner).  This captures both the
/// strictly-localized and the edge-extending vertex shape functions.
template <class T>
VertexFuncs<T> buildVertexFuncs(
    const gsMultiPatch<T>& mp,
    const std::vector<EdgeInfo<T> >& edges,
    const std::vector<EdgeEmbedding<T> >& edgeEmb,
    const VertexInfo<T>& vi,
    const AsG1MpOptions<T>& opts)
{
    VertexFuncs<T> vf;
    vf.vertexId = vi.id;

    // edgeId -> embedding
    std::map<index_t, const EdgeEmbedding<T>*> embOf;
    for (const EdgeEmbedding<T>& em : edgeEmb) embOf[em.edgeId] = &em;

    // --- gather incident interface edges ---
    std::vector<const EdgeInfo<T>*> incident;
    for (const EdgeInfo<T>& e : edges)
    {
        if (e.isBoundary) continue;
        const index_t ca = cornerIdxForPatch(vi, e.side1.patch);
        const index_t cb = cornerIdxForPatch(vi, e.side2.patch);
        if (ca < 0 || cb < 0) continue;
        const bool onA = (cornerAtSideEnd(e.side1,0)==ca) || (cornerAtSideEnd(e.side1,1)==ca);
        const bool onB = (cornerAtSideEnd(e.side2,0)==cb) || (cornerAtSideEnd(e.side2,1)==cb);
        if (onA && onB) incident.push_back(&e);
    }

    // --- candidate functions = dropped near-vertex columns ---
    typedef std::vector<std::pair<std::pair<index_t,index_t>, T> > CandFun; // ((patch,dof),val)
    std::vector<CandFun> cand;

    for (const EdgeInfo<T>* ep : incident)
    {
        const EdgeInfo<T>& e = *ep;
        auto itEm = embOf.find(e.id);
        if (itEm == embOf.end()) continue;
        const EdgeEmbedding<T>& em = *itEm->second;
        const index_t p1 = e.side1.patch, p2 = e.side2.patch;
        const index_t cv1 = cornerIdxForPatch(vi, p1);
        const int e1end = (cornerAtSideEnd(e.side1,0)==cv1) ? 0 : 1;

        // dropped index ranges at this vertex end (side1 indexing)
        const index_t tLo = (e1end==0) ? em.traceLo : 0;
        const index_t tHi = (e1end==1) ? em.traceHi : 0;
        const index_t dLo = (e1end==0) ? em.derivLo : 0;
        const index_t dHi = (e1end==1) ? em.derivHi : 0;
        std::vector<index_t> traceJ, derivJ;
        for (index_t j = 0; j < tLo; ++j)               traceJ.push_back(j);
        for (index_t j = em.nSm - tHi; j < em.nSm; ++j)  traceJ.push_back(j);
        for (index_t j = 0; j < dLo; ++j)               derivJ.push_back(j);
        for (index_t j = em.nLD - dHi; j < em.nLD; ++j)  derivJ.push_back(j);

        for (index_t jFull : traceJ)
        {
            CandFun f;
            const index_t e1col = em.nInt1 + em.nLD + jFull;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E1, e1col); it; ++it)
                f.push_back({{p1, it.row()}, it.value()});
            const index_t j2 = em.flipped ? (em.nSm - 1 - jFull) : jFull;
            const index_t e2col = em.nInt2 + em.nLD + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E2, e2col); it; ++it)
                f.push_back({{p2, it.row()}, it.value()});
            cand.push_back(f);
        }
        for (index_t jFull : derivJ)
        {
            CandFun f;
            const index_t e1col = em.nInt1 + jFull;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E1, e1col); it; ++it)
                f.push_back({{p1, it.row()}, it.value()});
            const index_t j2 = em.flipped ? (em.nLD - 1 - jFull) : jFull;
            const index_t e2col = em.nInt2 + j2;
            const T l2Sign = em.flipped ? T(1) : T(-1);
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E2, e2col); it; ++it)
                f.push_back({{p2, it.row()}, l2Sign * it.value()});
            cand.push_back(f);
        }
    }

    const index_t nCand = static_cast<index_t>(cand.size());
    if (nCand == 0) { vf.W = gsMatrix<T>(0,0); return vf; }

    // --- constraint matrix: physical C0 + G1 across every incident edge ---
    const index_t ns = opts.vtxSamples;
    const index_t nRows = 3 * ns * static_cast<index_t>(incident.size());
    gsMatrix<T> C = gsMatrix<T>::Zero(nRows, nCand);

    auto endRange = [&](const patchSide& s, index_t cv,
                        short_t& nd, short_t& td, T& ifcCoord, T& t0, T& tC)
    {
        nd = s.direction(); td = 1 - nd;
        gsMatrix<T> sup = mp.patch(s.patch).support();
        ifcCoord = sup(nd, s.parameter() ? 1 : 0);
        gsKnotVector<T> kv = asTB(mp, s.patch).component(td).knots();
        std::vector<T> uk = kv.unique();
        const int e = (cornerAtSideEnd(s,0)==cv) ? 0 : 1;
        const index_t K = std::min<index_t>(opts.vtxSpanCover,
                                            static_cast<index_t>(uk.size())-1);
        if (e == 0) { t0 = uk.front(); tC = uk[K]; }
        else        { t0 = uk.back();  tC = uk[uk.size()-1-K]; }
    };

    index_t rowBase = 0;
    for (const EdgeInfo<T>* ep : incident)
    {
        const EdgeInfo<T>& e = *ep;
        const index_t pa = e.side1.patch, pb = e.side2.patch;
        const index_t ca = cornerIdxForPatch(vi, pa);
        const index_t cb = cornerIdxForPatch(vi, pb);
        const gsTensorBSplineBasis<2,T>& tba = asTB(mp, pa);
        const gsTensorBSplineBasis<2,T>& tbb = asTB(mp, pb);

        short_t nda, tda, ndb, tdb;
        T ifcA, t0a, tCa, ifcB, t0b, tCb;
        endRange(e.side1, ca, nda, tda, ifcA, t0a, tCa);
        endRange(e.side2, cb, ndb, tdb, ifcB, t0b, tCb);

        for (index_t si = 0; si < ns; ++si)
        {
            const T sp = (ns > 1) ? T(si) / T(ns - 1) : T(0);
            const T ta = t0a + sp * (tCa - t0a);
            const T tb = t0b + sp * (tCb - t0b);
            gsMatrix<T> pA(2,1), pB(2,1);
            pA(nda,0) = ifcA; pA(tda,0) = ta;
            pB(ndb,0) = ifcB; pB(tdb,0) = tb;

            gsMatrix<T> dGa, dGb;
            mp.patch(pa).deriv_into(pA, dGa);
            mp.patch(pb).deriv_into(pB, dGb);
            gsMatrix<T> Ja(2,2), Jb(2,2);
            Ja(0,0)=dGa(0,0); Ja(0,1)=dGa(1,0); Ja(1,0)=dGa(2,0); Ja(1,1)=dGa(3,0);
            Jb(0,0)=dGb(0,0); Jb(0,1)=dGb(1,0); Jb(1,0)=dGb(2,0); Jb(1,1)=dGb(3,0);
            gsMatrix<T> JaIT = Ja.inverse().transpose();
            gsMatrix<T> JbIT = Jb.inverse().transpose();

            const index_t rVal = rowBase + 3*si, rGx = rVal+1, rGy = rVal+2;

            for (index_t i = 0; i < nCand; ++i)
            {
                T valA=0, gxA=0, gyA=0, valB=0, gxB=0, gyB=0;
                for (const auto& ent : cand[i])
                {
                    const index_t pp = ent.first.first, dof = ent.first.second;
                    const T w = ent.second;
                    if (pp == pa)
                    {
                        valA += w * tba.evalSingle(dof, pA)(0,0);
                        gsMatrix<T> dd = tba.derivSingle(dof, pA);
                        gsVector<T> gp(2); gp(0)=dd(0,0); gp(1)=dd(1,0);
                        gsVector<T> g = JaIT * gp; gxA += w*g(0); gyA += w*g(1);
                    }
                    if (pp == pb)
                    {
                        valB += w * tbb.evalSingle(dof, pB)(0,0);
                        gsMatrix<T> dd = tbb.derivSingle(dof, pB);
                        gsVector<T> gp(2); gp(0)=dd(0,0); gp(1)=dd(1,0);
                        gsVector<T> g = JbIT * gp; gxB += w*g(0); gyB += w*g(1);
                    }
                }
                C(rVal,i) += valA - valB;
                C(rGx, i) += gxA  - gxB;
                C(rGy, i) += gyA  - gyB;
            }
        }
        rowBase += 3 * ns;
    }

    // --- null space of constraints over candidate coefficients ---
    typename gsMatrix<T>::JacobiSVD svd(C, gsEigen::ComputeFullV);
    gsMatrix<T> V = svd.matrixV();
    gsVector<T> sv = svd.singularValues();
    const T smax = (sv.size() > 0) ? sv(0) : T(0);
    const T tol  = opts.vtxNullTol * std::max(smax, T(1));
    index_t rank = 0;
    for (index_t i = 0; i < sv.size(); ++i) if (sv(i) > tol) ++rank;
    const index_t nNull = nCand - rank;

    // --- realize each null combination as patch coefficients ---
    std::map<std::pair<index_t,index_t>, index_t> localIdx;
    std::vector<std::map<index_t,T> > cols(std::max<index_t>(nNull,0)); // localdof->val
    for (index_t q = 0; q < nNull; ++q)
    {
        std::map<std::pair<index_t,index_t>, T> acc;
        for (index_t i = 0; i < nCand; ++i)
        {
            const T x = V(i, rank + q);
            if (std::abs(x) < T(1e-14)) continue;
            for (const auto& ent : cand[i])
                acc[ent.first] += x * ent.second;
        }
        for (const auto& kv : acc)
        {
            if (std::abs(kv.second) < T(1e-12)) continue;
            std::map<std::pair<index_t,index_t>,index_t>::iterator f = localIdx.find(kv.first);
            index_t li;
            if (f == localIdx.end())
            {
                li = static_cast<index_t>(vf.localDofs.size());
                localIdx[kv.first] = li;
                vf.localDofs.push_back(kv.first);
            }
            else li = f->second;
            cols[q][li] = kv.second;
        }
    }

    const index_t nLocal = static_cast<index_t>(vf.localDofs.size());
    vf.W = gsMatrix<T>::Zero(nLocal, std::max<index_t>(nNull,0));
    for (index_t q = 0; q < nNull; ++q)
        for (const auto& kv : cols[q])
            vf.W(kv.first, q) = kv.second;

    if (opts.verbose)
        gsInfo << "  vertex " << vi.id << " (valence " << vi.valence
               << "): " << nNull << " functions from " << nCand
               << " dropped columns\n";
    return vf;
}



/// Collect, for patch p, the set of DOFs lying on any *interface*
/// first or second layer (these are owned by edges, not the interior).
template <class T>
std::set<index_t> interfaceLayerDofs(
    const gsMultiPatch<T>& mp,
    const std::vector<EdgeInfo<T> >& edges,
    index_t p)
{
    const gsTensorBSplineBasis<2,T>& tb = asTB(mp, p);
    std::set<index_t> layer;
    for (const EdgeInfo<T>& e : edges)
    {
        if (e.isBoundary) continue;
        for (int which = 0; which < 2; ++which)
        {
            const patchSide ps = (which == 0) ? e.side1 : e.side2;
            if (ps.patch != p) continue;
            gsMatrix<index_t> b0 = tb.boundary(ps.side());
            gsMatrix<index_t> b1 = tb.boundaryOffset(ps.side(), 1);
            for (index_t i = 0; i < b0.rows(); ++i) layer.insert(b0(i,0));
            for (index_t i = 0; i < b1.rows(); ++i) layer.insert(b1(i,0));
        }
    }
    return layer;
}

// ====================================================================
// Phase E -- reparametrization + KST vertex driver
// ====================================================================

/// Reparametrize `geo` so the given corner `cIdx` (1..4) maps to the
/// standard frame (vertex at (0,0), incident edges along +u, +v with a
/// consistent orientation), mirroring gsUnstructuredSplines'
/// gsG1AuxiliaryPatch.  Returns the rotated patch `geoRot` and a
/// permutation `perm` with `perm[rotDof] = origDof`.
template <class T>
void reparamCorner(const gsTensorBSpline<2,T>& geo, index_t cIdx,
                   gsTensorBSpline<2,T>& geoRot, std::vector<index_t>& perm)
{
    gsKnotVector<T> kvU = geo.basis().knots(0);
    gsKnotVector<T> kvV = geo.basis().knots(1);
    index_t n0 = geo.basis().size(0), n1 = geo.basis().size(1);
    const index_t N = n0 * n1;
    // M holds [x, y, origIndex] per row (tensor lexicographic i + j*n0).
    gsMatrix<T> M(N, 3);
    M.leftCols(2) = geo.coefs();
    for (index_t r = 0; r < N; ++r) M(r,2) = T(r);

    auto antiClock = [&]()
    {
        gsMatrix<T> Mn(N,3);
        for (index_t j = 0; j < n0; ++j)
            for (index_t i = 0; i < n1; ++i)
                Mn.row(i + j*n1) = M.row((n0-1-j) + n0*i);
        M = Mn; std::swap(kvU,kvV); std::swap(n0,n1);
    };
    auto clock = [&]()
    {
        gsMatrix<T> Mn(N,3);
        for (index_t j = n0-1; j >= 0; --j)
            for (index_t i = 0; i < n1; ++i)
                Mn.row(i + (n0-1-j)*n1) = M.row((n1*n0-1-j) - n0*i);
        M = Mn; std::swap(kvU,kvV); std::swap(n0,n1);
    };
    auto twice = [&]()
    {
        gsMatrix<T> Mn(N,3);
        for (index_t r = 0; r < N; ++r) Mn.row(r) = M.row(N-1-r);
        M = Mn; // knot order unchanged
    };
    auto swapAx = [&]()
    {
        gsMatrix<T> Mn(N,3);
        for (index_t j = 0; j < n0; ++j)
            for (index_t i = 0; i < n1; ++i)
                Mn.row(i + j*n1) = M.row(j + n0*i);
        M = Mn; std::swap(kvU,kvV); std::swap(n0,n1);
    };

    const bool switched = (geo.orientation() == -1);
    if (switched) swapAx();
    if (!switched)
    {
        if      (cIdx == 2) antiClock();
        else if (cIdx == 3) clock();
        else if (cIdx == 4) twice();
    }
    else
    {
        if      (cIdx == 2) clock();
        else if (cIdx == 3) antiClock();
        else if (cIdx == 4) twice();
    }

    gsMatrix<T> coefsRot = M.leftCols(2);
    geoRot = gsTensorBSpline<2,T>(kvU, kvV, coefsRot);
    perm.resize(N);
    for (index_t r = 0; r < N; ++r) perm[r] = static_cast<index_t>(std::lround(M(r,2)));
}

/// Is interface edge `e` incident to vertex `vi` (both endpoints'
/// corners belong to the vertex)?
template <class T>
bool edgeAtVertex(const EdgeInfo<T>& e, const VertexInfo<T>& vi)
{
    if (e.isBoundary) return false;
    const index_t ca = cornerIdxForPatch(vi, e.side1.patch);
    const index_t cb = cornerIdxForPatch(vi, e.side2.patch);
    if (ca < 0 || cb < 0) return false;
    const bool onA = (cornerAtSideEnd(e.side1,0)==ca) || (cornerAtSideEnd(e.side1,1)==ca);
    const bool onB = (cornerAtSideEnd(e.side2,0)==cb) || (cornerAtSideEnd(e.side2,1)==cb);
    return onA && onB;
}

/// Number of distinct interfaces meeting at vertex `vi`.
template <class T>
index_t vtxInterfaceCount(const VertexInfo<T>& vi,
                          const std::vector<EdgeInfo<T> >& edges)
{
    index_t c = 0;
    for (const EdgeInfo<T>& e : edges) if (edgeAtVertex(e, vi)) ++c;
    return c;
}

/// Find the vertex whose corner-list contains (patch, cornerIdx).
template <class T>
const VertexInfo<T>* vertexAtCorner(index_t patch, index_t cIdx,
                                    const std::vector<VertexInfo<T> >& verts)
{
    for (const VertexInfo<T>& v : verts)
        for (const patchCorner& c : v.corners)
            if (c.patch == patch && c.m_index == cIdx) return &v;
    return nullptr;
}

/// Build the 6 AS-G1 vertex shape functions for one interior vertex by
/// the (ported) gsUnstructuredSplines closed-form construction with
/// exact geometry-derived linear gluing data and interpolation at
/// anchors.  Output: shared columns written into every incident patch.
template <class T>
VertexFuncs<T> buildVertexFuncsKST(
    const gsMultiPatch<T>& mp,
    const std::vector<EdgeInfo<T> >& edges,
    const VertexInfo<T>& vi,
    const AsG1MpOptions<T>& opts)
{
    GISMO_UNUSED(edges);
    VertexFuncs<T> vf;
    vf.vertexId = vi.id;

    const index_t nInc = static_cast<index_t>(vi.corners.size());

    // --- reparametrize each incident patch to standard form ---
    std::vector<gsTensorBSpline<2,T> > geoRot(nInc);
    std::vector<std::vector<index_t> > perm(nInc);
    std::vector<index_t> patchOf(nInc);
    for (index_t k = 0; k < nInc; ++k)
    {
        patchOf[k] = vi.corners[k].patch;
        const gsTensorBSpline<2,T>& geo =
            dynamic_cast<const gsTensorBSpline<2,T>&>(mp.patch(patchOf[k]));
        reparamCorner<T>(geo, vi.corners[k].m_index, geoRot[k], perm[k]);
    }

    // --- sigma (uniform scaling) ---
    T p = 0, h_geo = 1;
    for (index_t k = 0; k < nInc; ++k)
    {
        const gsTensorBSplineBasis<2,T>& b = geoRot[k].basis();
        p = std::max<T>(p, std::max(b.degree(0), b.degree(1)));
    }
    for (index_t k = 0; k < nInc; ++k)
    {
        const gsTensorBSplineBasis<2,T>& b = geoRot[k].basis();
        for (index_t j = 0; j < 2; ++j)
            h_geo = std::min<T>(h_geo, b.knots(j).at(static_cast<size_t>(p)+1));
    }
    T sigmaRaw = 0;
    {
        gsMatrix<T> zero; zero.setZero(2,1);
        for (index_t k = 0; k < nInc; ++k)
        {
            gsMatrix<T> d = geoRot[k].deriv(zero);
            sigmaRaw += d.template lpNorm<gsEigen::Infinity>();
        }
    }
    const T sigmaScale = sigmaRaw * h_geo / (T(nInc) * p);
    const T sigma = (sigmaScale != T(0)) ? T(1)/sigmaScale : T(1);

    gsMatrix<T> Phi(6,6); Phi.setIdentity();
    Phi.col(1) *= sigma; Phi.col(2) *= sigma;
    Phi.col(3) *= sigma*sigma; Phi.col(4) *= sigma*sigma; Phi.col(5) *= sigma*sigma;

    // --- assign 6 global columns; build per-patch contributions ---
    std::map<std::pair<index_t,index_t>, index_t> localIdx;
    std::vector<std::map<index_t,T> > colAcc(6); // localdof -> value, per bfID

    gsMatrix<T> zero; zero.setZero(2,1);

    for (index_t k = 0; k < nInc; ++k)
    {
        gsTensorBSplineBasis<2,T>& tbRot =
            const_cast<gsTensorBSplineBasis<2,T>&>(geoRot[k].basis());

        // Determine kindOfEdge: is the South (u-edge, dir0) / West
        // (v-edge, dir1) of the reparametrized patch an interface?  Match
        // the rotated edge midpoints to the original corner's two sides.
        std::vector<bool> kindOfEdge(2, true);
        {
            const index_t pk = patchOf[k];
            std::vector<patchSide> csides;
            patchCorner(pk, vi.corners[k].m_index).getContainingSides(2, csides);
            const T Lu = dynamic_cast<const gsBSplineBasis<T>&>(tbRot.component(0)).support()(0,1);
            const T Lv = dynamic_cast<const gsBSplineBasis<T>&>(tbRot.component(1)).support()(0,1);
            gsMatrix<T> uvS(2,1); uvS << Lu/2, T(0);
            gsMatrix<T> uvW(2,1); uvW << T(0), Lv/2;
            gsMatrix<T> Ps = geoRot[k].eval(uvS);
            gsMatrix<T> Pw = geoRot[k].eval(uvW);
            for (int dir = 0; dir < 2; ++dir)
            {
                const gsMatrix<T>& Ptarget = (dir==0) ? Ps : Pw;
                T best = std::numeric_limits<T>::max(); bool bestI = true;
                for (const patchSide& s : csides)
                {
                    gsMatrix<T> mp_(2,1);
                    mp_(s.direction(),0) = s.parameter() ? T(1) : T(0);
                    mp_(1-s.direction(),0) = T(0.5);
                    gsMatrix<T> Pm = mp.patch(pk).eval(mp_);
                    const T d = (Pm - Ptarget).norm();
                    if (d < best) { best = d; bestI = mp.isInterface(s); }
                }
                kindOfEdge[dir] = bestI;
            }
        }

        // exact linear alpha/beta from the (bilinear) reparametrized geo
        std::vector<gsBSpline<T> > alpha(2), beta(2);
        std::vector<gsBSplineBasis<T> > bplus(2), bminus(2);
        for (index_t i = 0; i < 2; ++i)
        {
            const gsBSplineBasis<T>& comp =
                dynamic_cast<const gsBSplineBasis<T>&>(tbRot.component(i));
            const T L = comp.support()(0,1);
            T av[2], bv[2];
            const T ss[2] = { T(0), L };
            for (int q = 0; q < 2; ++q)
            {
                gsMatrix<T> uv(2,1); uv.setZero();
                uv(i,0) = ss[q];
                gsMatrix<T> J = geoRot[k].jacobian(uv);
                av[q] = J.determinant();
                const T d0n = J.col(i).squaredNorm();
                bv[q] = -(J.col(1).dot(J.col(0))) / d0n;
            }
            gsKnotVector<T> kvL(T(0), L, 0, 2);
            gsMatrix<T> ac(2,1); ac(0,0)=av[0]; ac(1,0)=av[1];
            gsMatrix<T> bc(2,1); bc(0,0)=bv[0]; bc(1,0)=bv[1];
            alpha[i] = gsBSpline<T>(kvL, ac);
            beta[i]  = gsBSpline<T>(kvL, bc);

            gsBSplineBasis<T> sm = comp; sm.elevateContinuity(1);
            gsBSplineBasis<T> ld = comp; ld.degreeReduce(1);
            bplus[i]  = sm;
            bminus[i] = ld;
        }

        gsMatrix<T> anchors = tbRot.anchors();
        for (index_t bfID = 0; bfID < 6; ++bfID)
        {
            gsVertexBasisMP<T> vb(geoRot[k], tbRot, alpha, beta, bplus, bminus,
                                  Phi, kindOfEdge, bfID);
            gsMatrix<T> vals = vb.eval(anchors);                 // 1 x nBasis
            typename gsGeometry<T>::uPtr g = tbRot.interpolateAtAnchors(give(vals));
            gsMatrix<T> rotCoefs = g->coefs();                  // nBasis x 1

            for (index_t r = 0; r < rotCoefs.rows(); ++r)
            {
                const T c = rotCoefs(r,0);
                if (std::abs(c) < T(1e-13)) continue;
                const index_t origDof = perm[k][r];
                std::pair<index_t,index_t> key(patchOf[k], origDof);
                if (localIdx.find(key) == localIdx.end())
                {
                    const index_t li = static_cast<index_t>(vf.localDofs.size());
                    localIdx[key] = li;
                    vf.localDofs.push_back(key);
                }
                colAcc[bfID][localIdx[key]] += c;
            }
        }
    }

    const index_t nLocal = static_cast<index_t>(vf.localDofs.size());
    vf.W = gsMatrix<T>::Zero(nLocal, 6);
    for (index_t bfID = 0; bfID < 6; ++bfID)
        for (const auto& kv : colAcc[bfID])
            vf.W(kv.first, bfID) = kv.second;

    if (opts.verbose)
        gsInfo << "  vertex " << vi.id << " (valence " << vi.valence
               << "): 6 KST functions, " << nLocal << " local DOFs\n";
    return vf;
}

// ====================================================================
// Build
// ====================================================================

template <class T>
AsG1Multipatch<T> buildAsG1Multipatch(
    const gsMultiPatch<T>& mpIn,
    const AsG1MpOptions<T>& opts = AsG1MpOptions<T>())
{
    AsG1Multipatch<T> out;
    out.mp = mpIn;
    gsMultiPatch<T>& mp = out.mp;
    mp.computeTopology();

    // ---- Degree elevation ----
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < opts.minDegree)
        mp.degreeElevate(opts.minDegree - inputDeg);

    // ---- Refinement (knot multiplicity = deg-1, C^1) ----
    if (opts.refinements > 0)
    {
        const short_t deg = mp.patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t i = 0; i < opts.refinements; ++i)
            mp.uniformRefine(1, mult);
    }

    const index_t nP = static_cast<index_t>(mp.nPatches());

    // ---- Topology ----
    out.edges = enumerateEdges(mp);
    out.verts = enumerateVertices(mp);

    // ---- Gluing data ----
    const index_t nFail = computeAllGluingData<T>(
        mp, out.edges, opts.gluingSource, opts.eps,
        opts.numGaussPerSpan, opts.verbose);
    GISMO_ENSURE(nFail == 0,
                 "buildAsG1Multipatch: gluing data failed on "
                     << nFail << " interface(s).");

    // ---- Per-patch interior DOFs ----
    out.interiorDofs.resize(nP);
    out.interiorOff.assign(nP + 1, 0);
    for (index_t p = 0; p < nP; ++p)
    {
        const gsTensorBSplineBasis<2,T>& tb = asTB(mp, p);
        std::set<index_t> layer = interfaceLayerDofs(mp, out.edges, p);
        std::vector<index_t> interior;
        interior.reserve(tb.size() - layer.size());
        for (index_t i = 0; i < tb.size(); ++i)
            if (layer.find(i) == layer.end())
                interior.push_back(i);
        out.interiorDofs[p] = interior;
        out.interiorOff[p + 1] =
            out.interiorOff[p] + static_cast<index_t>(interior.size());
    }
    const index_t totInterior = out.interiorOff[nP];
    out.edgeBlockOff = totInterior;

    // ---- Per-edge embeddings + global column allocation ----
    index_t gcol = out.edgeBlockOff;
    for (const EdgeInfo<T>& e : out.edges)
    {
        if (e.isBoundary) continue;

        EdgeEmbedding<T> em;
        em.edgeId  = e.id;
        em.flipped = e.flipped;

        const gsTensorBSplineBasis<2,T>& tb1 = asTB(mp, e.side1.patch);
        const gsTensorBSplineBasis<2,T>& tb2 = asTB(mp, e.side2.patch);

        const T a1_0 = e.gluingData(0), a1_1 = e.gluingData(1);
        const T b1_0 = e.gluingData(2), b1_1 = e.gluingData(3);
        const T a2_0 = e.gluingData(4), a2_1 = e.gluingData(5);
        const T b2_0 = e.gluingData(6), b2_1 = e.gluingData(7);
        const T tSign = e.flipped ? T(-1) : T(1);

        em.E1 = asg1v4::createGluingDataArgyrisBasis<T>(
            tb1, e.side1.side(), a1_0, a1_1, b1_0, b1_1, T(1e-12), tSign);
        em.E2 = asg1v4::createGluingDataArgyrisBasis<T>(
            tb2, e.side2.side(), a2_0, a2_1, b2_0, b2_1, T(1e-12), tSign);

        // Shared trace / deriv counts from the side bases.
        gsBSplineBasis<T> sb1 = *tb1.boundaryBasis(e.side1.side());
        gsBSplineBasis<T> smB1 = sb1; smB1.elevateContinuity(1);
        gsBSplineBasis<T> ldB1 = sb1; ldB1.degreeReduce(1);
        gsBSplineBasis<T> sb2 = *tb2.boundaryBasis(e.side2.side());
        gsBSplineBasis<T> smB2 = sb2; smB2.elevateContinuity(1);
        gsBSplineBasis<T> ldB2 = sb2; ldB2.degreeReduce(1);

        GISMO_ENSURE(smB1.size() == smB2.size(),
                     "Trace basis sizes differ across edge " << e.id);
        GISMO_ENSURE(ldB1.size() == ldB2.size(),
                     "Deriv basis sizes differ across edge " << e.id);

        em.nSm = smB1.size();
        em.nLD = ldB1.size();

        // E_i column layout: [interior | deriv(nLD) | trace(nSm) ]
        em.nInt1 = em.E1.cols() - em.nLD - em.nSm;
        em.nInt2 = em.E2.cols() - em.nLD - em.nSm;

        // Phase D: drop near-vertex cols at an end whose vertex is an
        // *extraordinary* point of the interface graph, i.e. where >= 2
        // interfaces meet (interior vertex, or a boundary vertex where two
        // interfaces meet).  At a one-interface end (regular boundary
        // corner) the columns are kept.
        auto endNeedsTrunc = [&](int end)->bool
        {
            if (!opts.truncate) return false;
            const index_t cIdx = cornerAtSideEnd(e.side1, end);
            const VertexInfo<T>* v =
                vertexAtCorner<T>(e.side1.patch, cIdx, out.verts);
            return v && vtxInterfaceCount<T>(*v, out.edges) >= 2;
        };
        const bool loInt = endNeedsTrunc(0);
        const bool hiInt = endNeedsTrunc(1);
        em.traceLo = loInt ? ASG1_TRACE_DROP : 0;
        em.traceHi = hiInt ? ASG1_TRACE_DROP : 0;
        em.derivLo = loInt ? ASG1_DERIV_DROP : 0;
        em.derivHi = hiInt ? ASG1_DERIV_DROP : 0;
        // Clamp so at least one trace/deriv column survives on very coarse
        // interfaces (truncation needs enough refinement, typically r>=2).
        if (em.traceLo + em.traceHi > em.nSm - 1)
        {
            const index_t avail = std::max<index_t>(em.nSm - 1, 0);
            em.traceLo = std::min(em.traceLo, avail);
            em.traceHi = std::min(em.traceHi, avail - em.traceLo);
        }
        if (em.derivLo + em.derivHi > em.nLD - 1)
        {
            const index_t avail = std::max<index_t>(em.nLD - 1, 0);
            em.derivLo = std::min(em.derivLo, avail);
            em.derivHi = std::min(em.derivHi, avail - em.derivLo);
        }
        em.nSmKept = em.nSm - em.traceLo - em.traceHi;
        em.nLDKept = em.nLD - em.derivLo - em.derivHi;
        GISMO_ENSURE(em.nSmKept >= 0 && em.nLDKept >= 0,
                     "Edge " << e.id << " too coarse for truncation.");

        // Allocate global columns: kept trace block then kept deriv block.
        em.gOffTrace = gcol;            gcol += em.nSmKept;
        em.gOffDeriv = gcol;            gcol += em.nLDKept;

        out.edgeEmb.push_back(em);
    }
    out.nGlobal = gcol;

    // ---- Allocate G[p] ----
    out.G.resize(nP);
    for (index_t p = 0; p < nP; ++p)
        out.G[p] = gsSparseMatrix<T>(asTB(mp, p).size(), out.nGlobal);

    // ---- Interior identity blocks ----
    for (index_t p = 0; p < nP; ++p)
    {
        const index_t off = out.interiorOff[p];
        const std::vector<index_t>& interior = out.interiorDofs[p];
        for (size_t j = 0; j < interior.size(); ++j)
            out.G[p].insert(interior[j], off + static_cast<index_t>(j)) = T(1);
    }

    // ---- Edge coupling blocks ----
    for (const EdgeEmbedding<T>& em : out.edgeEmb)
    {
        const EdgeInfo<T>& e = out.edges[em.edgeId];
        const index_t p1 = e.side1.patch, p2 = e.side2.patch;

        // -- Patch 1 (aligned) --  kept deriv cols -> global deriv block
        for (index_t k = 0; k < em.nLDKept; ++k)
        {
            const index_t jFull = em.derivLo + k;
            const index_t e1col = em.nInt1 + jFull;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E1, e1col); it; ++it)
                out.G[p1].insert(it.row(), em.gOffDeriv + k) = it.value();
        }
        // kept trace cols -> global trace block
        for (index_t k = 0; k < em.nSmKept; ++k)
        {
            const index_t jFull = em.traceLo + k;
            const index_t e1col = em.nInt1 + em.nLD + jFull;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E1, e1col); it; ++it)
                out.G[p1].insert(it.row(), em.gOffTrace + k) = it.value();
        }

        // -- Patch 2 (flipped-aware signs/index, as in v4) --
        const T l2Sign = em.flipped ? T(1) : T(-1);
        for (index_t k = 0; k < em.nLDKept; ++k)
        {
            const index_t jFull = em.derivLo + k;
            const index_t j2 = em.flipped ? (em.nLD - 1 - jFull) : jFull;
            const index_t e2col = em.nInt2 + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E2, e2col); it; ++it)
                out.G[p2].insert(it.row(), em.gOffDeriv + k) = l2Sign * it.value();
        }
        for (index_t k = 0; k < em.nSmKept; ++k)
        {
            const index_t jFull = em.traceLo + k;
            const index_t j2 = em.flipped ? (em.nSm - 1 - jFull) : jFull;
            const index_t e2col = em.nInt2 + em.nLD + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(em.E2, e2col); it; ++it)
                out.G[p2].insert(it.row(), em.gOffTrace + k) = it.value();
        }
    }

    for (index_t p = 0; p < nP; ++p)
        out.G[p].makeCompressed();

    out.vertexBlockOff = out.nGlobal;

    // ---- Phase E: interior vertex space ----
    if (opts.vertex)
    {
        index_t addCols = 0;
        for (const VertexInfo<T>& vi : out.verts)
        {
            // Build a vertex block wherever >= 2 interfaces meet
            // (interior vertex or 2-interface boundary vertex).
            if (vtxInterfaceCount<T>(vi, out.edges) < 2) continue;
            VertexFuncs<T> vf = buildVertexFuncsKST<T>(mp, out.edges, vi, opts);
            if (vf.W.cols() == 0) continue;
            vf.gColBase = out.nGlobal + addCols;
            addCols += vf.W.cols();
            out.vertexFuncs.push_back(vf);
        }

        if (addCols > 0)
        {
            const index_t newN = out.nGlobal + addCols;
            for (index_t p = 0; p < nP; ++p)
            {
                gsSparseEntries<T> ent;
                for (index_t c = 0; c < out.G[p].outerSize(); ++c)
                    for (typename gsSparseMatrix<T>::InnerIterator it(out.G[p], c); it; ++it)
                        ent.add(it.row(), it.col(), it.value());

                for (const VertexFuncs<T>& vf : out.vertexFuncs)
                    for (index_t li = 0; li < static_cast<index_t>(vf.localDofs.size()); ++li)
                    {
                        if (vf.localDofs[li].first != p) continue;
                        const index_t dof = vf.localDofs[li].second;
                        for (index_t q = 0; q < vf.W.cols(); ++q)
                        {
                            const T w = vf.W(li, q);
                            if (std::abs(w) > T(1e-13))
                                ent.add(dof, vf.gColBase + q, w);
                        }
                    }

                gsSparseMatrix<T> Gp(asTB(mp, p).size(), newN);
                Gp.setFrom(ent);
                Gp.makeCompressed();
                out.G[p] = Gp;
            }
            out.nGlobal = newN;
        }
    }

    return out;
}

// ====================================================================
// Smoothness check (per interface, physical gradient)
// ====================================================================

template <class T>
struct G1ReportMP
{
    T maxValErr  = T(0);
    T maxGradErr = T(0);
    bool pass    = false;
    std::vector<T> perEdgeGradErr;
};

/// Verify AS-G1 of every global basis function across every interior
/// edge, in physical (gradient) space, at `nCheck` sample points.
template <class T>
G1ReportMP<T> smoothnessCheckMP(const AsG1Multipatch<T>& B, index_t nCheck = 21)
{
    G1ReportMP<T> R;
    const gsMultiPatch<T>& mp = B.mp;

    for (const EdgeInfo<T>& e : B.edges)
    {
        if (e.isBoundary) continue;
        const patchSide ps1 = e.side1, ps2 = e.side2;
        const short_t nd1 = ps1.direction(), td1 = 1 - nd1;
        const short_t nd2 = ps2.direction(), td2 = 1 - nd2;
        const bool par1 = ps1.parameter(), par2 = ps2.parameter();

        const gsTensorBSplineBasis<2,T>& tb1 = asTB(mp, ps1.patch);
        const gsTensorBSplineBasis<2,T>& tb2 = asTB(mp, ps2.patch);

        gsMatrix<T> sup1 = mp.patch(ps1.patch).support();
        gsMatrix<T> sup2 = mp.patch(ps2.patch).support();
        const T ifc1 = sup1(nd1, par1 ? 1 : 0);
        const T ifc2 = sup2(nd2, par2 ? 1 : 0);
        const T t1a = sup1(td1, 0), t1b = sup1(td1, 1);
        const T t2a = sup2(td2, 0), t2b = sup2(td2, 1);

        // Only global DOFs that are nonzero near this interface can break
        // G1 there; interior identity DOFs away from it vanish (value and
        // gradient) on the interface.  Collect the relevant columns by
        // scanning the two patches' first/second-layer rows.
        std::set<index_t> layerRows1, layerRows2;
        {
            gsMatrix<index_t> b0 = tb1.boundary(ps1.side()), b1 = tb1.boundaryOffset(ps1.side(),1);
            for (index_t i=0;i<b0.rows();++i) layerRows1.insert(b0(i,0));
            for (index_t i=0;i<b1.rows();++i) layerRows1.insert(b1(i,0));
            gsMatrix<index_t> c0 = tb2.boundary(ps2.side()), c1m = tb2.boundaryOffset(ps2.side(),1);
            for (index_t i=0;i<c0.rows();++i) layerRows2.insert(c0(i,0));
            for (index_t i=0;i<c1m.rows();++i) layerRows2.insert(c1m(i,0));
        }
        std::set<index_t> activeCols;
        for (index_t c=0;c<B.G[ps1.patch].outerSize();++c)
            for (typename gsSparseMatrix<T>::InnerIterator it(B.G[ps1.patch],c); it; ++it)
                if (layerRows1.count(it.row())) { activeCols.insert(it.col()); break; }
        for (index_t c=0;c<B.G[ps2.patch].outerSize();++c)
            for (typename gsSparseMatrix<T>::InnerIterator it(B.G[ps2.patch],c); it; ++it)
                if (layerRows2.count(it.row())) { activeCols.insert(it.col()); break; }

        T edgeMaxGrad = T(0);

        for (index_t idx : activeCols)
        {
            gsVector<T> g = gsVector<T>::Zero(B.nGlobal);
            g(idx) = T(1);
            gsVector<T> c1 = B.G[ps1.patch] * g;
            gsVector<T> c2 = B.G[ps2.patch] * g;
            auto f1 = tb1.makeGeometry(c1);
            auto f2 = tb2.makeGeometry(c2);

            for (index_t i = 0; i < nCheck; ++i)
            {
                const T s = T(i) / T(nCheck - 1);
                const T t1 = t1a + s * (t1b - t1a);
                const T s2 = e.flipped ? (T(1) - s) : s;
                const T t2 = t2a + s2 * (t2b - t2a);

                gsMatrix<T> pt1(2,1), pt2(2,1);
                pt1(nd1,0) = ifc1; pt1(td1,0) = t1;
                pt2(nd2,0) = ifc2; pt2(td2,0) = t2;

                R.maxValErr = std::max(R.maxValErr,
                    std::abs(f1->eval(pt1)(0,0) - f2->eval(pt2)(0,0)));

                gsMatrix<T> df1, df2, dG1, dG2;
                f1->deriv_into(pt1, df1);
                f2->deriv_into(pt2, df2);
                mp.patch(ps1.patch).deriv_into(pt1, dG1);
                mp.patch(ps2.patch).deriv_into(pt2, dG2);

                gsMatrix<T> J1(2,2), J2(2,2);
                J1(0,0)=dG1(0,0); J1(0,1)=dG1(1,0); J1(1,0)=dG1(2,0); J1(1,1)=dG1(3,0);
                J2(0,0)=dG2(0,0); J2(0,1)=dG2(1,0); J2(1,0)=dG2(2,0); J2(1,1)=dG2(3,0);

                gsVector<T> pg1(2), pg2(2);
                pg1(0)=df1(0,0); pg1(1)=df1(1,0);
                pg2(0)=df2(0,0); pg2(1)=df2(1,0);

                gsVector<T> phys1 = J1.inverse().transpose() * pg1;
                gsVector<T> phys2 = J2.inverse().transpose() * pg2;
                const T ge = (phys1 - phys2).norm();
                R.maxGradErr = std::max(R.maxGradErr, ge);
                edgeMaxGrad  = std::max(edgeMaxGrad, ge);
            }
        }
        R.perEdgeGradErr.push_back(edgeMaxGrad);
    }
    R.pass = (R.maxValErr < T(1e-8)) && (R.maxGradErr < T(1e-3));
    return R;
}

// ====================================================================
// Export to gsMappedBasis (for PDE / biharmonic assemblers)
// ====================================================================

/// Stack the per-patch matrices `G[p]` (rows = local patch DOFs) into a
/// single (totalLocalDofs x nGlobal) sparse weight-mapper matrix, with
/// local DOFs ordered patch-by-patch (the gsMultiBasis convention).
template <class T>
gsSparseMatrix<T> mappingMatrix(const AsG1Multipatch<T>& B)
{
    index_t totalLocal = 0;
    for (size_t p = 0; p < B.G.size(); ++p)
        totalLocal += B.G[p].rows();

    gsSparseEntries<T> ent;
    index_t off = 0;
    for (size_t p = 0; p < B.G.size(); ++p)
    {
        for (index_t c = 0; c < B.G[p].outerSize(); ++c)
            for (typename gsSparseMatrix<T>::InnerIterator it(B.G[p], c); it; ++it)
                ent.add(off + it.row(), it.col(), it.value());
        off += B.G[p].rows();
    }
    gsSparseMatrix<T> M(totalLocal, B.nGlobal);
    M.setFrom(ent);
    M.makeCompressed();
    return M;
}

/// Build a gsMappedBasis from the AS-G1 multi-patch construction.  The
/// resulting basis has `B.nGlobal` global functions and plugs directly
/// into gsMappedSpline / the biharmonic assemblers.
template <class T>
void toMappedBasis(const AsG1Multipatch<T>& B, gsMappedBasis<2,T>& mbasis)
{
    gsMultiBasis<T> mb(B.mp);
    gsSparseMatrix<T> M = mappingMatrix(B);
    mbasis.init(mb, M);
}

} // namespace asg1mp
} // namespace gismo
