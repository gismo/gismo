/** @file gsMomentRule.h

    @brief Solve-free nodal moment fitting: compresses the points of an
    arbitrary wrapped quadrature rule onto the tensor Gauss grid of the element.

    A cut-cell (or otherwise adaptive) rule may deliver very many quadrature
    points per element -- adaptive subdivision multiplies the leaf rule by the
    size of the tree. gsMomentRule replaces that point cloud by the *fixed*
    tensor Gauss grid of the element, carrying modified weights, so that the
    assembly cost per element no longer depends on how hard the integrand (or
    the interface) was to resolve.

    Algorithm (per element, in the physical frame of [lower,upper]): let
    \f$ g_j \f$ be the \f$ \prod_d n_d \f$ tensor Gauss points of the element
    and \f$ L_j \f$ the nodal Lagrange basis of that grid. For an integrand f,
    \f[ \int f \approx \sum_j f(g_j)\, w_j, \qquad
        w_j = \sum_k \omega_k\, \alpha(x_k)\, L_j(x_k), \f]
    where \f$ (x_k,\omega_k) \f$ are the points and weights returned by the
    <em>wrapped</em> rule and \f$ \alpha \f$ is an optional multiplicative
    field. This is nodal moment fitting: the moment matrix
    \f$ A_{ij} = L_j(g_i) \f$ is the identity by the interpolation property of
    the nodal basis, so <b>no linear system is solved</b> (the "solve-free"
    variant). The Lagrange functions are evaluated in the numerically stable
    second (normalized) barycentric form.

    <b>The class carries no geometry knowledge</b>: no level set, no
    classification, no inside/outside/cut logic. It is a pure compressor around
    another rule. Elements that need no compression are passed through
    unchanged (see below), which makes uncut elements free without any
    classification at all.

    <b>Ownership.</b> The rule <b>takes ownership</b> of the wrapped rule
    (a gsQuadRule<T>::uPtr moved in). The gsQuadRule hierarchy has no clone(),
    so ownership transfer is the only way to hold a polymorphic rule without
    slicing it -- a by-value or by-const-reference copy would silently discard
    any overridden mapTo(). There is deliberately no constructor taking a
    const gsQuadRule<T>&.

    <b>Exactness (read this before choosing the order).</b> Moment fitting
    reproduces polynomials only up to tensor degree \f$ n-1 \f$ per direction
    (the interpolation degree of the nodal Lagrange basis), whereas Gauss with
    the same \f$ n \f$ is exact to degree \f$ 2n-1 \f$: moment fitting is
    therefore roughly twice as demanding in point count as Gauss for the same
    integrand. A mass-like integrand \f$ \phi_i\phi_j \f$ (tensor degree
    \f$ 2p \f$ for a degree-\f$ p \f$ basis) needs
    \f[ n \geq 2p+1 \quad\text{per direction,} \f]
    which is what exactnessOrder() returns. A stiffness integrand
    \f$ \nabla\phi_i\cdot\nabla\phi_j \f$ needs <b>the same</b> \f$ n \f$ in
    \f$ d \geq 2 \f$ despite differentiation lowering one factor: for a
    tensor-product basis \f$ \partial_x\phi_i\,\partial_x\phi_j =
    (N'_aN'_c)(x)\,(N_bN_d)(y) \f$ has degree \f$ 2p-2 \f$ in \f$ x \f$ but
    still \f$ 2p \f$ in \f$ y \f$. Only in 1D does the reduced requirement
    \f$ n \geq 2p-1 \f$ apply. Under-integrating here is not a small accuracy
    loss: measured on a mass-matrix probe with p = 2, the exactness order
    n = 5 stays SPD at every refinement, while n = 3 is <b>indefinite at every
    refinement</b> (a genuine negative diagonal entry), and a Poisson+Nitsche
    solve with the under-integrated order returns NaN on the coarsest mesh. At
    the exactness order the same solve reproduces the baseline accuracy with
    two orders of magnitude fewer cut-cell quadrature points.

    <b>Why the order is NOT read off the wrapped rule.</b> Two independent
    reasons:
    1. <em>Mechanically unavailable.</em> gsQuadRule<T>::numNodes() and dim()
       report the <em>base</em> members m_weights/m_nodes. Rules that compute
       their points inside mapTo() (the cut-cell rules this class is meant to
       wrap) never populate those members, so both queries return 0 for
       precisely the rules of interest; gsPatchRule additionally shadows dim()
       non-virtually and declares its own members, so it reports 0 too through
       a gsQuadRule<T>*.
    2. <em>Wrong even when available.</em> A Gauss rule's \f$ n \f$ is sized
       for exactness \f$ 2n-1 \f$; moment fitting reaches only \f$ n-1 \f$ and
       therefore needs roughly <b>double</b>. Copying the wrapped rule's order
       reproduces the under-integration failure above by construction.
    The order is therefore always supplied by the caller, ideally through
    exactnessOrder().

    <b>Negative weights are expected and harmless</b> at an adequate order.
    They occur for small cut fractions, or for interfaces that the Gauss grid
    resolves poorly. They do <b>not</b> by themselves imply an indefinite
    assembled operator: in the measurement quoted above the n = 5 configuration
    is SPD at every refinement while carrying hundreds of negative weights. Any
    caller relying on the exact positivity of individual matrix diagonals
    should nevertheless know that sliver degrees of freedom can push a diagonal
    marginally negative (observed magnitudes orders below eps*lambda_max, i.e.
    numerically harmless, but systematic in occurrence). The count and the most
    negative weight are recorded in Stats. Follow-up: non-negative moment
    fitting, arXiv:2604.15921.

    <b>Total mass is preserved identically</b>: \f$ \sum_j w_j =
    \sum_k \omega_k \sum_j L_j(x_k) = \sum_k \omega_k \f$ by the partition of
    unity of the Lagrange basis. Reproducing the area/volume of a cut cell is
    therefore <b>not</b> evidence that the compression is correct -- accuracy
    must be checked against non-constant integrands.

    <b>alpha is the integrand field, not the FCM blend.</b> The optional
    \f$ \alpha \f$ of this class is a multiplicative scalar field <em>inside</em>
    the integrand (a CT material field, a Young's modulus, an indicator). It
    multiplies \f$ \omega_k \f$ at the <em>underlying</em> points, because it is
    integrated by the wrapped rule. It is <b>not</b> the fictitious-domain blend
    \f$ w \leftarrow (1-a)w + a\,w^G \f$ of the finite cell method, which
    re-adds material that the wrapped rule dropped and needs the analytic
    full-cell Gauss weight \f$ w^G \f$ -- something \f$ \alpha \f$ evaluated at
    the underlying points cannot recover. That blend is not implemented here.

    <b>Deliberate capability gap: the FCM blend has no replacement.</b> The
    deleted gsAlgoimMomentFittingRule carried the blend as a scalar
    \c fictitiousWeight in [0,1] (default 0). It was dropped when moment fitting
    moved into the core, because a pure compressor never sees the geometry and
    so cannot form \f$ w^G \f$ -- reconstructing it would mean handing the rule
    a level set again, which is exactly the coupling this class exists to
    remove. The cost is real and is paid by the immersed drivers: their
    \c --alpha flags are gated by
    \c GISMO_ENSURE(alpha==0) rather than silently ignoring the value, and the
    former \c MomentFitting_Alpha unit test is gone. Anyone needing the blend
    should add it as an explicit scalar member here (blending the compressed
    weights against an analytic full-cell Gauss rule built from the element
    box, which this class does know), NOT by widening \f$ \alpha \f$ -- the two
    have different meanings and must stay separately named.

    The rule keeps mutable state (the lazily built 1D caches, an evaluation
    scratch buffer and the statistics), so it is <b>not thread-safe</b>: use one
    rule object per thread.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

// NOTE: gsQuadrature.h includes THIS header; never include it (or anything
// pulling it in) from here.
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsFunction.h>
#include <gsUtils/gsCombinatorics.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace gismo
{

/**
    \brief Solve-free nodal moment fitting rule: compresses the points of a
    wrapped (owned) quadrature rule onto the tensor Gauss grid of the element.

    See the file header for the algorithm, the exactness requirement
    \f$ n \geq 2p+1 \f$, and the semantics of \a alpha.

    \ingroup Assembler
*/
template<class T>
class gsMomentRule : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsMomentRule> uPtr;

    /// \brief Element and quadrature-point tallies, accumulated over all
    /// mapTo() calls since construction or the last resetStats().
    ///
    /// nElements counts every mapTo() call. The points actually emitted are
    /// <b>nOutputQPs + nPassThroughQPs</b>: an element whose wrapped rule
    /// returns nothing emits nothing and contributes to nElements only.
    /// nOutputQPs, nNegativeWeights and minWeight describe the *compressed*
    /// elements only; nPassThroughQPs the ones emitted unchanged.
    /// Comparing nOutputQPs alone against nUnderlyingQPs therefore overstates
    /// the compression, since nUnderlyingQPs covers every visited element,
    /// pass-through ones included.
    struct Stats
    {
        index_t nElements = 0;          ///< mapTo() calls
        index_t nUnderlyingQPs = 0;     ///< points returned by the wrapped rule
        index_t nOutputQPs = 0;         ///< points emitted by compressed elements
        index_t nPassThroughQPs = 0;    ///< points emitted unchanged
        index_t nPassThroughElements = 0;
        index_t nNegativeWeights = 0;   ///< compressed elements only
        T minWeight = std::numeric_limits<T>::max();

        Stats & operator+=(const Stats & other)
        {
            nElements += other.nElements;
            nUnderlyingQPs += other.nUnderlyingQPs;
            nOutputQPs += other.nOutputQPs;
            nPassThroughQPs += other.nPassThroughQPs;
            nPassThroughElements += other.nPassThroughElements;
            nNegativeWeights += other.nNegativeWeights;
            minWeight = std::min(minWeight, other.minWeight);
            return *this;
        }
    };

    /// \brief Anisotropic constructor.
    ///
    /// Takes OWNERSHIP of \a rule (the gsQuadRule hierarchy has no clone(), so
    /// a polymorphic rule cannot be copied without slicing).
    ///
    /// \param rule      the wrapped rule, moved in
    /// \param numNodes  output points per direction, one entry per direction
    ///                  (see exactnessOrder())
    /// \param alpha     optional multiplicative field inside the integrand,
    ///                  NOT owned: it must outlive this rule. It is the
    ///                  integrand's field, not the FCM fictitious blend.
    gsMomentRule(typename gsQuadRule<T>::uPtr rule,
                 gsVector<index_t> numNodes,
                 const gsFunction<T> * alpha = nullptr)
    :
    gsQuadRule<T>(),
    m_rule(give(rule)),
    m_numNodes(give(numNodes)),
    m_alpha(alpha),
    m_isotropic(0),
    m_initialized(false),
    m_hitTol(0)
    {
        GISMO_ENSURE(nullptr != m_rule.get(),
                     "gsMomentRule: the wrapped rule must not be null.");
        GISMO_ENSURE(m_numNodes.rows() > 0,
                     "gsMomentRule: numNodes must have at least one entry.");
        GISMO_ENSURE((m_numNodes.array() >= 2).all(),
                     "gsMomentRule: at least 2 output points per direction are "
                     "required; a one-point rule reproduces only constants. Got "
                     << m_numNodes.transpose() << ".");
        _checkAlpha();
    }

    /// \brief Isotropic convenience constructor: \a numNodes points in every
    /// direction, the dimension being taken from the first mapTo() call.
    ///
    /// Takes OWNERSHIP of \a rule; \a alpha is not owned (see above).
    gsMomentRule(typename gsQuadRule<T>::uPtr rule,
                 index_t numNodes,
                 const gsFunction<T> * alpha = nullptr)
    :
    gsQuadRule<T>(),
    m_rule(give(rule)),
    m_alpha(alpha),
    m_isotropic(numNodes),
    m_initialized(false),
    m_hitTol(0)
    {
        GISMO_ENSURE(nullptr != m_rule.get(),
                     "gsMomentRule: the wrapped rule must not be null.");
        GISMO_ENSURE(m_isotropic >= 2,
                     "gsMomentRule: at least 2 output points per direction are "
                     "required; a one-point rule reproduces only constants. Got "
                     << m_isotropic << ".");
        _checkAlpha();
    }

    /// Make function returning a smart pointer
    static uPtr make(typename gsQuadRule<T>::uPtr rule,
                     gsVector<index_t> numNodes,
                     const gsFunction<T> * alpha = nullptr)
    { return uPtr(new gsMomentRule(give(rule), give(numNodes), alpha)); }

    /// Make function returning a smart pointer
    static uPtr make(typename gsQuadRule<T>::uPtr rule,
                     index_t numNodes,
                     const gsFunction<T> * alpha = nullptr)
    { return uPtr(new gsMomentRule(give(rule), numNodes, alpha)); }

    /// \brief Points per direction making the compressed rule exact for an
    /// integrand of tensor degree 2*\a degree (mass-like, and stiffness-like
    /// too for \a dim >= 2): n = 2*degree+1 in each of \a dim directions.
    ///
    /// See the file header for why the order can neither be read off nor
    /// copied from the wrapped rule.
    static gsVector<index_t> exactnessOrder(short_t degree, short_t dim)
    {
        return gsVector<index_t>::Constant(dim,
                   2 * static_cast<index_t>(degree) + 1);
    }

    using gsQuadRule<T>::mapTo; // unhide the other overloads

    /// \brief Quadrature of the element [\a lower, \a upper], compressed onto
    /// its tensor Gauss grid (see the file header).
    ///
    /// Complexity per compressed element: O(K * prod_d n_d) products, with K
    /// the number of points of the wrapped rule -- no linear solve.
    ///
    /// If the wrapped rule returns no points, none are emitted (nodes has zero
    /// columns): gsExprAssembler skips such elements entirely.
    void mapTo(const gsVector<T> & lower,
               const gsVector<T> & upper,
               gsMatrix<T> & nodes,
               gsVector<T> & weights) const override
    {
        GISMO_ASSERT(lower.rows() == upper.rows(),
                     "gsMomentRule: lower and upper bounds must have the same "
                     "dimension.");
        const short_t d = static_cast<short_t>(lower.rows());
        _ensureInit(d);

        gsMatrix<T> uPts;
        gsVector<T> uWts;
        m_rule->mapTo(lower, upper, uPts, uWts);

        ++m_stats.nElements;
        m_stats.nUnderlyingQPs += uPts.cols();

        // Zero-column contract: an element without points is skipped by the
        // assembler. Emitting a grid of zero weights instead would be wrong
        // (and needlessly expensive).
        if (0 == uPts.cols())
        {
            nodes.resize(d, 0);
            weights.resize(0);
            return;
        }

        // Compressing few points onto more points only adds interpolation
        // error, so pass them through. Disabled whenever an alpha field is
        // set: alpha must be applied even when there are few points.
        if (uPts.cols() <= m_numNodes.prod() && nullptr == m_alpha)
        {
            nodes = give(uPts);
            weights = give(uWts);
            ++m_stats.nPassThroughElements;
            m_stats.nPassThroughQPs += nodes.cols();
            return;
        }

        // The tensor Gauss weights of the (uncut) element are not needed here,
        // only the grid points; they are discarded.
        gsVector<T> gWts;
        _outputGrid(lower, upper, nodes, gWts);

        if (nullptr != m_alpha)
            m_alpha->eval_into(uPts, m_alphaVals);

        weights.setZero(m_numNodes.prod());
        for (index_t k = 0; k != uPts.cols(); ++k)
            _addLagrange(lower, upper, uPts, k,
                         nullptr != m_alpha ? uWts[k] * m_alphaVals(0,k) : uWts[k],
                         weights);

        m_stats.nOutputQPs += weights.size();
        m_stats.minWeight = std::min(m_stats.minWeight, weights.minCoeff());
        m_stats.nNegativeWeights +=
            static_cast<index_t>((weights.array() < static_cast<T>(0)).count());
    }

    void mapTo(T, T, gsMatrix<T> &, gsVector<T> &) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// The wrapped rule.
    const gsQuadRule<T> & underlying() const { return *m_rule; }

    /// \brief Output points per direction. Empty until the first mapTo() call
    /// when the isotropic constructor was used (the dimension is only known
    /// then).
    const gsVector<index_t> & numNodesPerDirection() const { return m_numNodes; }

    /// Statistics accumulated over all mapTo() calls since construction or the
    /// last resetStats().
    const Stats & stats() const { return m_stats; }
    void resetStats() const { m_stats = Stats(); }

private:

    /// The optional alpha field must be scalar-valued.
    void _checkAlpha() const
    {
        if (nullptr != m_alpha)
            GISMO_ENSURE(1 == m_alpha->targetDim(),
                         "gsMomentRule: alpha must be a scalar field, got "
                         "targetDim() = " << m_alpha->targetDim() << ".");
    }

    /// \brief Builds, on the first call, the per-direction 1D Gauss data and
    /// the barycentric weights b_j = 1 / prod_{m!=j} (t_j - t_m).
    ///
    /// Lazy because the dimension is only known from mapTo()'s arguments.
    void _ensureInit(short_t d) const
    {
        if (m_initialized)
        {
            GISMO_ASSERT(m_numNodes.rows() == d,
                         "gsMomentRule: dimension " << d << " does not match the "
                         "dimension " << m_numNodes.rows() << " of the first call.");
            return;
        }

        if (0 == m_numNodes.rows())     // isotropic constructor
            m_numNodes = gsVector<index_t>::Constant(d, m_isotropic);
        else
            GISMO_ENSURE(m_numNodes.rows() == d,
                         "gsMomentRule: numNodes has " << m_numNodes.rows()
                         << " entries but the element is " << d << "-dimensional.");

        m_nodes1d.resize(d);
        m_gweights1d.resize(d);
        m_bary1d.resize(d);
        m_lagBuffer.resize(d);

        for (short_t i = 0; i != d; ++i)
        {
            const index_t n = m_numNodes[i];
            gsGaussRule<T> gauss(n);
            const gsVector<T> t = gauss.referenceNodes().row(0).transpose();

            gsVector<T> b(n);
            for (index_t j = 0; j != n; ++j)
            {
                T prod = 1;
                for (index_t m = 0; m != n; ++m)
                    if (m != j) prod *= (t[j] - t[m]);
                b[j] = static_cast<T>(1) / prod;
            }

            m_nodes1d[i] = t;
            m_gweights1d[i] = gauss.referenceWeights();
            m_bary1d[i] = give(b);
        }

        // Below this distance to a node the barycentric quotient blows up and
        // the Lagrange vector is the indicator of that node.
        m_hitTol = math::max(static_cast<T>(1e-13),
                             static_cast<T>(10) * std::numeric_limits<T>::epsilon());

        m_initialized = true;
    }

    /// Tensor Gauss rule of the (uncut) element from the cached 1D data, in
    /// the lexicographic order (direction 0 fastest) shared with
    /// _addLagrange(). Built here instead of via gsGaussRule::mapTo() so that
    /// node order and Lagrange index agree by construction.
    void _outputGrid(const gsVector<T> & lower,
                     const gsVector<T> & upper,
                     gsMatrix<T> & nodes,
                     gsVector<T> & weights) const
    {
        const short_t d = static_cast<short_t>(lower.rows());
        const index_t sz = m_numNodes.prod();
        nodes.resize(d, sz);
        weights.resize(sz);

        gsVector<index_t> cur = gsVector<index_t>::Zero(d);
        index_t flat = 0;
        do
        {
            T w = 1;
            for (short_t i = 0; i != d; ++i)
            {
                const T half = static_cast<T>(0.5) * (upper[i] - lower[i]);
                nodes(i,flat) = lower[i] + half * (m_nodes1d[i][cur[i]] + static_cast<T>(1));
                w *= half * m_gweights1d[i][cur[i]];
            }
            weights[flat++] = w;
        } while (nextLexicographic(cur, m_numNodes));
    }

    /// Second (normalized) barycentric form of the nodal Lagrange basis of
    /// direction \a dir at the reference coordinate \a u in [-1,1].
    void _lagrange1d(short_t dir, T u, gsVector<T> & l) const
    {
        const gsVector<T> & t = m_nodes1d[dir];
        const gsVector<T> & b = m_bary1d[dir];
        const index_t n = t.size();
        l.setZero(n);
        for (index_t j = 0; j != n; ++j)
            if (math::abs(u - t[j]) < m_hitTol) { l[j] = static_cast<T>(1); return; }

        T sum = 0;
        for (index_t j = 0; j != n; ++j)
        {
            l[j] = b[j] / (u - t[j]);
            sum += l[j];
        }
        l /= sum;   // scale invariant: the magnitude of b never matters
    }

    /// Adds \a weight * L_j(pts.col(k)) to \a result, for all output nodes j.
    void _addLagrange(const gsVector<T> & lower,
                      const gsVector<T> & upper,
                      const gsMatrix<T> & pts,
                      index_t k,
                      T weight,
                      gsVector<T> & result) const
    {
        const short_t d = static_cast<short_t>(lower.rows());
        for (short_t i = 0; i != d; ++i)
        {
            const T u = (static_cast<T>(2) * pts(i,k) - lower[i] - upper[i])
                      / (upper[i] - lower[i]);
            _lagrange1d(i, u, m_lagBuffer[i]);
        }

        gsVector<index_t> cur = gsVector<index_t>::Zero(d);
        index_t flat = 0;
        do
        {
            T v = weight;
            for (short_t i = 0; i != d; ++i)
                v *= m_lagBuffer[i][cur[i]];
            result[flat++] += v;
        } while (nextLexicographic(cur, m_numNodes));
    }

private:

    /// The wrapped rule, owned.
    typename gsQuadRule<T>::uPtr m_rule;

    /// Output points per direction; filled lazily when the isotropic
    /// constructor was used (hence mutable).
    mutable gsVector<index_t> m_numNodes;

    /// Optional multiplicative field inside the integrand; NOT owned.
    const gsFunction<T> * m_alpha;

    /// Points per direction of the isotropic constructor (0 otherwise).
    index_t m_isotropic;

    /// Output grid: per-direction 1D Gauss nodes, Gauss weights and
    /// barycentric weights on [-1,1], built on the first mapTo() call.
    mutable bool m_initialized;
    mutable std::vector<gsVector<T> > m_nodes1d;
    mutable std::vector<gsVector<T> > m_gweights1d;
    mutable std::vector<gsVector<T> > m_bary1d;
    mutable T m_hitTol;

    mutable std::vector<gsVector<T> > m_lagBuffer;
    mutable gsMatrix<T> m_alphaVals;
    mutable Stats m_stats;

}; // class gsMomentRule

} // namespace gismo
