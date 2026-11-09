/** @file gsQuasiInterpolate.h

    @brief Different Quasi-Interpolation Schemes based on the article
    "Spline methods (Lyche Morken)"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner, A. Mantzaflaris, H. Verhelst
**/


#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsHSplines/gsRationalTHBSplineBasis.h>

namespace gismo {

/** \brief Quasi-interpolation operators

    The struct gsQuasiInterpolate has three public member functions to
    use. These functions are implementations of different quasi
    interpolation methods, described in "Spline methods (Lyche
    Morken)" \cite splinemethods2008.  They take a function and
    approximate it via a B-Spline function, whose basis you have to
    provide. More details can be found in the description of the
    respective implementations and in \cite splinemethods2008.

    \tparam T coefficient type
 */

template <class T>
struct gsQuasiInterpolate
{

    /*
    enum method
    {
        Taylor     = 1, ///< Taylor
        Schoenberg = 2, ///< Schoenberg
        EvalBased  = 3
    };

    switch(type)
    {
    case(1): gsQuasiInterpolate<T>::Schoenberg(*basis, fun, coefs); expConvRate = 2.0; break;
    case(2): gsQuasiInterpolate<T>::Taylor(*basis, fun, deg, coefs); expConvRate = (deg+1); break;
    case(3): gsQuasiInterpolate<T>::EvalBased(*basis, fun, true, coefs); expConvRate = (deg+1); break;
    case(0): gsQuasiInterpolate<T>::EvalBased(*basis, fun, false, coefs); expConvRate = (deg+1); break;
    default: GISMO_ERROR("invalid option");
    }
    */

    static gsMatrix<T> localIntpl(const gsBasis<T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i,
                                  const gsMatrix<T> &ab);

    /// Per-index local interpolation. Handles tensor, THB (exact reproduction)
    /// and rational THB bases. A non-truncated hierarchical (HB) basis has no
    /// per-index quasi-interpolant and throws: use the bulk overload.
    static gsMatrix<T> localIntpl(const gsBasis<T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i);

    /// Per-index local interpolation on a hierarchical basis. Throws for a
    /// non-truncated hierarchical (HB) basis: use the bulk overload.
    template<short_t d>
    static gsMatrix<T> localIntpl(const gsHTensorBasis<d,T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i);

    /**
     * @brief Local-interpolation quasi-interpolant of \a fun on the whole basis \a b.
     *
     * Tensor, NURBS, THB and rational THB bases use the per-index local
     * interpolation (for THB: Speleers & Manni, Numer. Math. 132 (2016) 155-184),
     * evaluated in parallel over the indices. A non-truncated hierarchical (HB)
     * basis uses the level-by-level residual scheme of \ref localIntplHB, which
     * is not expressible per index. A rational basis over an HB basis is not
     * supported and throws.
     *
     * @param b      the basis
     * @param fun    the function to approximate
     * @param result the coefficients, size b.size() x fun.targetDim()
     */
    static void localIntpl(const gsBasis<T> &b,
                           const gsFunction<T> &fun,
                           gsMatrix<T> &result);
    
    /// \brief Dimension-independent per-coefficient Taylor QI on a
    /// tensor B-spline basis. Computes the \a j-th coefficient as the
    /// tensor product of the univariate Taylor quasi-interpolants (see
    /// \ref Taylor). \tparam d is the parameter-domain dimension.
    ///
    /// The supported bases are the tensor B-spline bases
    /// gsTensorBSplineBasis<d,T>, d = 1..4 (including gsBSplineBasis). The
    /// per-coefficient Taylor functional of Lyche-Morken Thm 8.5 is defined for a
    /// single tensor knot vector; applied to the level-l tensor basis of a
    /// hierarchical function it is not a dual functional of the hierarchical
    /// basis and gives a wrong result. Every other basis therefore throws, in
    /// the per-index overloads and in the bulk overload below (which throws
    /// before its parallel loop); the gsHTensorBasis overload always throws. Use
    /// localIntpl for hierarchical and rational bases.
    ///
    /// Trap: needs derivatives of \a fun of total order
    /// \f$\sum_{\mathrm{dir}}\min(r,p_{\mathrm{dir}})\f$ at the anchor;
    /// gsFunctionExpr provides only up to order 2.
    ///
    /// Cost per coefficient: one evalAllDers_into at one point plus
    /// \f$\prod_{\mathrm{dir}}(\min(r,p_{\mathrm{dir}})+1)\f$ terms.
    template<short_t d>
    static gsMatrix<T> localTaylor(const gsTensorBSplineBasis<d,T> &b,
                                  const gsFunction<T> &fun,
                                  const index_t &r,
                                  index_t j);

    static gsMatrix<T> localTaylor(const gsBasis<T> &b,
                                const gsFunction<T>  &fun,
                                const index_t &r,
                                index_t i);

    template<short_t d>
    static gsMatrix<T> localTaylor(const gsHTensorBasis<d,T> &b,
                                const gsFunction<T>  &fun,
                                const index_t &r,
                                index_t i);

    /// \brief Bulk Taylor QI on a tensor B-spline basis, d = 1..4; throws on
    /// any other basis before the parallel loop (see the per-coefficient
    /// overload above for the supported bases, the reason and the cost).
    static void localTaylor(const gsBasis<T> &b,
                        const gsFunction<T>  &fun,
                        const index_t &r,
                        gsMatrix<T> & result);    
    
    
    /// \brief Local L2 projection coefficient of function \a i on the element
    /// \a ab, with the active functions of \a b.
    static gsMatrix<T> localL2(const gsBasis<T> &b,   
                                const gsFunction<T>  &source,                                              
                                index_t i,
                                const gsMatrix<T> &ab);

    
    /// \brief Local L2 projection coefficient of function \a i, dispatching on
    /// the basis type (tensor, NURBS, THB, rational THB). Throws on a
    /// non-truncated hierarchical (HB) basis; use the bulk overload.
    static gsMatrix<T> localL2(const gsBasis<T> &b,
                                const gsFunction<T>  &source,
                                index_t i);

    /// \brief Local L2 projection coefficient of function \a i of a
    /// hierarchical basis, with the level-l tensor basis on
    /// elementInSupportOf(i). Throws on an HB basis.
    template<short_t d>
    static gsMatrix<T> localL2(const gsHTensorBasis<d,T> &b,
                                const gsFunction<T>  &source,
                                index_t i);

    /**
     * @brief Quasi-interpolation by local L2 projection.
     *
     * Supported bases and their paths:
     * - tensor and NURBS: local projection with the active functions on
     *   elementInSupportOf(i);
     * - THB: level-l tensor basis on elementInSupportOf(i), as in
     *   Speleers & Manni, Numer. Math. 132 (2016) 155-184 for the
     *   interpolation variant;
     * - rational THB: localL2Rational;
     * - HB (non-truncated): localL2HB, a level-by-level residual loop;
     * - rational over HB: throws.
     *
     * The local functional is a discrete L2 projection with a Gauss-Lobatto
     * rule of p+1 nodes per direction. That rule is exact only to degree 2p-1,
     * so the local mass matrix is not exactly integrated. Reproduction of the
     * space holds because the \f$(p+1)^d\f$ tensor Lobatto nodes are unisolvent
     * for the \f$(p+1)^d\f$ local tensor functions: the collocation matrix B is
     * square and invertible, so \f$M^{-1}BW f^T\f$ equals interpolation at the
     * Lobatto nodes.
     *
     * Cost: per function one dense LU of a \f$(p+1)^d\times(p+1)^d\f$ matrix,
     * \f$O((p+1)^{3d})\f$.
     *
     * @param b      the basis
     * @param source the function to approximate
     * @param[out] result the coefficients, size b.size() x source.targetDim()
     */
    static void localL2(const gsBasis<T> &b,
                        const gsFunction<T>  &source,
                        gsMatrix<T> & result);


    /** \brief A quasi-interpolation scheme based on the tayor expansion of the function to approximate.
     *  See Theorem 8.5 of "Spline methods (Lyche Morken)"
     *  Theorem: (Lyche, Morken: Thm 8.5, page 178)
     *  Let \f$p\f$ and \f$\boldsymbol{\tau}\f$ be the degree and knotvector of the quasi-interpolant, respectively.
     *   Futhermore let \f$r\f$ be an integer with \f$ 0 \le r \le p \f$ and let \f$x_j\f$ be a number in \f$[\tau_j,
     *   \tau_{j+p+1}]\f$ for \f$j=1,\dots,n\f$. Consider the quasi-interpolant
     *   \f[
     *   Q_{p,r}~f=\sum\limits_{j=1}^n{\lambda_j(f)B_{j,p}}, \quad \text{where} \quad
     *   \lambda_j(f) = \frac{1}{p!}\sum\limits_{k=0}^r{(-1)^kD^{p-k}\rho_{j,p}(x_j)D^kf(x_j)}
     *   \f]
     *   and \f$\rho_{j,p}(y) = (y-\tau_{j+1}) \cdots (y - \tau_{j+p})\f$.
     *   Then \f$Q_{p,r}\f$ reproduces all polynomials of degree \f$r\f$ and \f$Q_{p,p}\f$ reproduces all splines
     *   in \f$\mathbb{S}_{p,\tau}\f$.
     *
     * \param b     the B-spline basis of the interpolant (knots and degree)
     * \param fun   a function to approximate
     * \param r     an integer in [0,deg] (order of maximal derivatives of the function)
     *
     * 1-D gsBSplineBasis only; throws otherwise.
     *
     * \param[out] result   a B-spline function, that approximates the given function
     */
    static void Taylor(const gsBasis<T> &bb, const gsFunction<T> &fun, const index_t &r, gsMatrix<T> &result);

     /** \brief A quasi-interpolation scheme based on the tayor expansion of the function to approximate.
     *  See Theorem 8.5 of "Spline methods (Lyche Morken)"
     *  Theorem: (Lyche, Morken: Thm 8.5, page 178)
     *  Let \f$p\f$ and \f$\boldsymbol{\tau}\f$ be the degree and knotvector of the quasi-interpolant, respectively.
     *   Futhermore let \f$r\f$ be an integer with \f$ 0 \le r \le p \f$ and let \f$x_j\f$ be a number in \f$[\tau_j,
     *   \tau_{j+p+1}]\f$ for \f$j=1,\dots,n\f$. Consider the quasi-interpolant
     *   \f[
     *   Q_{p,r}~f=\sum\limits_{j=1}^n{\lambda_j(f)B_{j,p}}, \quad \text{where} \quad
     *   \lambda_j(f) = \frac{1}{p!}\sum\limits_{k=0}^r{(-1)^kD^{p-k}\rho_{j,p}(x_j)D^kf(x_j)}
     *   \f]
     *   and \f$\rho_{j,p}(y) = (y-\tau_{j+1}) \cdots (y - \tau_{j+p})\f$.
     *   Then \f$Q_{p,r}\f$ reproduces all polynomials of degree \f$r\f$ and \f$Q_{p,p}\f$ reproduces all splines
     *   in \f$\mathbb{S}_{p,\tau}\f$.
     *
     * \param b     the B-spline basis of the interpolant (knots and degree)
     * \param fun   a function to approximate
     * \param r     an integer in [0,deg] (order of maximal derivatives of the function)
     *
     * Forwards to localTaylor(bb, fun, r, result): tensor B-spline bases only.
     *
     * \param[out] result   a B-spline function, that approximates the given function
     */
    static void Taylor2D(const gsBasis<T> &bb, const gsFunction<T> &fun, const index_t &r, gsMatrix<T> &result);

    /**
     * @brief A quasi-interpolation scheme based on Schoenberg Variation Diminishing Spline Approximation.
     * See Exercise 9.1 of "Spline methods (Lyche Morken)"
     * @param b     the B-spline basis of the interpolant (knots and degree)
     * @param fun   a function to approximate
     * @param[out] result   a B-spline function, that approximates the given function
     */
    static void Schoenberg(const gsBasis<T> &b, const gsFunction<T> &fun,
                           gsMatrix<T> &result);

    static gsMatrix<T> Schoenberg(const gsBasis<T> &b, const gsFunction<T> &fun, index_t i);


    /**
     * @brief A quasi-interpolation scheme based on the evaluation of the function at certain points.
     * See sections 8.2.1, 8.2.2, 8.2.3 and Theorem 8.7 or Lemma 9.7 of "Spline methods (Lyche Morken)". The formulas for the special cases (degrees 1, 2 and 3) look like this:
     * \f[
     * P_{deg}~f(x) = \sum_{j=1}^n {\lambda_j(f) B_j(x)}
     * \f]
     * where the coefficients \f$ \lambda_j\f$ for \f$ j=1,\dots,n\f$ are given as:
     * \f[
     * \lambda_j(f) = f(\tau_{j+1})
     * \f]
     * for degree 1,
     * \f[
     * \lambda_j(f) = \begin{cases}
     * f(\tau_1) &\mbox{if } j=1; \\
     * \frac{1}{2} (-f(x_{j,0}) + 4f(x_{j,1}) - f(x_{j,2}) ), &\mbox{if } 1<j<n; \\
     * f(\tau_{n+1}) &\mbox{if } j=n; \end{cases}
     * \f]
     * where \f$ x_{j,0} = \tau_{j+1}, \quad x_{j,1} = \frac{\tau_{j+1}+\tau_{j+2}}{2}, \quad x_{j,2} = \tau_{j+2} \f$
     *
     * for degree 2 and
     * \f[
     * \lambda_j(f) = \begin{cases}
     * f(\tau_4) &\mbox{if} j=1; \\
     * \frac{1}{18}(-5f(\tau_4)+40f(\tau_{9/2})-24f(\tau_5)+8f(\tau_{11/2})-f(\tau_6)) &\mbox{if } j=2; \\
     * \frac{1}{6} (f(\tau_{j+1}) -8f(\tau_{j+3/2}) +20 f(\tau_{j+2}) -8f(\tau_{j+5/2})+f(\tau_{j+3})), &\mbox{if } 2<j<n-1; \\
     * \frac{1}{18}(-f(\tau_{n-1})+8f(\tau_{n-1/2})-24f(\tau_n)+40f(\tau_{n+1/2})-5f(\tau_{n+1})) &\mbox{if } j=n-1; \\
     * f(\tau_{n+1}) &\mbox{if } j=n; \end{cases}
     * \f]
     * where \f$ \tau_{j+k/2} = \frac{\tau_{j+(k-1)/2}+\tau_{j+(k+1)/2}}{2} \f$,
     *
     * for degree 3.
     *
     * Theorem 8.7:
     *
     * Let \f$ \mathbb{S}_{p,\mathbf{\tau}} \f$ be a spline space with a \f$p+1\f$-regular knot vector \f$ \tau = (\mathbf{\tau}_i)_{i=1}^{n+p+1} \f$.
     * Let \f$ (x_{j,k})_{k=0}^r \f$ be \f$r+1\f$ distinct points in \f$ [\tau_j,\tau_{j+p+1}] \f$ for \f$ j=1, \dots, n \f$ and let \f$\omega_{j,k} \f$
     * be the j-th B-spline coefficient of the polynomial
     * \f[
     * p_{j,k}(x) = \prod_{s=0, s\ne k}^r {\frac{x-x_{j,s}}{x_{j,k}-x_{j,s}}}.
     * \f]
     * Then \f$P_{p,p}~f = f \f$ for all \f$ f \in \tau_r \f$ and if \f$r=p\f$ and all the numbers \f$(x_{j,k})_{k=0}^r \f$ lie in one subinterval
     * \f[
     * \tau_j \le \tau_{\ell_j} \le x_{j,0} < x_{j,1} < \cdots < x_{j,r} \le \tau_{\ell_j +1} \le \tau_{j+p+1}
     * \f]
     * then \f$P_{p,p}~f = f\f$ for all \f$ f \in \mathbb{S}_{p,\mathbf{\tau}} \f$.
     * @param b     the B-spline basis of the interpolant (knots and degree)
     * @param fun   a function to approximate,
     * @param specialCase   if set to true, use the special implementations for degrees 1, 2 and 3;
     * if set to false, use the general implementation
     * @param result    a B-spline function, that approximates the given function
     */
    static void EvalBased(const gsBasis<T> &bb, const gsFunction<T> &fun, const bool specialCase, gsMatrix<T> &result);


    /*
    static void qiCwiseData(const gsTensorBSplineBasis<T,2> & tbsp,
                            const gsVector<index_t> & ind,
                            std::vector<gsMatrix<T> > & qiNodes,
                            std::vector<gsMatrix<T> > & qiWeights);

    static void compute(const gsTensorBSplineBasis<T,2> & tbsp,
                        const gsFunction<T> & fun,
                        gsTensorBSpline<T> & res);
    */

protected:

    /**
     * @brief Level-by-level residual quasi-interpolation on a non-truncated
     * hierarchical (HB) basis.
     *
     * Levels l = 0, 1, ... are processed in ascending order. With \f$s_{<l}\f$
     * the spline built from the coefficients of levels below l, the level-l
     * coefficient of function i is obtained by interpolating the residual
     * \f$f - s_{<l}\f$ with the level-l tensor basis at Gauss points on a
     * level-l cell Q of \f$\Omega^l\setminus\Omega^{l+1}\f$ inside supp(i).
     *
     * On Q every active HB function of level > l vanishes, since it is supported
     * in \f$\Omega^{l+1}\f$. Hence for f in the HB space,
     * \f$f - s_{<l}\f$ on Q equals \f$\sum_{\mathrm{level}(j)=l} c_j\beta_j|_Q\f$, a
     * level-l tensor polynomial on Q, and the local interpolation recovers
     * \f$c_i\f$ exactly. The operator is therefore a projector onto the HB space.
     * This scheme is derived here and verified numerically; no literature
     * reference is attached to it.
     *
     * For data outside the space the result coincides with the THB quasi-interpolant
     * on the same hierarchy; this is a measured fact (<= 3.4e-13 on the tested
     * hierarchies), not a derived one: two projectors onto the same space need
     * only agree on that space.
     *
     * Cost: per function one dense LU of a \f$(p+1)^d\times(p+1)^d\f$ matrix,
     * \f$O((p+1)^{3d})\f$, plus evaluations of f, the level-l tensor basis and
     * \f$s_{<l}\f$ at \f$(p+1)^d\f$ points. Per level one makeGeometry
     * (\f$O(n\cdot\mathrm{targetDim})\f$ copy and a basis clone). Levels are
     * sequential, the functions of one level are processed in parallel.
     */
    template<short_t d>
    static void localIntplHB(const gsHTensorBasis<d,T> & b,
                             const gsFunction<T> & fun,
                             gsMatrix<T> & result);

    /**
     * @brief Local interpolation coefficient on a rational THB-spline basis.
     *
     * A rational THB spline is \f$s = \sum_k c_k w_k T_k / W\f$ with
     * \f$W = \sum_k w_k T_k\f$ and T_k the source THB functions. Hence
     * \f$g = sW = \sum_k (c_k w_k) T_k\f$ is a THB spline, and by preservation of
     * coefficients its level-l coefficient on a cell in \f$\Omega^l\setminus\Omega^{l+1}\f$
     * of active function i equals \f$c_i w_i\f$. The function \a fun times \f$W\f$ is
     * interpolated on that cell with the level-l tensor basis (as in the hierarchical
     * quasi-interpolant of Speleers & Manni, Numer. Math. 132 (2016) 155-184; see also
     * Giannelli, Juettler, Speleers 2014) and the result is divided by \f$w_i\f$.
     * Weights are assumed nonzero.
     * Accuracy degrades with the spread of the weights, since the coefficient is
     * obtained by dividing by \f$w_i\f$ (on a p=3 test mesh: max/min weight ratio 10
     * gives ~1e-12, 100 gives ~6e-11, 1000 gives ~2e-9).
     */
    template<short_t d>
    static gsMatrix<T> localIntplRational(const gsRationalBasis<gsTHBSplineBasis<d,T> > & b,
                                          const gsFunction<T> & fun,
                                          index_t i);

    /**
     * @brief Level-by-level residual loop shared by the HB quasi-interpolants.
     *
     * \a L2 = false gives local interpolation (Gauss-Legendre points), \a L2 =
     * true the local L2 functional (Gauss-Lobatto points). The functions of one
     * level are processed in parallel, the levels sequentially. See
     * localIntplHB for the theory.
     */
    template<short_t d, bool L2>
    static void localHBLevelLoop(const gsHTensorBasis<d,T> & b,
                                 const gsFunction<T> & fun,
                                 gsMatrix<T> & result);

    /**
     * @brief Local L2 projection on a non-truncated hierarchical (HB) basis.
     *
     * Same scheme as localIntplHB with the local L2 functional of localL2 in
     * place of local interpolation. On a level-l cell
     * \f$Q\subset\Omega^l\setminus\Omega^{l+1}\f$ the HB functions of level
     * > l vanish, so \f$f - s_{<l}\f$ on Q is a level-l tensor polynomial and
     * the projection recovers \f$c_i\f$. The Lobatto nodes lie on
     * \f$\partial Q\f$, so the vanishing of the level > l functions there needs
     * continuity across \f$\partial Q\f$, i.e. p >= 1 in every direction.
     *
     * Cost: as localIntplHB.
     */
    template<short_t d>
    static void localL2HB(const gsHTensorBasis<d,T> & b,
                          const gsFunction<T> & fun,
                          gsMatrix<T> & result);

    /**
     * @brief Local L2 projection coefficient on a rational THB-spline basis.
     *
     * With \f$W = \sum_k w_k T_k\f$, the function \f$g = sW\f$ is a THB spline
     * with coefficients \f$c_k w_k\f$, so its level-l coefficient on a cell in
     * \f$\Omega^l\setminus\Omega^{l+1}\f$ of active function i equals
     * \f$c_i w_i\f$. \a fun times \f$W\f$ is projected on that cell with the
     * level-l tensor basis (local L2 functional of localL2) and the result is
     * divided by \f$w_i\f$. Weights must be nonzero; accuracy degrades with the
     * weight spread, since the coefficient is obtained by dividing by \f$w_i\f$.
     */
    template<short_t d>
    static gsMatrix<T> localL2Rational(const gsRationalBasis<gsTHBSplineBasis<d,T> > & b,
                                       const gsFunction<T> & fun,
                                       index_t i);

    /// True for a rational basis whose source is a non-truncated hierarchical (HB) basis.
    static bool isRationalOverHB(const gsBasis<T> & b);

    /**
     * @brief Compute the derivative of a certain order of a normalized polynomial (leading coefficient is 1) defined by its roots at a given point.
     *  \f$g(y) = (y-y_1) \cdots (y-y_n)\f$, where \f$y_1,\dots,y_n\f$ are the roots of the polynomial.
     * @param zeros roots of the polynomial
     * @param order the order of the derivative to compute
     * @param x     evaluation point
     * @return      the value of the derivative, at the given point, \f$D^\alpha g(x)\f$, where \f$\alpha\f$ is the given order.
     */
    static T derivProd(const std::vector<T> &zeros, const index_t &order, const T &x);


    /**
     * @brief Row index, within \c derivs[|alpha|] of \ref
     * gsFunctionSet::evalAllDers_into, of the mixed partial derivative
     * \f$ \partial^\alpha f^{(comp)} \f$ for a function of domain
     * dimension \a d.
     *
     * Encodes the packing convention of \c evalAllDers_into: per target
     * component the block holds, for order \f$m=|\alpha|\f$: the value
     * (m=0); the first derivatives \f$\partial_0,\dots,\partial_{d-1}\f$
     * (m=1); for m=2 the pure second derivatives
     * \f$\partial_{00},\dots,\partial_{d-1,d-1}\f$ first, then the mixed
     * ones \f$\partial_{ab}\f$ (a<b) in lexicographic order; and for
     * \f$m\ge 3\f$ the derivatives in composition (lexicographic) order,
     * see \ref nextComposition.
     *
     * @param alpha per-direction derivative orders (size \a d)
     * @param d     domain dimension
     * @param comp  target component index
     * @return      the row of \f$\partial^\alpha f^{(comp)}\f$
     */
    static index_t derivRow(const gsVector<index_t> &alpha, short_t d, index_t comp);


    /**
     * @brief Compute a number of equally distributed points in a given interval \f$[a,b]\f$.
     * You get a list of points \f$\{a, a+(b-a)\frac{1}{n-1}, \dots, a+(b-a)\frac{n-2}{n-1}, b\}\f$.
     * @param a start value of the interval
     * @param b end value of the interval
     * @param n number of points
     * @param[out] computed points
     */
    static void distributePoints(T a, T b, int n, gsMatrix<T> &points);


    /**
     * @brief To compute the control points \f$ \lambda_i(f) = \sum\limits_{k=0}^p{\omega_{i,k}f(x_{i,k})} \f$
     * of the quasi-interpolant one uses the function computeControlPoints. The weights \f$ \omega_{i,k} \f$ can be
     * computed as \f$\omega_{i,k} = \gamma_i(p_{i,k})\f$, for \f$k=0,1,\dots,p\f$, where
     * \f[ \gamma_i(g) = \frac{1}{p!}\sum\limits_{(j_1,\dots,j_p)\in \mathcal{P}_p}{(\tau_{i+j_1}-v_1)\cdots(\tau_{i+j_p}-vp)},\f]
     * for a polynomial \f$g(x) = (x-v_1) \cdots (x-v_p)\f$, where \f$\mathcal{P}_p\f$ is the set of all permutations of the intergers \f$\{1,2,\dots,p\}\f$.
     * @param points    the points \f$x_{i,k}\f$  of the above formula
     * @param knots     the knotvector of the quasi-interpolant
     * @param pos       the index i of the above formula
     * @param[out] weights   the computed weights \f$\omega_{i,k}\f$ of the above formula
     */
    static void computeWeights(const gsMatrix<T> &points, const gsKnotVector<T> &knots, const index_t &pos, gsMatrix<T> &weights);


    /**
     * @brief The quasi-interpolant is a spline function, in particular a linear combination of some controlpoints and the B-spline basis functions.
     *  \f$Q_p~f(x) = \sum\limits_{i=1}^n{\lambda_i(f)B_{i,p(x)}}\f$ where the controlpoints can be computed as \f$\lambda_i(f) = \sum\limits_{k=0}^p{\omega_{i,k}f(x_i,k)}\f$.
     *  The points \f$x_{i,k}\f$ are equally distributed points in the largest subinterval of \f$[\tau_{i+1}, \tau_{i+p}]\f$.
     * @param weights   the weights \f$\omega_i,k\f$ of the above formula
     * @param fun       the function to approximate, \f$f\f$ of the above formula
     * @param xik       the points \f$x_{i,k}\f$ of the above formula
     * @return      the computed control point  \f$\lambda_i(f)\f$
     */
    static gsMatrix<T> computeControlPoints(const gsMatrix<T> &weights, const gsFunction<T> &fun, const gsMatrix<T> &xik);

    /**
     * @brief This function finds the greatest knot interval in a given range in a knot vector.
     * @param knots     the knot vector
     * @param posStart  the index of the left knot of the first interval to be considered
     * @param posEnd    the index of the right knot of the last interval to be considers
     * @return          the index of the left knot of the largest knot interval
     */
    static int greatestSubInterval(const gsKnotVector<T> &knots, const index_t &posStart, const index_t &posEnd);


}; //struct

#ifdef GISMO_WITH_PYBIND11

    /**
     * @brief Initializes the Python wrapper for the class: gsQuasiInterpolate
     */
    void pybind11_init_gsQuasiInterpolate(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // gismo
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsQuasiInterpolate.hpp)
#endif
