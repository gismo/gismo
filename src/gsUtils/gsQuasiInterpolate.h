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

    static gsMatrix<T> localIntpl(const gsBasis<T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i);

    template<short_t d>
    static gsMatrix<T> localIntpl(const gsHTensorBasis<d,T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i);


    static void localIntpl(const gsBasis<T> &b,
                           const gsFunction<T> &fun,
                           gsMatrix<T> &result);


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
     * \param[out] result   a B-spline function, that approximates the given function
     */
    static void Taylor(const gsBasis<T> &bb, const gsFunction<T> &fun, const int &r, gsMatrix<T> &result);


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
    //*/

protected:

    /**
     * @brief Compute the derivative of a certain order of a normalized polynomial (leading coefficient is 1) defined by its roots at a given point.
     *  \f$g(y) = (y-y_1) \cdots (y-y_n)\f$, where \f$y_1,\dots,y_n\f$ are the roots of the polynomial.
     * @param zeros roots of the polynomial
     * @param order the order of the derivative to compute
     * @param x     evaluation point
     * @return      the value of the derivative, at the given point, \f$D^\alpha g(x)\f$, where \f$\alpha\f$ is the given order.
     */
    static T derivProd(const std::vector<T> &zeros, const int &order, const T &x);


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
    static void computeWeights(const gsMatrix<T> &points, const gsKnotVector<T> &knots, const int &pos, gsMatrix<T> &weights);


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
    static int greatestSubInterval(const gsKnotVector<T> &knots, const int &posStart, const int &posEnd);


}; //struct

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsQuasiInterpolate.hpp)
#endif
