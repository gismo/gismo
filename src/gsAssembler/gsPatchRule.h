/** @file gsPatchRule.h

    @brief Provides patch-wise quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsAssembler/gsQuadRule.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>

namespace gismo
{

/**
    \brief Class that represents the (tensor) patch quadrature rule

    The class uses the \a basis to construct the quadrature. The rule uses (and distributes) quadrature points according to a Newton scheme and the number of points and weights is based on the number of points and weights in a target space with a specific \a order and \a regularity.

    Optionally, one can \a over-integrate the quadrature rule, meaning that extra quadrature points are added on the boundary elements. The number of points for over integration is equal to the degree of the knot vector of the basis.

    The method works as follows(see Johannessen 2017 for more information):
    Given a vector of quadrature points \f$\mathbf{\xi}=[ \xi_1, \xi_2,...\xi_n ]\f$ (where \f$n\f$ depends on the order and regularity of the target space, the fact that this space provides an even/odd number of points and on over-integration or not) with corresponding weights \f$\mathbf{w}=[ w_1, w_2,...w_n ]\f$. Initially, we know the integrals of our basis functions \f$ \int_\mathbb{R} N_{j,p}(\xi)\:\text{d}\xi\f$, \f$j=\{1,2,...,m\}\f$, where \f$m\f$ is the number of basis functions in our space. Then, we want to minimize the residual
    \f{eqnarray*}{
    \mathbf{F}
    \left(
        \begin{bmatrix} \mathbf{w} \\ \mathbf{\xi} \end{bmatrix}
    \right)
     =
     \begin{bmatrix}
        N_{1,p}(\xi_1) & N_{1,p}(\xi_2) & \dots & N_{1,p}(\xi_n)\\
        N_{2,p}(\xi_1) & N_{2,p}(\xi_2) & \dots & N_{2,p}(\xi_n)\\
        \vdots &  & \ddots & \vdots \\
        N_{m,p}(\xi_1) & N_{m,p}(\xi_2) & \dots & N_{m,p}(\xi_n)
     \end{bmatrix}
     \begin{bmatrix}
        w_1 \\ w_2 \\ \vdots \\ w_n
     \end{bmatrix}
     -
     \begin{bmatrix}
        \int N_{1,p}(\xi)\:\text{d}\xi\\
        \int N_{2,p}(\xi)\:\text{d}\xi\\
        \vdots \\
        \int N_{m,p}(\xi)\:\text{d}\xi
     \end{bmatrix}
    \f}
    For the Newton iterations, we compute the Jacobian according to \f$\mathbf{F}\f$:
    \f{eqnarray*}{
        \frac{\partial\mathbf{F}}{\partial \mathbf{w}} =
     \begin{bmatrix}
        N_{1,p}(\xi_1) & N_{1,p}(\xi_2) & \dots & N_{1,p}(\xi_n)\\
        N_{2,p}(\xi_1) & N_{2,p}(\xi_2) & \dots & N_{2,p}(\xi_n)\\
        \vdots &  & \ddots & \vdots \\
        N_{m,p}(\xi_1) & N_{m,p}(\xi_2) & \dots & N_{m,p}(\xi_n)
     \end{bmatrix}
    \f}
    \f{eqnarray*}{
        \frac{\partial\mathbf{F}}{\partial \mathbf{\xi}} =
     \begin{bmatrix}
        w_1 N'_{1,p}(\xi_1) & w_1 N'_{1,p}(\xi_2) & \dots & w_1 N'_{1,p}(\xi_n)\\
        w_2 N'_{2,p}(\xi_1) & w_2 N'_{2,p}(\xi_2) & \dots & w_2 N'_{2,p}(\xi_n)\\
        \vdots &  & \ddots & \vdots \\
        w_m N'_{m,p}(\xi_1) & w_m N'_{m,p}(\xi_2) & \dots & w_m N'_{m,p}(\xi_n)
     \end{bmatrix}
    \f}
    Such that \f$ \partial \mathbf{F} = \left[ \frac{\partial\mathbf{F}}{\partial \mathbf{w}} , \frac{\partial\mathbf{F}}{\partial \mathbf{\xi}} \right] \in \mathbb{R}^{m\times m}\f$.

    The weights and quadrature points are updated using the Newton update. As an initial guess, we use \f$ w_i = \int N_{2i,p}(\xi) + N_{2i+1,p}(\xi) \: \text{d}\xi\f$ and \f$ \xi_i=\frac{\tau_{2i}+\tau_{2i+1}}{2}\f$ with \f$\tau_i\f$ the Greville abcissa of basis function \f$i\f$

    Reference:
    Johannessen, K. A. (2017). Optimal quadrature for univariate and tensor product splines.
    Computer Methods in Applied Mechanics and Engineering, 316, 84–99.
    https://doi.org/10.1016/j.cma.2016.04.030

    \ingroup Assembler
*/
template<class T>
class gsPatchRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsPatchRule> uPtr;

    /// Default empty constructor
    gsPatchRule()
    :
    m_basis(nullptr),
    m_deg(0),
    m_reg(0),
    m_over(false),
    m_fixDir(-1)
    {};

    /**
     * @brief      Initialize a (tensor-product) patch-rule with based on the \a basis
     *
     * @param[in]  basis          The basis
     * @param[in]  degree         The degree of the target space
     * @param[in]  regularity     The regularity of the target space
     * @param[in]  overintegrate  Over-integrate (i.e. add p points in boundary elements) when true
     */
    gsPatchRule(const gsBasis<T> & basis,
                const index_t degree,
                const index_t regularity,
                const bool overintegrate,
                const short_t fixDir = -1);

    /// Make function returning a smart pointer

    /**
     * @brief      Construct a smart-pointer to the quadrature rule
     *
     * @param[in]  basis          The basis
     * @param[in]  degree         The degree of the target space
     * @param[in]  regularity     The regularity of the target space
     * @param[in]  overintegrate  Over-integrate (i.e. add p points in boundary elements) when true
     *
     * @return     QuadRule pointer
     */
    static uPtr make(   const gsBasis<T> & basis,
                        const index_t degree,
                        const index_t regularity,
                        const bool overintegrate,
                        const short_t fixDir = -1)
    { return uPtr( new gsPatchRule(basis,degree,regularity,overintegrate,fixDir) ); }


    /**
     * @brief      Destructor
     */
    ~gsPatchRule() { };

public:

    /**
     * @brief      Maps the points in the d-dimensional cube with points \a lower and \a upper
     *
     * See \ref gsQuadRule for documentation
     *
     * @param[in]  lower    The lower corner
     * @param[in]  upper    The upper corner
     * @param      nodes    Quadrature points
     * @param      weights  Quadrature weights
     */
    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const;

    /// See \ref gsQuadRule for documentation
    void mapTo(  T startVal, T endVal,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const;

    /// Not implemented! See \ref gsQuadRule for documentation
    void mapToAll( const std::vector<T> & breaks,
                   gsMatrix<T> & nodes, gsVector<T> & weights ) const
    { GISMO_NO_IMPLEMENTATION }

    index_t dim() const { return m_basis->dim(); }

    gsVector<T> nodes(const index_t d=0) const { return m_nodes.at(d); }
    gsVector<T> weights(const index_t d=0) const { return m_weights.at(d); }

protected:

    /**
     * @brief      Computes the number of quadrature points to exactly integrate the space
     *
     * @param[in]  knots  knot vector
     *
     * @return     Integer number of quadrature points
     */
    index_t _numQuads(const gsKnotVector<T> knots) const
    {
        index_t ndof = (m_deg + 1)*knots.numElements() - (m_reg+1)*(knots.numElements()-1);
        ndof = static_cast<index_t>(std::ceil(ndof/2.0));
        return ndof;
    }

    /**
     * @brief      Determines whether the specified knot vector is symmetric.
     *
     * @param[in]  knots  knot vector
     *
     * @return     True if the specified knots is symmetric, False otherwise.
     */
    bool _isSymmetric(const gsKnotVector<T> knots) const
    {
        gsKnotVector<T> copy = knots;
        copy.reverse();
        std::vector<T> result;
        knots.symDifference(copy,result);
        return result.size()==0 ? true : false;
    }


    /**
     * @brief      Initializes the knot vector for the patch-rule implementation
     *
     * The knot vector for the basis in a specific direction is modified:
     * 1) the number of knots is made even. If this is not the case, an extra knot is added in the middle
     * 2) over-integration is applied. In this case, p knots are added to the boundary elements
     *
     * @param[in]  Bbasis  a gsBSplineBasis
     *
     * @return     knot vector
     */
    gsKnotVector<T> _init(const gsBSplineBasis<T> * Bbasis) const;

    /**
     * @brief      Integrates all basis functions from a knot vector (numerically with Gauss)
     *
     * This computes the integrals of the basis function for the newton integration.
     *
     * @param      knots  The knots
     *
     * @return     Greville points and Exact integrals of the basis
     */
    std::pair<gsMatrix<T>,gsVector<T> > _integrate(const gsKnotVector<T> & knots ) const;

    /**
     * @brief      Computes the points and weights for the patch rule using Newton iterations
     *
     * @param      knots      The knots of the basis
     * @param      greville   The greville point
     * @param      integrals  The exact integrals
     * @param[in]  tol        The tolerance
     *
     * @return     knots and weights
     */
    std::pair<gsVector<T>,gsVector<T> > _compute(const gsKnotVector<T> & knots,
                                                const gsMatrix<T> & greville,
                                                const gsVector<T> & integrals,
                                                const T tol = 1e-10) const;

private:
    const gsBasis<T> * m_basis;
    const index_t m_deg,m_reg;
    const bool m_over;
    const short_t m_fixDir;

    std::vector<gsVector<T> > m_nodes;
    std::vector<gsVector<T> > m_weights;

    gsVector<T> m_end;

    mutable typename gsSparseSolver<T>::QR m_solver;

    size_t m_dim;

    std::vector<std::map<T,T> > m_maps;
}; // class gsPatchRule

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPatchRule.hpp)
#endif
