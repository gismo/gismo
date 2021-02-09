/** @file gsPatchRule.h

    @brief Provides patch-wise quadrature rule

    Reference:
    Johannessen, K. A. (2017). Optimal quadrature for univariate and tensor product splines.
    Computer Methods in Applied Mechanics and Engineering, 316, 84–99.
    https://doi.org/10.1016/j.cma.2016.04.030

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
    m_over(false)
    {};

    /**
     * @brief      Initialize a (tensor-product) patch-rule with based on the \a basis
     *
     * @param[in]  basis          The basis
     * @param[in]  degree         The degree of the target space
     * @param[in]  regularity     The regularity of the target space
     * @param[in]  overintegrate  Over-integrate or not?
     */
    gsPatchRule(const gsBasis<T> & basis,
                const index_t degree,
                const index_t regularity,
                const bool overintegrate);

    /// Make function returning a smart pointer

    /**
     * @brief      Construct a smart-pointer to the quadrature rule
     *
     * @param[in]  basis          The basis
     * @param[in]  degree         The degree of the target space
     * @param[in]  regularity     The regularity of the target space
     * @param[in]  overintegrate  Over-integrate or not?
     *
     * @return     QuadRule pointer
     */
    static uPtr make(   const gsBasis<T> & basis,
                        const index_t degree,
                        const index_t regularity,
                        const bool overintegrate)
    { return uPtr( new gsPatchRule(basis,degree,regularity,overintegrate) ); }


    /**
     * @brief      Destructor
     */
    ~gsPatchRule() { };

public:

    /**
     * @brief      Maps the points in the d-dimensional cube with points lower and upper
     *
     * @param[in]  lower    The lower corner
     * @param[in]  upper    The upper corner
     * @param      nodes    Quadrature points
     * @param      weights  Quadrature weights
     */
    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const override;

    index_t dim() const { return m_basis->dim(); }

    gsVector<T> nodes(const index_t d=0) const { return m_nodes.at(d); }
    gsVector<T> weights(const index_t d=0) const { return m_weights.at(d); }

protected:

    /**
     * @brief      Computes the knot vector.
     *
     * @param[in]  Bbasis  a gsBSplineBasis
     *
     * @return     knot vector
     */
    gsKnotVector<T> _init(const gsBSplineBasis<T> * Bbasis) const;

    /**
     * @brief      Integrates all basis functions from a knot vector (numerically with Gauss)
     *
     * @param      knots  The knots
     *
     * @return     Greville points and Exact integrals of the basis
     */
    std::pair<gsMatrix<T>,gsMatrix<T>> _integrate(const gsKnotVector<T> & knots ) const;

    /**
     * @brief      Computes the points and weights for the patch rule
     *
     * @param      knots      The knots of the basis
     * @param      greville   The greville point
     * @param      integrals  The exact integrals
     * @param[in]  tol        The tolerance
     *
     * @return     knots and weights
     */
    std::pair<gsVector<T>,gsVector<T>> _compute(const gsKnotVector<T> & knots,
                                                const gsMatrix<T> & greville,
                                                const gsVector<T> & integrals,
                                                const T tol = 1e-10) const;

private:
    const gsBasis<T> * m_basis;
    const index_t m_deg,m_reg;
    const bool m_over;

    std::vector<gsVector<T>> m_nodes;
    std::vector<gsVector<T>> m_weights;

    mutable typename gsSparseSolver<T>::QR m_solver;

    mutable size_t m_dim;

    std::vector<std::map<T,T>> m_maps;
}; // class gsPatchRule

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPatchRule.hpp)
#endif
