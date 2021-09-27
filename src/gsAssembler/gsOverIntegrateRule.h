/** @file gsOverIntegrateRule.h

    @brief Over-integrates a Gauss-Legendre or Gauss-Lobatto rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsAssembler/gsQuadRule.h>
#include <gsCore/gsBasis.h>

namespace gismo
{

/**
    \brief Class that defines a mixed quadrature rule with different rules for the interior and the boundaries

    This class is defined using two other quadrature rules (only works for \ref gsGaussRule or \ref gsLobattoRule) for the interior and for the boundary. Depending on the location of the considered element, it will either use the interior rule or the boundary rule. In this way, one can for example use full (exact) integration on the boundary points and reduced integration in the interior.

    \ingroup Assembler
*/
template<class T>
class gsOverIntegrateRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsOverIntegrateRule> uPtr;

    /// Default empty constructor
    gsOverIntegrateRule()
    :
    m_basis(nullptr),
    m_interior(),
    m_boundary()
    {}

    /**
     * @brief      Constructor
     *
     * @param[in]  basis         The basis
     * @param[in]  quadInterior  vector with the univariate INTERIOR quadrature rules for each dimension of the basis.
     * @param[in]  quadBoundary  vector with the univariate BOUNDARY quadrature rules for each dimension of the basis.
     *
     * @note: Only works for QuadRules that do not re-implement mapTo (slicing occurs)
     */
    gsOverIntegrateRule(const  gsBasis<T> & basis,
                        const  std::vector<gsQuadRule<T> > & quadInterior,
                        const  std::vector<gsQuadRule<T> > & quadBoundary)
    :
    m_basis(&basis),
    m_interior(quadInterior),
    m_boundary(quadBoundary)
    {
        std::vector< gsVector<T> > nodes(m_basis->dim());
        m_start = m_basis->support().col(0);
        m_end = m_basis->support().col(1);
    };

    /**
     * @brief      Construct a smart-pointer to the rule
     *
     * @param[in]  basis         The basis
     * @param[in]  quadInterior  The rule used for the interior
     * @param[in]  quadBoundary  The rule used for the boundary
     *
     * @note: Only works for QuadRules that do not re-implement mapTo (slicing occurs)
     */
    static uPtr make(   const  gsBasis<T> & basis,
                        const  std::vector<gsQuadRule<T> > & quadInterior,
                        const  std::vector<gsQuadRule<T> > & quadBoundary )
    { return uPtr( new gsOverIntegrateRule(basis,quadInterior,quadBoundary) ); }

    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsOverIntegrateRule() { }

public:
    /// \brief Dimension of the rule
    index_t dim() const { return m_basis->dim(); }

    /**
     * @brief      Maps the points in the d-dimensional cube with points lower and upper
     *
     * @param[in]  lower    The lower corner
     * @param[in]  upper    The upper corner
     * @param      nodes    Quadrature points
     * @param      weights  Quadrature weights
     */
    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {

        gsVector<T> bot = lower-m_start;
        gsVector<T> top = upper-m_end;
        std::vector<gsVector<T> > elNodes(m_basis->dim());
        std::vector<gsVector<T> > elWeights(m_basis->dim());
        gsMatrix<T> tmp;
        for (index_t d = 0; d!=dim(); d++)
        {
            if (bot[d]==0.0 || top[d]==0.0) // is the coordinate on a side?
            {
                m_boundary[d].mapTo( lower[d], upper[d], tmp, elWeights[d]);
                GISMO_ASSERT(tmp.rows()==1,"Dimension of the nodes is wrong!");
                elNodes[d] = tmp.transpose();
            }
            else
            {
                m_interior[d].mapTo( lower[d], upper[d], tmp, elWeights[d]);
                GISMO_ASSERT(tmp.rows()==1,"Dimension of the nodes is wrong!");
                elNodes[d] = tmp.transpose();
            }
        }
        this->computeTensorProductRule_into(elNodes,elWeights,nodes,weights);
    }

private:
    const gsBasis<T> * m_basis;
    std::vector<gsQuadRule<T> > m_interior, m_boundary;
    gsVector<T> m_start,m_end;

}; // class gsOverIntegrateRule

} // namespace gismo
