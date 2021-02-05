/** @file gsOverIntegrateRule.h

    @brief Provides the Gauss-Legendre quadrature rule

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
    \brief Class that represents the (tensor) patch quadrature rule

    \ingroup Assembler
*/
template<class T>
class gsOverIntegrateRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsOverIntegrateRule()
    :
    m_basis(nullptr),
    m_interior(),
    m_boundary()
    {}

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsOverIntegrateRule(const  gsBasis<T> & basis,
                        const  gsQuadRule<T> & quadInterior,
                        const  gsQuadRule<T> & quadBoundary)
    :
    m_basis(&basis),
    m_interior(quadInterior),
    m_boundary(quadBoundary)
    {
        m_start = m_basis->support().col(0);
        m_end = m_basis->support().col(1);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsOverIntegrateRule() { }

public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes,
                   unsigned digits = 0 )
    {
        m_interior.setNodes(numNodes,digits);
        m_boundary.setNodes(numNodes,digits);
    };

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return m_basis->dim(); }

    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {
        if ((lower-m_start).prod()==0 || (upper-m_end).prod()==0)
            m_boundary.mapTo( lower, upper, nodes, weights);
        else
            m_interior.mapTo( lower, upper, nodes, weights);
    }

private:
    const gsBasis<T> * m_basis;
    gsQuadRule<T> m_interior, m_boundary;
    mutable gsVector<T> m_start,m_end;

}; // class gsOverIntegrateRule

} // namespace gismo
