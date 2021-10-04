/** @file gsQuadRule.hpp

    @brief Provides implementation of base class for a quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{


template<class T> void
gsQuadRule<T>::mapTo( T startVal, T endVal,
                      gsMatrix<T> & nodes, gsVector<T> & weights ) const
{
    GISMO_ASSERT( 1 == m_nodes.rows(), "Inconsistent quadrature mapping");

    const T h = (endVal-startVal) / T(2);

    // Linear map from [-1,1]^d to [startVal,endVal]
    nodes = (h * (m_nodes.array()+1)) + startVal;

    // Adjust the weights (multiply by the Jacobian of the linear map)
    weights.noalias() = (0==h?T(0.5):h) * m_weights;
}

// Note: left here for inlining
template<class T> void
gsQuadRule<T>::mapTo( const gsMatrix<T>& ab,
                      gsMatrix<T> & nodes) const
{
    GISMO_ASSERT( ab.rows() == m_nodes.rows(), "Inconsistent quadrature mapping");
    nodes.resize( m_nodes.rows(), m_nodes.cols() );
    nodes.setZero();
    const gsVector<T> h = (ab.col(1)-ab.col(0)) / T(2) ;
    nodes.noalias() = ( h.asDiagonal() * (m_nodes.array()+1).matrix() ).colwise() + ab.col(0);
}


template<class T> void
gsQuadRule<T>::mapToAll( const std::vector<T> & breaks,
                         gsMatrix<T> & nodes, gsVector<T> & weights ) const
{
    GISMO_ASSERT( 1 == m_nodes.rows(), "Inconsistent quadrature mapping.");
    GISMO_ASSERT( breaks.size()>1, "At least 2 breaks are needed.");

    const size_t nint    = breaks.size() - 1;
    const index_t nnodes = numNodes();

    nodes  .resize(1, nint*nnodes);
    weights.resize( nint*nnodes );

    for ( size_t i = 0; i!=nint; ++i)
    {
        const T startVal = breaks[i ];
        const T endVal   = breaks[i+1];
        const T h = (endVal-startVal) / T(2);

        // Linear map from [-1,1]^d to [startVal,endVal]
        nodes.middleCols(i*nnodes,nnodes) = (h * (m_nodes.array()+1)) + startVal;

        // Adjust the weights (multiply by the Jacobian of the linear map)
        weights.segment(i*nnodes,nnodes)  = (0==h?T(0.5):h) * m_weights;
    }
}


template<class T> void
gsQuadRule<T>::computeTensorProductRule(const std::vector<gsVector<T> > & nodes,
                                        const std::vector<gsVector<T> > & weights)
{
    this->computeTensorProductRule_into(nodes,weights,m_nodes,m_weights);
}

template<class T> void
gsQuadRule<T>::computeTensorProductRule_into(const std::vector<gsVector<T> > & nodes,
                                             const std::vector<gsVector<T> > & weights,
                                             gsMatrix<T> & targetNodes,
                                             gsVector<T> & targetWeights) const
{
    const short_t d  = static_cast<short_t>(nodes.size());
    GISMO_ASSERT( static_cast<size_t>(d) == weights.size(),
                  "Nodes and weights do not agree." );

    // compute the tensor quadrature rule
    gsPointGrid(nodes, targetNodes);

    gsVector<index_t> numNodes(d);
    for( short_t i=0; i<d; ++i )
        numNodes[i] = weights[i].rows();

    GISMO_ASSERT( targetNodes.cols() == numNodes.prod(),
                  "Inconsistent sizes in nodes and weights.");

    // Compute weight products
    targetWeights.resize( targetNodes.cols() );
    size_t r = 0;
    gsVector<index_t> curr(d);
    curr.setZero();
    do {
        targetWeights[r] = weights[0][curr[0]];
        for (short_t i=1; i<d; ++i)
            targetWeights[r] *= weights[i][curr[i]];
        ++r;
    } while (nextLexicographic(curr, numNodes));
}


} // namespace gismo
