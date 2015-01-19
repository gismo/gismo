/** @file gsQuadRule.hpp

    @brief Provides implementation of base class for a quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/
 
#pragma once


namespace gismo
{


template<class T> void
gsQuadRule<T>::mapTo( T startVal, T endVal,
                      gsMatrix<T> & nodes, gsVector<T> & weights ) const
{
    GISMO_ASSERT( 1 == m_nodes.rows(), "Inconsistent quadrature mapping");
    
    // the factor 0.5 is due to the fact that the one-dimensional reference interval is [-1,1].
    const T h = ( startVal != endVal ? 0.5 * (endVal-startVal) : T(0.5) );
  
    // Linear map from [-1,1]^d to [lower,upper]
    nodes  = (h * m_nodes).array() + 0.5*(startVal+endVal);
    
    // Adjust the weights (multiply by the Jacobian of the linear map)
    weights.noalias() = h * m_weights;
}


template<class T> void
gsQuadRule<T>::computeTensorProductRule(const std::vector<gsVector<T> > & nodes, 
                                        const std::vector<gsVector<T> > & weights)
{
    const int d  = nodes.size();
    GISMO_ASSERT( static_cast<std::size_t>(d) == weights.size(), "Nodes and weights do not agree." );

    // compute the tensor quadrature rule
    gsPointGrid(nodes, m_nodes);
    
    gsVector<index_t> numNodes(d);
    for( int i=0; i<d; ++i )
        numNodes[i] = weights[i].rows();

    GISMO_ASSERT( m_nodes.cols() == numNodes.prod(), "Inconsistent sizes in nodes and weights.");
    
    // Compute weight products
    m_weights.resize( m_nodes.cols() );
    unsigned r = 0;
    gsVector<index_t> curr(d);
    curr.setZero();
    do {
        m_weights[r] = weights[0][curr[0]];
        for (int i=1; i<d; ++i)
            m_weights[r] *= weights[i][curr[i]];
        ++r;
    } while (nextLexicographic(curr, numNodes));
}


} // namespace gismo
