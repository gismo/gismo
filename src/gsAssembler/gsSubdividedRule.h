/** @file gsSubdividedQuadRule.h

    @brief Allows to apply another quadrature to subsets of the overall integration domain (-1,1)^d

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsAssembler/gsQuadRule.h>

namespace gismo
{

/** 
    \brief Allows to apply another quadrature to subsets of the overall integration domain (-1,1)^d
    
    \ingroup Assembler
*/  
template<class T, class Rule>
class gsSubdividedRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsSubdividedRule() {}

    /// Initialize
    ///
    /// @param rule     Quadrature rule to be applied
    /// @param subdiv   Number of subdomains per direction
    gsSubdividedRule(const Rule& rule, index_t subdiv)
    { init(rule, subdiv); }

    /// Initialize
    ///
    /// @param rule     Quadrature rule to be applied
    /// @param subdiv   Number of subdomains per direction
    void init(const Rule& rule, index_t subdiv)
    {
        m_weights = rule.referenceWeights();
        m_nodes = rule.referenceNodes();
        const index_t d = m_nodes.rows();
        for (index_t k=0; k<d; ++k)
        {
            const index_t n = m_nodes.cols();
            GISMO_ASSERT ( n==m_weights.rows(), "Number of weigts and number of quadrature nodes do not agree." );

            gsVector<T> weights( n * subdiv );
            gsMatrix<T> nodes( d, n * subdiv );

            for (index_t i=0; i<subdiv; ++i)
            {
                for (index_t j=0; j<n; ++j)
                {
                    weights[i*n+j] = m_weights[j] / subdiv;
                    for (index_t l=0; l<d; ++l)
                        if (l==k)
                            nodes(l,i*n+j) = (((T)1+m_nodes(l,j)) / (2*subdiv) + (T)i / subdiv)*2-1;
                        else
                            nodes(l,i*n+j) = m_nodes(l,j);
                }
            }
            m_weights.swap(weights);
            m_nodes.swap(nodes);
        }
    }

    using gsQuadRule<T>::m_weights;
    using gsQuadRule<T>::m_nodes;

};

} // namespace gismo


