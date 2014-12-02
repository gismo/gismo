/** @file gsQuadRule.h

    @brief Provides a base class for a quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

/** 
    Class for a reference quadrature rule
*/
  
template<class T>
class gsQuadRule
{

public:

    /// Default empty constructor
    gsQuadRule() 
    { }

    // Construct a quadrature rule with \a numNodes number of quadrature
    // points per integration variable
    // gsQuadRule(gsVector<index_t> const & numNodes, 
    //            unsigned digits = std::numeric_limits<T>::digits10 )
    // { 
    //     setNodes(numNodes, digits);
    // }
    
    virtual ~gsQuadRule() { }

    // /// Setup a quadrature rule with \a numNodes number of quadrature
    // /// points per integration variable

    /**
     * @brief Initialize quadrature rule with \a numNodes number of
     * quadrature points per integration variable
     *
     * The call of this function initializes the quadrature rule for
     * further use, i.e., the quadrature points and weights on the
     * reference cube are set up and stored.  The dimension \a d of
     * the reference cube is specified by the length of the vector \a
     * numNodes.
     *
     * Example: numNodes = [2,5] corresponds to quadrature in 2D where
     * two quadrature points are used for the first coordinate and
     * five quadrature points for the second coordinate. In this case,
     * 2*5 = 10 quadrature nodes and weights are set up.
     *
     * @param numNodes vector of length \a d containing numbers of
     * nodes to be used (per integration variable).
     *
     * @param digits accuracy of nodes and weights.
     *
     */
    virtual void setNodes( gsVector<index_t> const & numNodes,
                           unsigned digits = std::numeric_limits<T>::digits10 )
    { GISMO_NO_IMPLEMENTATION }

    /**
     * @brief Returns reference nodes for the currently kept rule.
     *
     * @param[out] v Vector of length \a d = dim(), where each entry
     * \a v_i of \a v is again a vector. \a v_i contains the reference
     * quadrature points for the <em>i</em>-th coordinate direction.
     *
     */
    const gsMatrix<T> & referenceNodes() const { return m_nodes; }

    /**
     * @brief Returns reference weights for the currently kept rule.
     *
     * @param[out] v Vector of length \a d = dim(), where each entry
     * \a v_i of \a v is again a vector. \a v_i contains the reference
     * quadrature weights for the <em>i</em>-th coordinate direction.
     *
     */
    const gsVector<T> & referenceWeights() const { return m_weights; }

    // Reference element
    //const gsMatrix<T> & referenceElement() { }

    /// Number of nodes in the currently kept rule
    index_t numNodes() const { return m_weights.size(); }

    /// Dimension of the rule
    index_t dim() const { return m_nodes.rows(); }


    /// Maps quadrature rule (i.e., points and weights) from the reference interval to an element.
    /**
    * The currently kept rule (which is initialized on the reference hypercube [-1,1]^<em>d</em> by calling setNodes()) is mapped
    * to the <em>d</em>-dimensional hypercube specified by \a lower and \a upper.\n
    * For example, for <em>d=2</em>, the square <em>[a,b]x[c,d]</em> is defined by <em>lower = [a,c]</em>, <em>upper = [b,c]</em>.
    * \param[in] lower vector of length \a d, defining the coordinates of the lower corner of the hypercube.
    * \param[in] upper vector of length \a d, defining the coordinates of the upper corner of the hypercube.
    * \param[in,out] nodes will be overwritten with the coordinates of the quadrature nodes.\n
    * Size of the matrix \a nodes = <em>d</em> x <em>n</em>, where \n
    * \a d is the dimension of the element, and\n
    * \a n is the number of quadrature points.
    * \param[in,out] weights will be overwritten with the corresponding Gauss quadrature weights.\n
    * Length of the vector \a weights = number of quadrature nodes.
    */
    inline void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const;

    inline void mapTo( T startVal, T endVal,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const;
    
protected:

    /// Reference quadrature nodes (on the interval [-1,1]).
    gsMatrix<T> m_nodes;

    /// Reference quadrature weights (corresponding to interval [-1,1]).
    gsVector<T> m_weights;

}; // class gsQuadRule


//////////////////////////////////////////////////
//////////////////////////////////////////////////


template<class T> void
gsQuadRule<T>::mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                      gsMatrix<T> & nodes, gsVector<T> & weights ) const
{
    const index_t d = lower.size();
    GISMO_ASSERT( d == m_nodes.rows(), "Inconsistent quadrature mapping");
    
    nodes.resize( m_nodes.rows(), m_nodes.cols() );
    weights.resize( m_weights.size() );
    nodes.setZero();
    weights.setZero();

    gsVector<T> h(d);
    T hprod(1.0); // for the computation of the size of the cube.

    for ( index_t i = 0; i<d; ++i)
    {
        // the factor 0.5 is due to the fact that the one-dimensional reference interval is [-1,1].
        h[i] = ( lower[i] != upper[i] ? 0.5 * (upper[i]-lower[i]) : T(0.5) );
        hprod *= h[i];
    }
  
    // Linear map from [-1,1]^d to [lower,upper]
    nodes.noalias()   = ( h.asDiagonal() * m_nodes ).colwise() + 0.5*(lower+upper);

    // Adjust the weights (multiply by the Jacobian of the linear map)
    weights.noalias() = hprod * m_weights;
}

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


} // namespace gismo
