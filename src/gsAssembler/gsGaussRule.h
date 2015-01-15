/** @file gsGaussRule.h

    @brief Provides the Gauss-Legendre quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsQuadRule.h>

namespace gismo
{

/** 
    Class that represents the (tensor) Gauss-Legendre quadrature rule
    
    \ingroup Assembler
*/ 
 
template<class T>
class gsGaussRule : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsGaussRule() { }
    
    gsGaussRule(gsVector<index_t> const & numNodes, 
                unsigned digits = std::numeric_limits<T>::digits10 ): 
        gsQuadRule<T>() 
        { 
            setNodes(numNodes, digits);
        }

    gsGaussRule(index_t numNodes, 
                unsigned digits = std::numeric_limits<T>::digits10 ): 
        gsQuadRule<T>() 
        { 
            setNodes(numNodes, digits);
        }

    
    ~gsGaussRule() { }
    
public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes, 
                   unsigned digits = std::numeric_limits<T>::digits10 );

    void setNodes( index_t numNodes, 
                   unsigned digits = std::numeric_limits<T>::digits10 )
    {
        gsVector<index_t> nn(1);
        nn[0] = numNodes;
        setNodes(nn, digits);
    }
    
private:
    
    /**
     * @brief Computes the Gauss quadrature rule with \a n nodes in the interval [-1,1].
     *
     * This function is called by setNodes(), if lookupReference() (which is called first) returned \a false.
     */
    static void computeReference(index_t n, gsVector<T> & x, gsVector<T> & w, 
                                 unsigned digits =  std::numeric_limits<T>::digits10);

    /** 
     *@brief  Look up function for the Gauss quadrature rule in the interval [-1,1].
     *
     * When the member function setNodes() is called, it will first try to look up
     *the corresponding Gauss rule. If this look up was not successful, the function computeReference() will be called.
     *\return \a true if the look up was successful
     */
    static bool lookupReference (index_t n, gsVector<T> & x, gsVector<T> & w);

}; // class gsGaussRule


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGaussRule.hpp)
#endif
