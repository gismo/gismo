/** @file gsLobattoRule.h

    @brief Provides the Gauss-Lobatto quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/**
    \brief Class that represents the (tensor) Gauss-Lobatto quadrature rule

    \ingroup Assembler
*/

template<class T>
class gsLobattoRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsLobattoRule> uPtr;

    /// Default empty constructor
    gsLobattoRule() { }

    /// Initialize a tensor-product Lobatto quadrature rule with \a numNodes (direction-wise)
    gsLobattoRule(gsVector<index_t> const & numNodes,
                const unsigned digits = 0 )
    {
        gsLobattoRule::setNodes(numNodes, digits);
    }

    /// Make function returning a smart pointer
    static uPtr make(gsVector<index_t> const & numNodes,
                        const unsigned digits = 0 )
    { return uPtr( new gsLobattoRule(numNodes,digits) ); }

    /// Initialize a 1D Lobatto quadrature rule with \a numNodes
    gsLobattoRule(index_t numNodes, const unsigned digits = 0 )
    {
        this->setNodes(numNodes, digits);
    }

    ~gsLobattoRule() { }

public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes,
                   unsigned digits = 0 );

    using gsQuadRule<T>::setNodes;// unhide base

private:

    /**
     * @brief Computes the Lobatto quadrature rule with \a n nodes in the interval [-1,1].
     *
     * This function is called by setNodes(), if lookupReference() (which is called first) returned \a false.
     */
    static void computeReference(index_t n, gsVector<T> & x, gsVector<T> & w,
                                 unsigned digits = 0 );

    /**
     *@brief  Look up function for the Lobatto quadrature rule in the interval [-1,1].
     *
     * When the member function setNodes() is called, it will first try to look up
     *the corresponding Lobatto rule. If this look up was not successful, the function computeReference() will be called.
     *\return \a true if the look up was successful
     */
    static bool lookupReference (index_t n, gsVector<T> & x, gsVector<T> & w);

}; // class gsLobattoRule


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsLobattoRule.hpp)
#endif
