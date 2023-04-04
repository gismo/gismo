/** @file gsNewtonCotesRule.h

    @brief Provides the Newton-Cotes quadrature rule

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
    \brief Class that represents the (tensor) Newton-Cotes quadrature rule

    \ingroup Assembler
*/
template<class T>
class gsNewtonCotesRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsNewtonCotesRule> uPtr;

    /// Default empty constructor
    gsNewtonCotesRule() { }

    /// Initialize a tensor-product Newton-Cotes quadrature rule with \a numNodes (direction-wise)
    gsNewtonCotesRule(gsVector<index_t> const & numNodes,
                const unsigned digits = 0 )
    {
        gsNewtonCotesRule::setNodes(numNodes, digits);
    }

    /// Make function returning a smart pointer
    static uPtr make(gsVector<index_t> const & numNodes,
                        const unsigned digits = 0 )
    { return uPtr( new gsNewtonCotesRule(numNodes,digits) ); }

    /// Initialize a 1D Newton-Cotes quadrature rule with \a numNodes
    gsNewtonCotesRule(index_t numNodes, const unsigned digits = 0 )
    {
        this->setNodes(numNodes, digits);
    }

    /// Initialize a tensor-product Newton-Cotes quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsNewtonCotesRule(const gsBasis<T> & basis, const T quA, const index_t quB, short_t fixDir = -1);
    //const unsigned digits = std::numeric_limits<T>::digits10 );

    /// Initialize a tensor-product Newton-Cotes quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise). Values of quA
    /// and quB are taken from the \a options
    gsNewtonCotesRule(const gsBasis<T> & basis, const gsOptionList & options, short_t fixDir = -1);
    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsNewtonCotesRule() { }
    
public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes, 
                   unsigned digits = 0 );

    using gsQuadRule<T>::setNodes; // unhide base

private:

    void init(const gsBasis<T> & basis, const T quA,
              const index_t quB, short_t fixDir);
    
    /**
     * @brief Computes the Newton-Cotes quadrature rule with \a n
     * nodes in the interval [-1,1].
     *
     */
    static void computeReference(index_t n, gsVector<T> & x, gsVector<T> & w);

}; // class gsNewtonCotesRule


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNewtonCotesRule.hpp)
#endif
