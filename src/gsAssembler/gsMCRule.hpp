/** @file gsGaussRule.hpp

    @brief Provides implementation of the Gauss-Legendre quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsDomain/gsDomain.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template<class T> void
gsGaussRule<T>::init(const gsDomain<T> & domain, const T quA, const index_t quB, short_t fixDir)
//const unsigned digits)
{
    const short_t d  = domain.dim();
    GISMO_ASSERT( fixDir < d && fixDir>-2, "Invalid input fixDir = "<<fixDir);

    std::vector<gsVector<T> > nodes(d);
    std::vector<gsVector<T> > weights(d);
    if (-1==fixDir)
        fixDir = d;
    else
    {
        nodes  [fixDir].setZero(1); // numNodes == 1
        weights[fixDir].setConstant(1, 2.0);
    }

    // Note: skipping accuracy and lookup tests here (commented)

    //if (digits <= 30 )
    //{
    short_t i;
    for(i=0; i!=fixDir; ++i )
    {
        //note: +0.5 for rounding
        const index_t numNodes = cast<T,index_t>(quA * static_cast<T>(domain.degree(i)) + static_cast<T>(quB) + static_cast<T>(0.5));
        computeReference(numNodes, nodes[i], weights[i]);
    }
    ++i;// skip fixed direction
    for(; i<d; ++i )
    {
        const index_t numNodes = cast<T,index_t>(quA * static_cast<T>(domain.degree(i)) + static_cast<T>(quB) + static_cast<T>(0.5));
        computeReference(numNodes, nodes[i], weights[i]);
    }

    //}
    //else
    //{
    //    for( short_t i=0; i<d; ++i )
    //    {
    //        const index_t numNodes = quA * basis.degree(i) + quB;
    //        computeReference(numNodes, nodes[i], weights[i], digits);
    //    }
    //}

    this->computeTensorProductRule(nodes, weights);
}

template<class T>
gsGaussRule<T>::gsGaussRule(const gsDomain<T> & domain,
                            const T quA, const index_t quB,
                            const short_t fixDir)
//const unsigned digits)
{
    init(domain, quA, quB, fixDir);
}

template<class T>
gsGaussRule<T>::gsGaussRule(const gsDomain<T> & domain,
                            const gsOptionList & options,
                            const short_t fixDir)
{
    const T       quA = options.getReal("quA");
    const index_t quB = options.getInt ("quB");
    init(domain, quA, quB, fixDir);
}

template<class T>
gsGaussRule<T>::gsGaussRule(const gsBasis<T> & basis,
                            const T quA, const index_t quB,
                            const short_t fixDir)
:
gsGaussRule(*basis.domain(), quA, quB, fixDir)
{
}

template<class T>
gsGaussRule<T>::gsGaussRule(const gsBasis<T> & basis,
                            const gsOptionList & options,
                            const short_t fixDir)
:
gsGaussRule(*basis.domain(),options,fixDir)
{
}


template<class T> void
gsGaussRule<T>::setNodes( gsVector<index_t> const & numNodes,
                          unsigned digits)
{
    const index_t d = numNodes.rows();

    // Get base rule nodes and weights
    std::vector<gsVector<T> > nodes(d);
    std::vector<gsVector<T> > weights(d);

    if (digits == 0)
    {
        for (index_t i = 0; i < d; ++i)
        {
            computeReference(numNodes[i], nodes[i], weights[i], REAL_DIG);
        }
    }
    else
    {
        for (index_t i = 0; i < d; ++i)
            computeReference(numNodes[i], nodes[i], weights[i], digits);
    }

    this->computeTensorProductRule(nodes, weights);
}

template<class T> void
gsGaussRule<T>::computeReference(index_t n,       // Number of points
                                 gsVector<T> & x, // Quadrature points
                                 gsVector<T> & w, // Quadrature weights
                                 unsigned digits) // Number of exact decimal digits
{
    // Allocate space for points and weights.
    x.resize(n);
    w.resize(n);
}

} // namespace gismo
