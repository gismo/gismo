/** @file gsNewtonCotesRule.hpp

    @brief Provides implementation of the Newton-Cotes quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template<class T> void
gsNewtonCotesRule<T>::init(const gsBasis<T> & basis, const T quA,
                           const index_t quB, short_t fixDir)
//const unsigned digits)
{
    const short_t d  = basis.dim();
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

    //if (digits <= 30 )
    //{
    short_t i;
    for(i=0; i!=fixDir; ++i )
    {
        //note: +0.5 for rounding
        const index_t numNodes = cast<T,index_t>(quA * static_cast<T>(basis.degree(i)) + static_cast<T>(quB) + static_cast<T>(0.5));
        //const bool found = 
        computeReference(numNodes, nodes[i], weights[i]);
    }
    ++i;// skip fixed direction
    for(; i<d; ++i )
    {
        const index_t numNodes = cast<T,index_t>(quA * static_cast<T>(basis.degree(i)) + static_cast<T>(quB) + static_cast<T>(0.5));
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
gsNewtonCotesRule<T>::gsNewtonCotesRule(const gsBasis<T> & basis, 
                            const T quA, const index_t quB,
                            const short_t fixDir)
//const unsigned digits)
{
    init(basis, quA, quB, fixDir);
}

template<class T>
gsNewtonCotesRule<T>::gsNewtonCotesRule(const gsBasis<T> & basis, 
                            const gsOptionList & options,
                            const short_t fixDir)
//const unsigned digits)
{
    const T       quA = options.getReal("quA");
    const index_t quB = options.getInt ("quB");
    init(basis, quA, quB, fixDir);
}


template<class T> void
gsNewtonCotesRule<T>::setNodes( gsVector<index_t> const & numNodes, 
                          unsigned digits)
{
    GISMO_UNUSED(digits);
    const index_t d = numNodes.rows();

    // Get base rule nodes and weights
    std::vector<gsVector<T> > nodes(d);
    std::vector<gsVector<T> > weights(d);

    for (index_t i = 0; i < d; ++i)
        computeReference(numNodes[i], nodes[i], weights[i]);

    this->computeTensorProductRule(nodes, weights);
}



template<class T> void
gsNewtonCotesRule<T>::computeReference(index_t n,       // Number of points
                                      gsVector<T> & x, // Quadrature points
                                      gsVector<T> & w) // Quadrature weights
{

    index_t i,j,k;
    x.resize(n);
    
    if ( 1 == n )
    {
        x[0] = (T)(0);
        w.setConstant(1,(T)(2));
        return;
    }

    //todo: if (closed) ... else (open)
    x[0]   = (T)(-1);
    x[n-1] = (T)( 1);
    for (i = 1; i < n-1; ++i)
        x[i] = ( T ) (2*i-n+1) / (T)(n-1);

    // todo: if weights_required
    // weights
    w.resize(n);
    gsVector<T> d(n);
    T a, b;
    for (i = 0; i < n; ++i)
    {
        //  Compute the Lagrange basis polynomial
        d.setZero();
        d[i] = (T)1.0;

        for (j = 2; j <= n; ++j)
            for (k = j; k <= n; ++k)
                d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );

        for ( j = 1; j <= n - 1; ++j)
            for ( k = 1; k <= n - j; ++k)
                d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];

        //  Evaluate the antiderivative of the polynomial at the endpoints
        a = b = d[n-1] / (T)( n );
        for ( j = n - 2; 0 <= j; --j)
        {
            const T dj1= d[j] / (T)(j+1);
            a+=      dj1;
            b = -b + dj1;
        }
        w[i] = a + b;
    }
    x[0] += 1e-12;
    x[n-1] -= 1e-12;
}

} // namespace gismo
