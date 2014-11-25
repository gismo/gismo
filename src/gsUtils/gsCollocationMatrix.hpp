/** @file gsCollocationMatrix.hpp

    @brief Provides implementation of routines to compute the collocation matrix.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>


namespace gismo
{

template<class T> inline
void gsCollocationMatrix_into (gsBasis<T> const& basis, gsMatrix<T> const& u, 
                       gsSparseMatrix<T> & res)
{
    res.resize( u.cols(), basis.size() );
    res.setZero();

    // Evaluate basis functions on u
    gsMatrix<T> ev;
    basis.eval_into(u, ev);
    // Get indices of nonzero functions
    gsMatrix<unsigned> act;
    basis.active_into(u, act);

    gsSparseEntries<T> entries;
    entries.reserve( ev.cols() * act.rows() );

    //Construct matrix :  
    //rows= samples 1,..,n - cols= basis functions 1,..,n
    for (index_t k=0; k!= ev.cols(); ++k)
        for (index_t i=0; i!=act.rows(); ++i)
            entries.add(k , act(i,k), ev(i,k));

    res.setFrom(entries);
    res.makeCompressed();
}

template<class T> inline
gsSparseMatrix<T> * gsCollocationMatrix (gsBasis<T> const& basis, gsMatrix<T> const& u)
{
    gsSparseMatrix<T> * res = new gsSparseMatrix<T>;
    gsCollocationMatrix_into<T>(basis,u,*res);
    return res;
}


}; // namespace gismo
