/** @file gsTensorTools.h

    @brief Utility functions related to tensor-structured objects

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <gsUtils/gsCombinatorics.h>
#include <gsTensor/gsGridIterator.h>

namespace gismo
{

///
/// @brief Helper function to compute a lexicographically numbered index from tensor indices.
template <int d>
int fromTensorIndex(const gsVector<unsigned, d>& idx, const gsVector<unsigned, d>& sz)
{
    int result = 0;
    unsigned stride = 1;
    for (index_t i = 0; i < idx.size(); ++i)
    {
        result += stride * idx[i];
        stride *= sz[i];
    }
    return result;
}

/** @brief Combine component-wise transfer matrices into a transfer matrix for the tensor product basis.

    Given some kind of transformation (e.g., knot insertion/refinement) which transforms each component
    basis separately, this function computes a joint transfer matrix \a transfer which describes the
    transfer on the whole tensor product basis.

    The component transformations are allowed to change the size of the basis.
    
    \ingroup Tensor
 */
template <short_t d, typename T>
void tensorCombineTransferMatrices(
    gsSparseMatrix<T,RowMajor> B[d],
    gsSparseMatrix<T,RowMajor> & transfer)
{
    typedef typename gsSparseMatrix<T,RowMajor>::iterator InnerIt;
    gsSparseEntries<T> entries;
    gsVector<unsigned, d> oldSize, newSize, v, v_old;

    for (unsigned i = 0; i < d; ++i)
    {
        oldSize[i] = static_cast<unsigned>(B[i].innerSize());
        newSize[i] = static_cast<unsigned>(B[i].outerSize());
    }

    std::vector<InnerIt> it(d);

    // iterate over all new tensor indices
    v.setZero();
    index_t newIdx = 0;
    do {
        // set up iterators over component contributions
        for (unsigned i = 0; i < d; ++i)
            it[i] = B[i].begin(v[i]);

        // iterate over the component contributions
        bool more;
        do {
            // get the old tensor index
            T contrib = T(1);
            for (unsigned i = 0; i < d; ++i)
            {
                contrib *= it[i].value();
                v_old[i] = it[i].index();
            }

            const index_t oldIdx = fromTensorIndex<d>(v_old, oldSize);
            entries.add( newIdx, oldIdx, contrib );

            // advance iterators
            more = true;
            for (unsigned i = 0; i < d; ++i)
            {
                ++it[i];                  // increase current dimension
                if (it[i])                  // current dimension not yet exhausted?
                    break;
                else                        // current dimension exhausted
                {
                    if (i == d - 1)         // was it the last one?
                        more = false;       // then all elements exhausted
                    else
                        it[i] = B[i].begin(v[i]); // otherwise, reset this to start and increase the next dimension
                }
            }
        } while (more);
    } while (++newIdx, nextLexicographic(v, newSize));

    GISMO_ASSERT( newIdx == (index_t) newSize.prod(), "Iteration did not complete as expected." );

    transfer.resize( newSize.prod(), oldSize.prod() );
    transfer.setFrom( entries );
    transfer.makeCompressed();
}

/// \brief Helper to compute the strides of a d-tensor
template<typename VectIn, typename VectOut> inline
void tensorStrides(const VectIn & sz, VectOut & strides) 
{
    strides.derived().resize(sz.size());
    strides[0] = 1;
    for ( index_t i=1; i != sz.size(); ++i )
        strides[i] = strides[i-1] * sz[i-1];
}

/// Reorders (inplace) the given tensor \a coefs vector (regarded as a
/// \a sz.prod() x \a d matrix arranged as a flattened \a sz tensor,
/// so that the rows are re-arranged so that \a k1 and \a k2 are swapped
/// permutation \a perm. The \a sz is updated to the new ordering.
/// \ingroup Tensor
template <typename T, int d>
void swapTensorDirection( int k1, int k2,
                          gsVector<index_t,d> & sz, 
                          gsMatrix<T> & coefs)
{
    GISMO_ASSERT( sz.prod()  == coefs.rows(), 
                  "Input error, sizes do not match: "<<sz.prod()<<"!="<< coefs.rows() );
    GISMO_ASSERT( k1<d && k2 < d && k1>=0 && k2>=0,
                  "Invalid directions: "<< k1 <<", "<< k2 );
    
    if ( k1 == k2 )
        return; //Nothing to do
    
    /*
    gsVector<int,d> perm = gsVector<int,d>::LinSpaced(d,0,d-1);
    std::swap(perm[k1],perm[k2] );
    permuteTensorVector(perm,sz,coefs);
    return;
    */

    gsGridIterator<index_t,CUBE,d> it(sz);
    
    std::swap( sz[k1], sz[k2] );
    gsVector<index_t,d> perstr;
    tensorStrides(sz, perstr);
    std::swap(perstr[k1], perstr[k2] );
    
    gsMatrix<T> tmp(coefs.rows(), coefs.cols() );

    for(index_t r=0; it; ++it, ++r)
        tmp.row(perstr.dot(*it)) = coefs.row(r);

    coefs.swap(tmp);
}

/// Reorders (inplace) the given tensor \a coefs vector (regarded as a
/// \a sz.prod() x \a d matrix arranged as a flattened \a sz tensor,
/// so that the rows are re-arranged according to the input
/// permutation \a perm. The \a sz is updated to the new ordering.
/// \ingroup Tensor
template <typename T, int d>
void permuteTensorVector( const gsVector<index_t,d> & perm, 
                          gsVector<index_t,d> & sz, 
                          gsMatrix<T> & coefs)
{
    GISMO_ASSERT( sz.prod()  == coefs.rows(), 
                  "Input error, sizes do not match: "<<sz.prod()<<"!="<< coefs.rows() );
    GISMO_ASSERT( perm.sum() == sz.size()*(sz.size()-1)/2,
                  "Error in the permutation: "<< perm.transpose());

    if ( perm == gsVector<index_t>::LinSpaced(sz.size(),0,sz.size()-1) )
        return; //Nothing to do
        
    typename gsVector<index_t>::PermutationWrap P(perm);
    gsGridIterator<index_t,CUBE,d> it(sz);

    sz = P * sz;
    gsVector<index_t,d> perstr;
    tensorStrides(sz, perstr);
    perstr = P * perstr;
    
    // check: is it better to create a big permutation to apply to coefs ?
    gsMatrix<T> tmp(coefs.rows(), coefs.cols() );

    for(index_t r=0; it; ++it, ++r)
        tmp.row(perstr.dot(*it)) = coefs.row(r);
    
    coefs.swap(tmp);
}

/// \brief Flips tensor directions in place
/// \ingroup Tensor
template <typename T, int d>
void flipTensorVector(const int dir,
                      const gsVector<index_t,d> & sz, 
                      gsMatrix<T> & coefs)
{
    GISMO_ASSERT( sz.prod()  == coefs.rows(), 
                  "Input error, sizes do not match: "<<sz.prod()<<"!="<< coefs.rows() );

    gsVector<index_t,d> perstr = sz;
    perstr[dir] /= 2;
    gsGridIterator<index_t,CUBE,d> it(perstr);
    tensorStrides(sz, perstr);//reuse
    const index_t cc = sz[dir] - 1; 

    for(; it; ++it)
    {
        const index_t i1 = perstr.dot(*it);
        const index_t i2 = i1 + (cc - 2 * it->at(dir)) * perstr[dir];
        coefs.row( i1 ).swap( coefs.row( i2 ) );
    } 
}

} // namespace gismo
