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

namespace gismo
{

/// Helper function to compute a lexicographically numbered index from tensor indices.
template <unsigned d>
int fromTensorIndex(const gsVector<unsigned, d>& idx, const gsVector<unsigned, d>& sz)
{
    int result = 0;
    unsigned stride = 1;
    for (unsigned i = 0; i < d; ++i)
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
template <unsigned d, typename T>
void tensorCombineTransferMatrices(
    gsSparseMatrix<T,RowMajor> B[d],
    gsSparseMatrix<T,RowMajor> & transfer)
{
    typedef typename gsSparseMatrix<T,RowMajor>::InnerIterator InnerIt;
    gsSparseEntries<T> entries;
    gsVector<unsigned, d> oldSize, newSize, v, v_old;

    for (unsigned i = 0; i < d; ++i)
    {
        oldSize[i] = (unsigned) B[i].innerSize();
        newSize[i] = (unsigned) B[i].outerSize();
    }

    // iterate over all new tensor indices
    v.setZero();
    index_t newIdx = 0;
    do {
        assert( newIdx == fromTensorIndex<d>(v, newSize) );

        // set up iterators over component contributions
        //InnerIt it[d];
        std::vector<InnerIt> it;
        it.reserve(d);
        for (unsigned i = 0; i < d; ++i)
            it.push_back(  InnerIt( B[i], v[i] ) );

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
            index_t oldIdx = fromTensorIndex<d>(v_old, oldSize);

            entries.add( newIdx, oldIdx, contrib );

            // advance iterators
            more = true;
            for (unsigned i = 0; i < d; ++i)
            {
                ++it[i];                    // increase current dimension
                if (it[i])                  // current dimension not yet exhausted?
                    break;
                else                        // current dimension exhausted
                {
                    if (i == d - 1)         // was it the last one?
                        more = false;       // then all elements exhausted
                    else
                        it[i] = InnerIt( B[i], v[i] ); // otherwise, reset this to start and increase the next dimension
                }
            }
        } while (more);
    } while (++newIdx, nextLexicographic(v, newSize));

    assert( newIdx == (index_t) newSize.prod() );

    transfer.resize( newSize.prod(), oldSize.prod() );
    transfer.setFrom( entries );
    transfer.makeCompressed();
}

/// Helper to compute the strides of a d-tensor
template<int d>
void tensorStrides(const gsVector<int,d> & sz, gsVector<int,d> & strides) 
{
    strides.resize(sz.size());
    strides[0] = 1;
    for ( index_t i=1; i< sz.size(); ++i )
        strides[i] = strides[i-1] * sz[i-1];
}


/// Reorders (inplace) the given tensor \a coefs vector (regarded as a
/// \a sz.prod() x \a d matrix arranged as a flattened \a sz tensor,
/// so that the rows are re-arranged according to the input
/// permutation \a perm. The \a sz is updated to the new ordering.
template <typename T, int d>
void permuteTensorVector( const gsVector<int,d> & perm, 
                          gsVector<int,d> & sz, 
                          gsMatrix<T> & coefs)
{
    const int dd = sz.size();

    GISMO_ASSERT( perm.sum() == dd*(dd-1)/2, "Error in the permutation: "<< perm.transpose());

    Eigen::PermutationMatrix<d> P(perm);
    const gsVector<int,d> oldsz(sz);

    gsVector<int,d> perstr(dd);
    tensorStrides<d>(P*sz, perstr);

    // check: is it better to create a big permutation to apply to coefs ?
    //        otherwise, is the swapping possible without the temporary ?
    gsMatrix<T> tmp(coefs.rows(), coefs.cols() );

    index_t r = 0;
    gsVector<int,d> v(dd);
    v.setZero();
    do 
    {
        tmp.row( perstr.dot(P*v) ) = coefs.row( r++ );
    } while (nextLexicographic(v, sz));

    coefs.swap(tmp);
    sz = P * sz;
}

} // namespace gismo
