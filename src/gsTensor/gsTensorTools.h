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

///
/// @brief Helper function to compute a lexicographically numbered index from tensor indices.
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
/// so that the rows are re-arranged so that \a k1 and \a k2 are swapped
/// permutation \a perm. The \a sz is updated to the new ordering.
/// \ingroup Tensor
template <typename T, int d>
void swapTensorDirection( int k1, int k2,
                          gsVector<int,d> & sz, 
                          gsMatrix<T> & coefs)
{
    const int dd = sz.size();
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

    gsMatrix<T> tmp(coefs.rows(), coefs.cols() );
    gsVector<int,d> perstr;
    std::swap( sz[k1], sz[k2] );
    tensorStrides<d>(sz, perstr);
    std::swap( sz[k1], sz[k2] );
    
    index_t r = 0;
    gsVector<int,d> v(dd);
    v.setZero();
    do 
    {
        std::swap( v[k1], v[k2] );
        tmp.row(perstr.dot(v)) = coefs.row(r++);
        std::swap( v[k1], v[k2] );
    } 
    while (nextLexicographic(v, sz));

    coefs.swap(tmp);
    std::swap( sz[k1], sz[k2] );
}


/// Reorders (inplace) the given tensor \a coefs vector (regarded as a
/// \a sz.prod() x \a d matrix arranged as a flattened \a sz tensor,
/// so that the rows are re-arranged according to the input
/// permutation \a perm. The \a sz is updated to the new ordering.
/// \ingroup Tensor
template <typename T, int d>
void permuteTensorVector( const gsVector<index_t,d> & perm, 
                          gsVector<int,d> & sz, 
                          gsMatrix<T> & coefs)
{
    const int dd = sz.size();
    GISMO_ASSERT( sz.prod()  == coefs.rows(), 
                  "Input error, sizes do not match: "<<sz.prod()<<"!="<< coefs.rows() );
    GISMO_ASSERT( perm.sum() == dd*(dd-1)/2, "Error in the permutation: "<< perm.transpose());

    Eigen::PermutationMatrix<d> P(perm);

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
    } 
    while (nextLexicographic(v, sz));

    coefs.swap(tmp);
    sz = P * sz;
}

/// \brief Flips tensor directions in place
/// \ingroup Tensor
template <typename T, int d>
void flipTensorVector(const int dir,
                      const gsVector<int,d> & sz, 
                      gsMatrix<T> & coefs)
{
    GISMO_ASSERT( sz.prod()  == coefs.rows(), 
                  "Input error, sizes do not match: "<<sz.prod()<<"!="<< coefs.rows() );

    const int dd = sz.size();

    // compute strides
    gsVector<int,d> perstr(dd);
    tensorStrides<d>(sz, perstr);

    gsVector<int,d> v(dd), vend(dd);
    const int cc = sz[dir] - 1; 
    vend = sz; vend[dir] /= 2;
    v.setZero();
    do
    {
        const int i1 = perstr.dot(v);
        const int i2 = i1 + (cc-2*v[dir])*perstr[dir];

        coefs.row( i1 ).swap( coefs.row( i2 ) );
    } 
    while (nextLexicographic(v, vend));
}

/** \brief Computes the sparse Kronecker product of sparse matrix blocks.

    The sparse matrices \a m1 and \a m2 must have sizers n1 x k*n1 and
    n2 x k*n2 respectively.
    
    Let \f$ c_{1,k},\, c_{2,k}\f$ be the two blocks of m1 and m2, the result is 

    \f$ \sum_k c_{1,k} \prod c_{2,k} \f$

    \param[in]  m1     
    \param[in]  m2
    \param[out] result
    \param[in] nzPerCol

    \ingroup Tensor
*/
template<class T>
void gsSparseKroneckerProduct(const gsSparseMatrix<T> & m1, const gsSparseMatrix<T> & m2,
                              gsSparseMatrix<T> & result, index_t nzPerCol = 10)
{
    typedef typename gsSparseMatrix<T>::InnerIterator cIter;
    typedef typename std::vector<cIter>::iterator     vIter;

    // Assumes square coordinate matrices
    const index_t s2 = m2.rows(),
                  s1 = m1.rows();
    const index_t rk = m1.cols() / s1;

    result.resize (s1*s2, s1*s2);
    result.reserve(gsVector<index_t>::Constant(result.cols(), (nzPerCol+1)/2) );

    std::vector<cIter> it1(rk, cIter(m1,0) ), 
                       it2(rk, cIter(m2,0) );

    for (index_t k1=0; k1 != s1; ++k1) // for all cols of m1
        for (index_t k2=0; k2 != s2; ++k2) // for all cols of m2
        {
            for (index_t i=0; i != rk; ++i)
                it1[i] = cIter(m1, i*s1 + k1);
            
            for (; it1[0];) // for all rows of m1
            {
                for (index_t i=0; i != rk; ++i)
                    it2[i] = cIter(m2, i*s2 + k2);

                for (; it2[0];) // for all rows of m2
                {
                    const index_t i = it2[0].index() * s1 + it1[0].index(),
                                  j = k2             * s1 + k1   ;

                    // Lower triangular part only ?
                    //if ( j <= i )
                    {
                        T tmp = it1[0].value() * it2[0].value();
                        for (index_t r=1; r < rk; ++r)
                            tmp += it1[r].value() * it2[r].value();
                        result.insert(i,j) = tmp;
                    }

                    for (vIter i=it2.begin(); i != it2.end(); ++i) ++(*i);
                }

                for (vIter i=it1.begin(); i != it1.end(); ++i) ++(*i);
            }
        }

/*  // Equivalent Dense matrix version:
    T tmp;
    for (index_t c = 0; c != s1; c++) // for all cols of m1
        for (index_t j = 0; j != s2; j++) // for all cols of m2
            for (index_t r = 0; r != s1; r++)  // for all rows of m1
                for (index_t i = 0; i != s2; i++) // for all rows of m2
                {
                    tmp = m1(r,c)*m2(i,j);
                    for (index_t t = 1; t != rk; ++t)
                        tmp +=  m1(r,t*s1 + c)*m2(i,t*s2+j);
                    result(i*s1+r, j*s1+c) = tmp;
                }
// */
    //result.makeCompressed();
}


} // namespace gismo
