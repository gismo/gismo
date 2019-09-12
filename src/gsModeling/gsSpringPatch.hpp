/** @file gsSpringPatch.hpp

    @brief Provides spring patch construction from boundary data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsModeling/gsSpringPatch.h>
#include<gsNurbs/gsTensorBSplineBasis.h>
#include<gsNurbs/gsBSpline.h>
#include<gsTensor/gsGridIterator.h>

namespace gismo
{
template <typename T>
const gsGeometry<T> & gsSpringPatch<T>::compute()
{
    const short_t dim = m_boundary.dim();

    delete m_result;
    m_result = NULL;

    switch ( dim ) // dispatch to implementation
    {
    case 1:
        compute_impl<2>();
        break;
    case 2:
        compute_impl<3>();
        break;
    case 3:
        compute_impl<4>();
        break;
    default:
        GISMO_ERROR("Dimension "<< dim << "is invalid.");
        break;
    }
    return *m_result;
}

template <typename T>
template <unsigned d>
void gsSpringPatch<T>::compute_impl()
{
    gsTensorBSplineBasis<d,T> resultBasis; // Basis for the Coon's patch
    gsMatrix<T> coefs;                     // CPs of the Coon's patch

    // Resolve boundary configuration and set boundary coefficients
    this->preparePatch(resultBasis, coefs);

    gsVector<index_t,d> stride; // tensor strides

    // Check whether there are any interior points to fill in
    resultBasis.size_cwise(stride);
    if ( (stride.array() < 3).all() )
    {
        gsWarn<<"There where no interior control points.\n";
        m_result = resultBasis.makeGeometry( give(coefs) ).release();
        return;
    }

    // Compute the tensor strides
    resultBasis.stride_cwise(stride);

    // Initialize mapper
    const int sz  = resultBasis.size();
    gsDofMapper mapper(resultBasis);
    // boundary control point indices
    gsMatrix<index_t> bnd = resultBasis.allBoundary();
    mapper.markBoundary(0,bnd);
    mapper.finalize();

    // Sparse system
    gsSparseMatrix<T> A(mapper.freeSize(), mapper.freeSize() );
    gsMatrix<T> b(mapper.freeSize(), coefs.cols() );
    A.reserve( gsVector<int>::Constant(A.cols(), 2*d) );
    A.setZero();
    b.setZero();

    const T dd = 2*d;

    // Fill in system matrix A and right-hand side b
    // 2 d c_i - sum neib(c_i) = 0
    for ( int i = 0; i<sz ; i++ )
    {
        const index_t ii = mapper.index(i,0);// get i-index in the matrix

        if ( !mapper.is_free_index(ii))
            continue;

        A(ii,ii) = dd;

        for ( unsigned k = 0; k<d; k++ ) // for all neighbors
        {
            for ( int s = -1; s<2; s+=2 ) // +/- (up or down)
            {
                const unsigned j = i + s * stride[k];
                const index_t jj = mapper.index(j,0);// get j-index in the matrix

                if ( mapper.is_free_index(jj) ) // interior node ?
                    A(ii,jj) = -1;
                else // boundary node
                    b.row(ii) += coefs.row(j);
            }
        }
    }

    A.makeCompressed();

    // Solve system
    typename gsSparseSolver<T>::QR  solver(A);
    gsMatrix<T> solution = solver.solve ( b );

    // Fill in interior coefficients
    for ( int i = 0; i<sz; i++ )
    {
        const index_t ii = mapper.index(i,0);

        if ( mapper.is_free_index(ii) ) // interior node?
        {
            coefs.row(i) = solution.row(ii);
        }
    }

    // return the spring patch
    m_result = resultBasis.makeGeometry( give(coefs) ).release();
}


}// namespace gismo
