/** @file gsCrossApPatch.hpp

    @brief Provides cross approximation parameterizations from boundary data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsModeling/gsCrossApPatch.h>
#include<gsNurbs/gsTensorBSplineBasis.h>
#include<gsNurbs/gsBSpline.h>
#include<gsTensor/gsGridIterator.h>

namespace gismo
{

template <typename T>
const gsGeometry<T> & gsCrossApPatch<T>::compute()
{
    const short_t dim = m_boundary.dim();

    delete m_result;
    m_result = NULL;

    switch ( dim ) // dispatch to implementation
    {
    case 1:
        compute_impl<2>();
        break;
    /*
      case 2:
      compute_impl<3>();
      break;
      case 3:
      compute_impl<4>();
      break;
    */
    default:
        GISMO_ERROR("Dimension "<< dim << "is invalid.");
        break;
    }
    return *m_result;
}

template <typename T>
template <unsigned d>
void gsCrossApPatch<T>::compute_impl()
{
    GISMO_STATIC_ASSERT(d==2,"gsCrossApPatch<T>::compute_impl only works for d=2 dimensions.");

    gsTensorBSplineBasis<d,T> resultBasis; // Basis
    gsMatrix<T> coefs;                     // CPs

    // Resolve boundary configuration and set boundary coefficients
    this->preparePatch(resultBasis, coefs);

    gsVector<index_t,d> sz;
    // Check whether there are any interior points to fill in
    resultBasis.size_cwise(sz);
    if ( (sz.array() < 3).all() )
    {
        gsWarn<<"There where no interior control points.\n";
        m_result = resultBasis.makeGeometry( give(coefs) ).release();
        return;
    }

    gsMatrix<T> tmp0(sz[0],2), tmp1(sz[1],2), tmp;
    gsMatrix<T,d,d> cross;

    for (index_t i = 0; i!=coefs.cols(); ++i)
    {
        tmp0.col(0) = m_boundary[0].coefs().col(i);
        tmp0.col(1) = m_boundary[1].coefs().col(i);
        tmp1.col(0) = m_boundary[2].coefs().col(i);
        tmp1.col(1) = m_boundary[3].coefs().col(i);
        cross(0,0) = tmp0(0      ,0);
        cross(1,0) = tmp0(sz[0]-1,0);
        cross(0,1) = tmp0(0      ,1);
        cross(1,1) = tmp0(sz[0]-1,1);

        if ( math::abs(cross.determinant()) < 1e-11 )
        {
            gsWarn <<"Corner data is rank-deficient ("<<i<<")\n";
            //gsDebugVar(cross);
            tmp = tmp0 *
                cross.colPivHouseholderQr().solve(gsMatrix<T,2,2>::Identity())
                * tmp1.transpose();
        }
        else
        {
            tmp = tmp0 * cross.inverse() * tmp1.transpose();
        }

        coefs.col(i) = tmp.asVector();
    }

    // return the patch
    m_result = resultBasis.makeGeometry( give(coefs) ).release();
}


}// namespace gismo
