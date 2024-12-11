/** @file gsLowRankCorrectedOp.hpp

    @brief Provides the inverse of \f$Ainv^{-1} + U Q^{-1} V^T\f$ with \f$U Q^{-1} V^T\f$ being a low rank matrix using the Sherman Morrisson Woodburry formula

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#include <gsSolver/gsLowRankCorrectedOp.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

template<typename T>
gsLowRankCorrectedOp<T>::gsLowRankCorrectedOp(const BasePtr& Ainv, const gsSparseMatrix<T> & Q, const gsSparseMatrix<T> & U, const gsSparseMatrix<T> & V)
    : m_U( U ), m_V( V )
{
    GISMO_ASSERT
    (
        Ainv->cols() == Ainv->rows()
        && Ainv->rows() == U.rows()
        && U.cols() == Q.cols()
        && Q.cols() == Q.rows()
        && Q.rows() == V.cols()
        && V.rows() == Ainv->cols(),
        "Dimensions do not fit."
    );

    m_Ainv = Ainv;
    
    // Here, we implement the following:
    // W = Q - V.transpose() * inv(A) * U
    // tmp = inv(A) * U
    gsMatrix<T> tmp( Ainv->rows(), m_U.cols() );
    m_Ainv->apply( m_U, tmp );
    
    // W = Q - V.transpose() * tmp
    gsMatrix<T> W( Q );
    W.noalias() -= m_V.transpose() * tmp;

    m_Winv = makePartialPivLUSolver( W ); 

}
    
template<typename T>
void gsLowRankCorrectedOp<T>::apply(const gsMatrix<T>& x, gsMatrix<T>& result) const
{
    // Here, we could make permanently allocated vectors
    gsMatrix<T> tmp1;
    gsMatrix<T> tmp2;
    gsMatrix<T> tmp3;
    
    // compute    Ainv * ( I + U * Winv * V.transpose() * Ainv ) * x
 
    m_Ainv->apply( x, result );

    if( m_U.cols() == 0 )
        return;

    tmp1.noalias() = m_V.transpose() * result;
    m_Winv->apply( tmp1, tmp2 );
    tmp3.noalias() = m_U * tmp2;
    tmp3 += x;
    m_Ainv->apply( tmp3, result );

}


}

