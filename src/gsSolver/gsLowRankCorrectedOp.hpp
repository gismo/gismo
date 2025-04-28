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

namespace {
template<typename T>
gsMatrix<T> operator*( const typename gsLinearOperator<T>::Ptr& a, const gsMatrix<T>& b )
{
    gsMatrix<T> tmp( a->rows(), b.cols() );
    a->apply(b, tmp);
    return tmp;
}
}
    
template<typename T>
gsLowRankCorrectedOp<T>::gsLowRankCorrectedOp(const BasePtr& Ainv, const gsSparseMatrix<T> & Q, const gsSparseMatrix<T> & U, const gsSparseMatrix<T> & V)
    : m_U(U), m_V(V), m_Ainv(Ainv)
{
    GISMO_ASSERT( Ainv->cols() == Ainv->rows(), "Not quadratic.");
    GISMO_ASSERT( Q.cols() == Q.rows(), "Not quadratic.");
    GISMO_ASSERT( Ainv->rows() == U.rows(), "Dimensions do not fit: "<<Ainv->rows()<<"=="<<U.rows());
    GISMO_ASSERT( U.cols() == Q.cols(), "Dimensions do not fit: "<<U.cols()<<"=="<<Q.cols());
    GISMO_ASSERT( Q.rows() == V.cols(), "Dimensions do not fit: "<<Q.rows()<<"=="<<V.cols());
    GISMO_ASSERT( V.rows() == Ainv->cols(), "Dimensions do not fit: "<<V.rows()<<"=="<<Ainv->cols());
   
    gsMatrix<T> W = Q + m_V.transpose() * (m_Ainv * gsMatrix<T>(m_U));
    m_Winv = makePartialPivLUSolver( W ); 

}
    
template<typename T>
void gsLowRankCorrectedOp<T>::apply(const gsMatrix<T>& x, gsMatrix<T>& result) const
{
   
    // compute    ( I - Ainv * U * Winv * V.transpose() ) * Ainv * x
 
    m_Ainv->apply( x, result );

    if( m_U.cols() == 0 )
        return;

    result -= m_Ainv * gsMatrix<T>( m_U * (m_Winv * gsMatrix<T>(m_V.transpose() * result)) );

}


}

