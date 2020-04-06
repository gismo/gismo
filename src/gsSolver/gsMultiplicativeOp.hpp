/** @file gsMultiplicativeOp.hpp

    @brief Allows to set up multiplicative Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Sogn
*/

namespace gismo
{

template<typename T>
void gsMultiplicativeOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying->rows() == x.rows() && x.rows() == f.rows() && m_underlying->cols() == m_underlying->rows() && x.cols() == f.cols(),
        "Dimensions do not match."<<m_underlying->rows()<<"=="<<x.rows()<<"&&"<<x.rows()<<"=="<<f.rows()<<"&&"<<m_underlying->cols()<<"=="<<m_underlying->rows()<<"&&"<<x.cols()<<"=="<<f.cols());

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    //Use optimized step function if available  
    if (m_optReady)
    {
        gsMultiplicativeOp<T>::stepOpt(f, x);
    }
    else
    {
        gsMatrix<T> p(f.rows(),f.cols());
        const size_t sz = m_transfers.size();
        for (size_t i = 0; i < sz; ++i)
        {
            m_ops[i]->apply( m_transfers[i].transpose() * (f - (*m_underlying)*x), p );        
            x += m_transfers[i] * p;
        }
    }
}



template<typename T>
void gsMultiplicativeOp<T>::stepT(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying->rows() == x.rows() && x.rows() == f.rows() && m_underlying->cols() == m_underlying->rows() && x.cols() == f.cols(),
        "Dimensions do not match.");

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    //Use optimized stepT function if available
    if (m_optReady)
    {
        gsMultiplicativeOp<T>::stepTOpt(f, x);
    }
    else
    {
        gsMatrix<T> p(f.rows(),f.cols());
        const size_t sz = m_transfers.size();
        for (size_t i = sz-1; i != (size_t)-1; --i)
        {
            m_ops[i]->apply( m_transfers[i].transpose() * (f - (*m_underlying)*x), p );        
            x += m_transfers[i] * p;
        }
    }    

}


template<typename T>
void gsMultiplicativeOp<T>::stepOpt(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    gsMatrix<T> p(f.rows(),f.cols());
    gsMatrix<T> rest_r(f.rows(),f.cols());
    gsMatrix<T> r(f.rows(),f.cols());     

    r = f - (*m_underlying)*x; //inital residual

    const size_t sz = m_transfersVec.size();
    for (size_t k = 0; k < sz; ++k)
    {
        const size_t szME = m_transfersVec[k].rows(); //Size macro element

        rest_r.setZero(szME,1);         //restricted residual
        for(size_t j=0; j < szME; ++j)
        {
            rest_r(j,0) = r(m_transfersVec[k](j),0); 
        }
            
        m_ops[k]->apply(rest_r, p);
                
        //Update residual and x
        for(size_t j=0; j < szME; ++j)
        {
            for(typename gsSparseMatrix<T, ColMajor>::InnerIterator it(m_underlyingColMajorMat,m_transfersVec[k](j)); it; ++it)
            {
                r(it.row(),0) -= it.value()*p(j,0);
            }
            x(m_transfersVec[k](j),0) += p(j,0);
        }
    }
}


template<typename T>
void gsMultiplicativeOp<T>::stepTOpt(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    gsMatrix<T> p(f.rows(),f.cols());
    gsMatrix<T> rest_r(f.rows(),f.cols());
    gsMatrix<T> r(f.rows(),f.cols());
    
    r = f - (*m_underlying)*x; //inital residual
    
    const size_t sz = m_transfersVec.size();
    
    for (size_t k = sz-1; k != (size_t)-1; --k)
    {
        const size_t szME = m_transfersVec[k].rows(); //Size macro element

        rest_r.setZero(szME,1);         //restricted residual
        for(size_t j=0; j < szME; ++j)
        {
            rest_r(j,0) = r(m_transfersVec[k](j),0); 
        }
            
        m_ops[k]->apply(rest_r, p);
                
        //Update residual and x
        for(size_t j=0; j < szME; ++j)
        {
            for(typename gsSparseMatrix<T, ColMajor>::InnerIterator it(m_underlyingColMajorMat,m_transfersVec[k](j)); it; ++it)
            {
                r(it.row(),0) -= it.value()*p(j,0);
            }
            x(m_transfersVec[k](j),0) += p(j,0);
        }
    }
}

} // namespace gismo
