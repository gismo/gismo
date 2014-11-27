/**  gsAffineFunction.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on:  2014-11-27
*/


#include  <gsCore/gsAffineFunction.h>


namespace gismo {



template <typename T>
gsAffineFunction<T>::gsAffineFunction(const gsMatrix<T> mat, const gsVector<T> trans)
    : m_mat(mat), m_trans(trans)
{
    GISMO_ASSERT(m_mat.rows()==m_trans.rows(),"INCOMPATIBLE LINEAR MAP AND TRANSLATION IN AFFINE MAP");
}

template <typename T>
int gsAffineFunction<T>::domainDim() const
{
    return m_mat.cols();
}

template <typename T>
int gsAffineFunction<T>::targetDim() const
{
    return m_mat.rows();
}

template <typename T>
void gsAffineFunction<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result = (m_mat*u).colwise()+m_trans;
}

template <typename T>
void gsAffineFunction<T>::eval_component_into(const gsMatrix<T>& u,
                                                      const index_t comp,
                                                      gsMatrix<T>& result) const
{
    gsMatrix<T> temp;
    eval_into(u,temp);
    result=temp.row(comp);
}

template <typename T>
void gsAffineFunction<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    result.resize(m_mat.rows()* m_mat.cols(), u.cols());
    if (!u.cols())
    {
        return;
    }
    for (index_t col=0;col<m_mat.cols();++col)
    {
        result.block(m_mat.rows()*col,0,m_mat.rows(),1)=m_mat.col(col);
    }
    for (index_t col=1; col<result.cols();++col)
    {
        result.col(col)=result.col(0);
    }
}

template <typename T>
void gsAffineFunction<T>::deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const
{
    const index_t rows = domainDim()*domainDim()*targetDim();
    result.setZero(rows,u.cols());
}



} // namespace gismo
