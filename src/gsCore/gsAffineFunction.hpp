/**  @file gsAffineFunction.hpp

    @brief Implements an affine function.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on:  2014-11-27
*/

namespace gismo {

template <typename T>
gsAffineFunction<T>::gsAffineFunction(const gsVector<index_t> &dir, const gsVector<bool> &o, const gsMatrix<T> &box1, const gsMatrix<T> &box2)
{
    GISMO_ASSERT(box1.rows()==box2.rows(),
                 "The two boxes must be subset of Rn for the same n (same number of rows)");
    GISMO_ASSERT(box1.cols()==2 && 2==box2.cols(),
                 "The two boxes must be described by the lower and upper corner, the matrices must have two columns");

    const index_t dim = box1.rows();
    const gsVector<T> size1 = box1.col(1) - box1.col(0);
    const gsVector<T> size2 = box2.col(1) - box2.col(0);
    m_mat.setZero(dim,dim);
    m_trans.resize(dim);
    for (index_t i=0; i<dim; ++i)
    {
        const T ratio = size1[i]==0 ? T(1) : size2(dir[i])/size1[i];
        m_mat(dir(i),i) = o[i] ? ratio : -ratio;
        m_trans(dir(i)) = o[i] ? box2(dir[i],0) : box2(dir[i],1);
    }
    m_trans -= m_mat * box1.col(0);
}


template <typename T>
gsAffineFunction<T>::gsAffineFunction(const gsMatrix<T> mat, const gsVector<T> trans)
    : m_mat(mat), m_trans(trans)
{
    GISMO_ASSERT(m_mat.rows()==m_trans.rows(),
                 "Incompatible linear map and translation in affine map");
}

template <typename T>
short_t gsAffineFunction<T>::domainDim() const
{
    return m_mat.cols();
}

template <typename T>
short_t gsAffineFunction<T>::targetDim() const
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
