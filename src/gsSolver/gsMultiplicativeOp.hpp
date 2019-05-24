/** @file gsMultiplicativeOp.hpp

    @brief Allows to set up multiplicative Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

namespace gismo
{

template<typename T>
void gsMultiplicativeOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying.rows() == x.rows() && x.rows() == f.rows() && m_underlying.cols() == m_underlying.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    // TODO: This is totally inefficient
    gsMatrix<T> p;
    const size_t sz = m_transfers.size();
    for (size_t i = 0; i < sz; ++i)
    {
        m_ops[i]->apply( m_transfers[i].transpose() * (f - m_underlying*x), p );        
        x += m_transfers[i] * p;
    }
}

template<typename T>
void gsMultiplicativeOp<T>::stepT(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying.rows() == x.rows() && x.rows() == f.rows() && m_underlying.cols() == m_underlying.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    // TODO: This is totally inefficient
    gsMatrix<T> p;
    const size_t sz = m_transfers.size();
    for (size_t i = sz-1; i != (size_t)-1; --i)
    {
        m_ops[i]->apply( m_transfers[i].transpose() * (f - m_underlying*x), p );        
        x += m_transfers[i] * p;
    }
}

} // namespace gismo
