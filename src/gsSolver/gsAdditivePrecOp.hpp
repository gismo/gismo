/** @file gsAdditivePrecOp.hpp

    @brief Allows to set up additive Schwarz type smoothers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

template<typename T>
void gsAdditivePrecOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == f.rows() && m_A->cols() == m_A->rows() && x.cols() == f.cols(),
        "The dimensions do not fit." );

    const index_t n = m_ops.size();
    m_A->apply( x, m_res );
    m_res -= f;

    for (index_t i=0; i<n; ++i)
    {
        m_res_local.noalias() = m_transfers[i].transpose()*m_res;
        m_ops[i]->apply(m_res_local, m_corr_local);
        m_corr_local *= m_damping;
        x.noalias() -= m_transfers[i]*m_corr_local;
    }
}

} // namespace gismo
