/** @file gsAdditiveOp.hpp

    @brief Allows to set up additive Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

namespace gismo
{

template<typename T>
void gsAdditiveOp<T>::apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
{
    GISMO_ASSERT( this->rows() == input.rows(), "The dimensions do not fit." );

    x.setZero( input.rows(), input.cols() );

    const index_t n = m_ops.size();
    gsMatrix<T> res_local, corr_local;

    for (index_t i=0; i<n; ++i)
    {
        res_local.noalias() = m_transfers[i]->transpose()*input;
        m_ops[i]->apply(res_local, corr_local);
        x.noalias() += *(m_transfers[i])*corr_local;
    }
}

} // namespace gismo
