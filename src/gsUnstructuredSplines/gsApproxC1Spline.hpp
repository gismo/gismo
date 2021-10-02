/** @file gsApproxC1Spline.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

namespace gismo
{

template<short_t d,class T>
void gsApproxC1Spline<d,T>::compute()
{
    getOptions();
}

template<short_t d,class T>
void gsApproxC1Spline<d,T>::matrix_into(gsSparseMatrix<T> & matrix) const
{
    m_matrix = gsSparseMatrix<T>();

}

template<short_t d,class T>
gsMultiBasis<T> gsApproxC1Spline<d,T>::localBasis() const
{
    gsMultiBasis<T> basis;
    return basis;
}




} // namespace gismo
