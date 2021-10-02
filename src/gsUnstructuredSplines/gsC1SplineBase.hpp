/** @file gsC1SplineBase.hpp

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
gsC1SplineBase<d,T>::defaultOptions()
{
    m_options.addInt("test","a test option",0);
}

template<short_t d,class T>
gsC1SplineBase<d,T>::getOptions()
{
    index_t test = m_options.getInt("test");
}


} // namespace gismo
