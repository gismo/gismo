/** @file gsC1Basis.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsC1Basis.h>
#include <gsIO/gsFileData.h>

namespace gismo
{

template<short_t d,class T>
gsC1Basis<d,T>::defaultOptions()
{
    m_options.addInt("test","a test option",0);
}

template<short_t d,class T>
gsC1Basis<d,T>::getOptions()
{
    index_t test = m_options.getInt("test");
}

template<short_t d,class T>
gsC1Basis<d,T>::gsC1Basis(  gsMultiPatch<T> const & mp,
                                    index_t patchID )
:
m_patches(mp),
m_patchID(patchID)
{
    defaultOptions();
    getOptions();
}

// template<short_t d,class T>
// gsC1Basis<d,T>::gsC1Basis( gsC1Basis<d,T> other )
// :
// m_patches(other.m_patches),
// m_patchID(other.m_patchID),
// m_options(other.m_options)
// {

// }




} // namespace gismo
