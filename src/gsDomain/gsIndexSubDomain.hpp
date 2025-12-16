/** @file gsIndexSubDomain.hpp

    @brief Implementation of gsIndexSubDomain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsIndexSubDomain.h>

namespace gismo
{

template<class T>
typename gsIndexSubDomain<T>::domainIter
gsIndexSubDomain<T>::beginAll() const
{
    return domainIter(new gsIndexSubDomainIterator<T>(*this));
}

} // namespace gismo
