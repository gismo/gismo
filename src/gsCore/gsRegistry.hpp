/** @file gsRegistry.hpp

    @brief Implementation of gsRegistry (the per-Base singleton accessor).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsRegistry.h>

namespace gismo
{

template <class Base>
gsRegistry<Base> & gsRegistry<Base>::get()
{
    // Heap-allocated and never destroyed: registrations may fire from
    // static initializers of other libraries, so this object must
    // outlive every static-destruction order.
    static gsRegistry * r = new gsRegistry;
    return *r;
}

} // namespace gismo
