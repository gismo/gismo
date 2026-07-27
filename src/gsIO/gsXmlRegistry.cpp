/** @file gsXmlRegistry.cpp

    @brief The per-process XML reader/writer registry instance.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsIO/gsXmlRegistry.h>

#include <algorithm>

namespace gismo
{
namespace internal
{

gsXmlRegistry & gsXmlRegistry::get()
{
    // Heap-allocated, never destroyed: registrations fire from static
    // initializers of other libraries (see gsRegistry for rationale)
    static gsXmlRegistry * r = new gsXmlRegistry;
    return *r;
}

} // namespace internal
} // namespace gismo
