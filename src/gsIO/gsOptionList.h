/** @file gsOptionList.h

    @brief Forwarding header. gsOptionList moved to gsCore/gsOptionList.h
    (it is the universal option container, not an I/O concern). Its XML
    binding remains in gsIO/gsOptionListXml.h; both are included here to
    preserve the previous semantics of this header.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsOptionList.h>
#include <gsIO/gsOptionListXml.h>
