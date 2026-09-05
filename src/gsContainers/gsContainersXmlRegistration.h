/** @file gsContainersXmlRegistration.h

    @brief Header-only-mode hookup of the gsContainers XML specializations
    (gsMultiPatch, gsMultiBasis). These types are read by tag directly,
    not through an abstract-base dispatch, so nothing registers here.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#ifndef GISMO_BUILD_LIB
#include <gsContainers/gsContainersXml.hpp>
#endif
