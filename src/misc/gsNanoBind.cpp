/** @file gsNanoBind.cpp

    @brief NanoBind main module file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_NANOBIND

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.attr("__version__") = GISMO_VERSION;
  m.doc() = "G+Smo (Geometry + Simulation Modules)";
}

#endif // GISMO_WITH_NANOBIND
