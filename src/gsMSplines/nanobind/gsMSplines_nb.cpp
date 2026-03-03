/** @file gsMSplines_nb.cpp

    @brief NanoBind bindings for the gsMSplines module

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_NANOBIND

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <gismo.h>

namespace nb = nanobind;
using namespace gismo;

NB_MODULE(msplines, m) {
    m.doc() = "G+Smo MSplines module";

    using MBasis2 = gsMappedBasis<2, real_t>;
    nb::class_<MBasis2, gsFunctionSet<real_t>>(m, "gsMappedBasis2")
        .def("size", static_cast<index_t (MBasis2::*)() const>(&MBasis2::size))
        .def("nPatches", &MBasis2::nPatches)
        .def("__repr__", [](const MBasis2& b) {
            std::ostringstream os; os << b; return os.str();
        });
}

#endif // GISMO_WITH_NANOBIND
