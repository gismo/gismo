/** @file gsAssembler_nb.cpp

    @brief NanoBind bindings for the gsAssembler module

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

NB_MODULE(assembler, m) {
    m.doc() = "G+Smo Assembler module";

    using BiHarm = gsExprAssembler<real_t>;
    nb::class_<BiHarm>(m, "gsExprAssembler")
        .def(nb::init<index_t, index_t>(), nb::arg("rows") = 1, nb::arg("cols") = 1)
        .def("numDofs", &BiHarm::numDofs)
        .def("matrix", &BiHarm::matrix, nb::rv_policy::reference_internal)
        .def("rhs", &BiHarm::rhs, nb::rv_policy::reference_internal)
        .def("options", static_cast<gsOptionList& (BiHarm::*)()>(&BiHarm::options),
             nb::rv_policy::reference_internal)
        .def("__repr__", [](const BiHarm& a) {
            std::ostringstream os;
            os << "gsExprAssembler(dofs=" << a.numDofs() << ")";
            return os.str();
        });
}

#endif // GISMO_WITH_NANOBIND
