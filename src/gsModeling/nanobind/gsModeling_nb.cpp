/** @file gsModeling_nb.cpp

    @brief NanoBind bindings for the gsModeling module

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

NB_MODULE(gsModeling, m) {
    m.doc() = "G+Smo Modeling module";

    using Fitting = gsFitting<real_t>;
    nb::class_<Fitting>(m, "gsFitting")
        .def(nb::init<gsMatrix<real_t>, gsMatrix<real_t>, gsBasis<real_t>&>(),
             nb::arg("params"), nb::arg("points"), nb::arg("basis"))
        .def("compute", static_cast<void (Fitting::*)(real_t)>(&Fitting::compute),
             nb::arg("lambda") = 0)
        .def("computeErrors", static_cast<void (Fitting::*)()>(&Fitting::computeErrors))
        .def("minPointError", &Fitting::minPointError)
        .def("maxPointError", &Fitting::maxPointError)
        .def("result", &Fitting::result, nb::rv_policy::reference_internal)
        .def("__repr__", [](const Fitting& f) {
            return std::string("gsFitting(...)");
        });

    nb::class_<gsCoonsPatch<real_t>>(m, "gsCoonsPatch")
        .def(nb::init<const gsMultiPatch<real_t>&>())
        .def("compute", &gsCoonsPatch<real_t>::compute)
        .def("result", &gsCoonsPatch<real_t>::result, nb::rv_policy::reference_internal);

    nb::class_<gsSpringPatch<real_t>>(m, "gsSpringPatch")
        .def(nb::init<const gsMultiPatch<real_t>&>())
        .def("compute", &gsSpringPatch<real_t>::compute)
        .def("result", &gsSpringPatch<real_t>::result, nb::rv_policy::reference_internal);
}

#endif // GISMO_WITH_NANOBIND
