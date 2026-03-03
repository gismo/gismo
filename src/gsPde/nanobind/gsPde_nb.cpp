/** @file gsPde_nb.cpp

    @brief NanoBind bindings for the gsPde module

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

NB_MODULE(gsPde, m) {
    m.doc() = "G+Smo PDE module";

    nb::enum_<condition_type::type>(m, "bctype")
        .value("dirichlet", condition_type::dirichlet)
        .value("neumann", condition_type::neumann)
        .value("robin", condition_type::robin)
        .value("clamped", condition_type::clamped)
        .value("collapsed", condition_type::collapsed)
        .export_values();

    using BC = gsBoundaryConditions<real_t>;
    nb::class_<BC>(m, "gsBoundaryConditions")
        .def(nb::init<>())
        .def("addCondition", static_cast<void (BC::*)(int, boxSide, condition_type::type, gsFunctionSet<real_t>*, short_t, bool, int)>(
                 &BC::addCondition),
             nb::arg("patch"), nb::arg("side"), nb::arg("type"),
             nb::arg("f"), nb::arg("unknown") = 0,
             nb::arg("parametric") = false, nb::arg("comp") = -1)
        .def("numConditions", [](const BC& bc) { return bc.size(); })
        .def("clear", &BC::clear)
        .def("__repr__", [](const BC& bc) {
            std::ostringstream os; os << bc; return os.str();
        });

    using PL = gsPointLoads<real_t>;
    nb::class_<PL>(m, "gsPointLoads")
        .def(nb::init<>())
        .def("addLoad", static_cast<void (PL::*)(const gsVector<real_t> &, const gsVector<real_t> &, int, bool)>(&PL::addLoad),
             nb::arg("point"), nb::arg("value"), nb::arg("patch") = 0,
             nb::arg("parametric") = true)
        .def("addLoad", static_cast<void (PL::*)(const gsVector<real_t> &, real_t, int, bool)>(&PL::addLoad),
             nb::arg("point"), nb::arg("value"), nb::arg("patch") = 0,
             nb::arg("parametric") = true)
        .def("numLoads", &PL::numLoads)
        .def("clear", &PL::clear);
}

#endif // GISMO_WITH_NANOBIND
