/** @file gsHSplines_nb.cpp

    @brief NanoBind bindings for the gsHSplines module

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

NB_MODULE(gsHSplines, m) {
    m.doc() = "G+Smo Hierarchical Splines module";

    nb::class_<gsHTensorBasis<2, real_t>, gsBasis<real_t>>(m, "gsHTensorBasis2")
        .def("maxLevel", &gsHTensorBasis<2, real_t>::maxLevel)
        .def("treeSize", &gsHTensorBasis<2, real_t>::treeSize);
    nb::class_<gsHTensorBasis<3, real_t>, gsBasis<real_t>>(m, "gsHTensorBasis3")
        .def("maxLevel", &gsHTensorBasis<3, real_t>::maxLevel)
        .def("treeSize", &gsHTensorBasis<3, real_t>::treeSize);

    nb::class_<gsTHBSplineBasis<2, real_t>, gsHTensorBasis<2, real_t>>(m, "gsTHBSplineBasis2");
    nb::class_<gsTHBSplineBasis<3, real_t>, gsHTensorBasis<3, real_t>>(m, "gsTHBSplineBasis3");

    nb::class_<gsHBSplineBasis<2, real_t>, gsHTensorBasis<2, real_t>>(m, "gsHBSplineBasis2");
    nb::class_<gsHBSplineBasis<3, real_t>, gsHTensorBasis<3, real_t>>(m, "gsHBSplineBasis3");

    nb::class_<gsTHBSpline<2, real_t>, gsGeometry<real_t>>(m, "gsTHBSpline2")
        .def(nb::init<>());
    nb::class_<gsTHBSpline<3, real_t>, gsGeometry<real_t>>(m, "gsTHBSpline3")
        .def(nb::init<>());

    nb::class_<gsHBSpline<2, real_t>, gsGeometry<real_t>>(m, "gsHBSpline2")
        .def(nb::init<>());
    nb::class_<gsHBSpline<3, real_t>, gsGeometry<real_t>>(m, "gsHBSpline3")
        .def(nb::init<>());
}

#endif // GISMO_WITH_NANOBIND
