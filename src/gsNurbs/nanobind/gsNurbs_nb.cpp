/** @file gsNurbs_nb.cpp

    @brief NanoBind bindings for the gsNurbs module

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_NANOBIND

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gismo.h>

namespace nb = nanobind;
using namespace gismo;

NB_MODULE(nurbs, m) {
    m.doc() = "G+Smo Nurbs module";

    nb::class_<gsKnotVector<real_t>>(m, "gsKnotVector")
        .def(nb::init<>())
        .def(nb::init<real_t, real_t, unsigned, index_t, index_t, short_t>(),
             nb::arg("u0"), nb::arg("u1"), nb::arg("interior"),
             nb::arg("mult_ends") = 1, nb::arg("mult_interior") = 1,
             nb::arg("degree") = -1)
        .def("size", &gsKnotVector<real_t>::size)
        .def("uSize", &gsKnotVector<real_t>::uSize)
        .def("first", &gsKnotVector<real_t>::first)
        .def("last", &gsKnotVector<real_t>::last)
        .def("degree", &gsKnotVector<real_t>::degree)
        .def("numElements", &gsKnotVector<real_t>::numElements)
        .def("uniformRefine", &gsKnotVector<real_t>::uniformRefine,
             nb::arg("numKnots") = 1, nb::arg("mul") = 1)
        .def("degreeElevate", &gsKnotVector<real_t>::degreeElevate, nb::arg("i") = 1)
        .def("__repr__", [](const gsKnotVector<real_t>& kv) {
            std::ostringstream os; os << kv; return os.str();
        });

    using BSplineBasis = gsBSplineBasis<real_t>;
    nb::class_<BSplineBasis, gsBasis<real_t>>(m, "gsBSplineBasis")
        .def(nb::init<real_t, real_t, unsigned, int, unsigned>(),
             nb::arg("start") = 0.0, nb::arg("end") = 1.0,
             nb::arg("interior") = 0, nb::arg("degree") = 2,
             nb::arg("mult_interior") = 1)
        .def("knots", static_cast<const gsKnotVector<real_t>& (BSplineBasis::*)(int) const>(&BSplineBasis::knots),
             nb::arg("i") = 0, nb::rv_policy::reference_internal);

    using BSpline = gsBSpline<real_t>;
    nb::class_<BSpline, gsGeometry<real_t>>(m, "gsBSpline")
        .def(nb::init<>())
        .def("knots", static_cast<const gsKnotVector<real_t>& (BSpline::*)(int) const>(&BSpline::knots),
             nb::arg("i") = 0, nb::rv_policy::reference_internal)
        .def("insertKnot", static_cast<void (BSpline::*)(real_t, index_t)>(&BSpline::insertKnot),
             nb::arg("knot"), nb::arg("i") = 1);

    using TBBasis2 = gsTensorBSplineBasis<2, real_t>;
    using TBBasis3 = gsTensorBSplineBasis<3, real_t>;
    nb::class_<TBBasis2, gsBasis<real_t>>(m, "gsTensorBSplineBasis2");
    nb::class_<TBBasis3, gsBasis<real_t>>(m, "gsTensorBSplineBasis3");

    using TBSpline2 = gsTensorBSpline<2, real_t>;
    using TBSpline3 = gsTensorBSpline<3, real_t>;
    using TBSpline4 = gsTensorBSpline<4, real_t>;
    nb::class_<TBSpline2, gsGeometry<real_t>>(m, "gsTensorBSpline2")
        .def(nb::init<>());
    nb::class_<TBSpline3, gsGeometry<real_t>>(m, "gsTensorBSpline3")
        .def(nb::init<>());
    nb::class_<TBSpline4, gsGeometry<real_t>>(m, "gsTensorBSpline4")
        .def(nb::init<>());

    nb::class_<gsNurbsCreator<real_t>>(m, "gsNurbsCreator")
        .def_static("BSplineSquare", [](real_t r, real_t x, real_t y) {
                return *gsNurbsCreator<real_t>::BSplineSquare(r, x, y);
            }, nb::arg("r") = 1, nb::arg("x") = 0, nb::arg("y") = 0)
        .def_static("BSplineCube", [](real_t r, real_t x, real_t y, real_t z) {
                return *gsNurbsCreator<real_t>::BSplineCube(r, x, y, z);
            }, nb::arg("r") = 1, nb::arg("x") = 0, nb::arg("y") = 0, nb::arg("z") = 0);
}

#endif // GISMO_WITH_NANOBIND
