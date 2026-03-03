/** @file gsCore_nb.cpp

    @brief NanoBind bindings for the gsCore module

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
#include <nanobind/stl/pair.h>

#include <gismo.h>

namespace nb = nanobind;

NB_MODULE(gsCore, m) {
    m.doc() = "G+Smo Core module";

    using namespace gismo;

    nb::enum_<boundary::side>(m, "side")
        .value("west",  boundary::west)
        .value("east",  boundary::east)
        .value("south", boundary::south)
        .value("north", boundary::north)
        .value("front", boundary::front)
        .value("back",  boundary::back)
        .value("stime", boundary::stime)
        .value("etime", boundary::etime)
        .value("none",  boundary::none)
        .export_values();

    nb::class_<boxSide>(m, "boxSide")
        .def(nb::init<>())
        .def(nb::init<short_t>())
        .def("direction", &boxSide::direction)
        .def("parameter", &boxSide::parameter)
        .def("index", static_cast<short_t (boxSide::*)() const>(&boxSide::index))
        .def("opposite", &boxSide::opposite)
        .def("__repr__", [](const boxSide& s) {
            std::ostringstream os; os << "boxSide(" << s.index() << ")"; return os.str();
        });

    nb::class_<gsFunctionSet<real_t>>(m, "gsFunctionSet")
        .def("domainDim", &gsFunctionSet<real_t>::domainDim)
        .def("targetDim", &gsFunctionSet<real_t>::targetDim)
        .def("size", &gsFunctionSet<real_t>::size)
        .def("nPieces", &gsFunctionSet<real_t>::nPieces)
        .def("support", &gsFunctionSet<real_t>::support)
        .def("eval", &gsFunctionSet<real_t>::eval)
        .def("deriv", &gsFunctionSet<real_t>::deriv)
        .def("deriv2", &gsFunctionSet<real_t>::deriv2)
        .def("active", &gsFunctionSet<real_t>::active)
        .def("__repr__", [](const gsFunctionSet<real_t>& f) {
            std::ostringstream os; f.print(os); return os.str();
        });

    nb::class_<gsFunction<real_t>, gsFunctionSet<real_t>>(m, "gsFunction")
        .def("jacobian", &gsFunction<real_t>::jacobian)
        .def("hessian", &gsFunction<real_t>::hessian, nb::arg("u"), nb::arg("coord") = 0)
        .def("laplacian", &gsFunction<real_t>::laplacian);

    nb::class_<gsFunctionExpr<real_t>, gsFunction<real_t>>(m, "gsFunctionExpr")
        .def(nb::init<const std::string &, short_t>(),
             nb::arg("expression"), nb::arg("ddim"))
        .def(nb::init<const std::string &, const std::string &, short_t>(),
             nb::arg("expr1"), nb::arg("expr2"), nb::arg("ddim"))
        .def(nb::init<const std::string &, const std::string &, const std::string &, short_t>(),
             nb::arg("expr1"), nb::arg("expr2"), nb::arg("expr3"), nb::arg("ddim"))
        .def("addComponent", &gsFunctionExpr<real_t>::addComponent);

    nb::class_<gsBasis<real_t>, gsFunctionSet<real_t>>(m, "gsBasis")
        .def("dim", &gsBasis<real_t>::dim)
        .def("numElements", [](const gsBasis<real_t>& b) { return b.numElements(); })
        .def("degree", &gsBasis<real_t>::degree, nb::arg("i") = 0)
        .def("anchors", static_cast<gsMatrix<real_t> (gsBasis<real_t>::*)() const>(&gsBasis<real_t>::anchors))
        .def("evalSingle", &gsBasis<real_t>::evalSingle)
        .def("derivSingle", &gsBasis<real_t>::derivSingle)
        .def("uniformRefine", &gsBasis<real_t>::uniformRefine,
             nb::arg("numKnots") = 1, nb::arg("mul") = 1, nb::arg("dir") = -1)
        .def("degreeElevate", &gsBasis<real_t>::degreeElevate,
             nb::arg("elevationSteps") = 1, nb::arg("dir") = -1)
        .def("degreeReduce", &gsBasis<real_t>::degreeReduce,
             nb::arg("reductionSteps") = 1, nb::arg("dir") = -1);

    nb::class_<gsGeometry<real_t>, gsFunction<real_t>>(m, "gsGeometry")
        .def("basis", static_cast<const gsBasis<real_t>& (gsGeometry<real_t>::*)() const>(&gsGeometry<real_t>::basis),
             nb::rv_policy::reference_internal)
        .def("coefs", static_cast<const gsMatrix<real_t>& (gsGeometry<real_t>::*)() const>(&gsGeometry<real_t>::coefs),
             nb::rv_policy::reference_internal)
        .def("setCoefs", &gsGeometry<real_t>::setCoefs)
        .def("coefsSize", &gsGeometry<real_t>::coefsSize)
        .def("coefDim", &gsGeometry<real_t>::coefDim)
        .def("geoDim", &gsGeometry<real_t>::geoDim)
        .def("parDim", &gsGeometry<real_t>::parDim)
        .def("coDim", &gsGeometry<real_t>::coDim)
        .def("parameterRange", &gsGeometry<real_t>::parameterRange)
        .def("uniformRefine", &gsGeometry<real_t>::uniformRefine,
             nb::arg("numKnots") = 1, nb::arg("mul") = 1, nb::arg("dir") = -1)
        .def("degreeElevate", &gsGeometry<real_t>::degreeElevate,
             nb::arg("elevationSteps") = 1, nb::arg("dir") = -1)
        .def("translate", &gsGeometry<real_t>::translate)
        .def("scale", static_cast<void (gsGeometry<real_t>::*)(real_t, int)>(&gsGeometry<real_t>::scale),
             nb::arg("s"), nb::arg("coord") = -1)
        .def("toggleOrientation", &gsGeometry<real_t>::toggleOrientation);

    nb::class_<gsBoxTopology>(m, "gsBoxTopology")
        .def(nb::init<>())
        .def("nBoxes", &gsBoxTopology::nBoxes)
        .def("dim", &gsBoxTopology::dim)
        .def("nInterfaces", &gsBoxTopology::nInterfaces)
        .def("nBoundary", &gsBoxTopology::nBoundary)
        .def("__repr__", [](const gsBoxTopology& t) {
            std::ostringstream os;
            os << "gsBoxTopology(boxes=" << t.nBoxes() << ", interfaces=" << t.nInterfaces()
               << ", boundary=" << t.nBoundary() << ")";
            return os.str();
        });

    using MultiPatch = gsMultiPatch<real_t>;
    nb::class_<MultiPatch, gsBoxTopology>(m, "gsMultiPatch")
        .def(nb::init<>())
        .def("domainDim", &MultiPatch::domainDim)
        .def("targetDim", &MultiPatch::targetDim)
        .def("nPatches", &MultiPatch::nPatches)
        .def("patch", static_cast<gsGeometry<real_t>& (MultiPatch::*)(size_t) const>(&MultiPatch::patch),
             nb::rv_policy::reference_internal)
        .def("addPatch", static_cast<index_t (MultiPatch::*)(const gsGeometry<real_t>&)>(&MultiPatch::addPatch))
        .def("degreeElevate", &MultiPatch::degreeElevate,
             nb::arg("elevationSteps") = 1, nb::arg("dir") = -1)
        .def("uniformRefine", &MultiPatch::uniformRefine,
             nb::arg("numKnots") = 1, nb::arg("mul") = 1, nb::arg("dir") = -1)
        .def("basis", static_cast<gsBasis<real_t>& (MultiPatch::*)(const size_t)>(&MultiPatch::basis),
             nb::rv_policy::reference_internal)
        .def("computeTopology", static_cast<bool (MultiPatch::*)(real_t, bool, bool)>(&MultiPatch::computeTopology))
        .def("fixOrientation", &MultiPatch::fixOrientation)
        .def("clear", &MultiPatch::clear);

    using MultiBasis = gsMultiBasis<real_t>;
    nb::class_<MultiBasis, gsFunctionSet<real_t>>(m, "gsMultiBasis")
        .def(nb::init<>())
        .def(nb::init<const MultiPatch&, bool>(),
             nb::arg("mp"), nb::arg("numeratorOnly") = false)
        .def("nBases", &MultiBasis::nBases)
        .def("dim", &MultiBasis::dim)
        .def("degree", &MultiBasis::degree, nb::arg("i") = 0, nb::arg("comp") = 0)
        .def("maxDegree", &MultiBasis::maxDegree)
        .def("minDegree", &MultiBasis::minDegree)
        .def("basis", static_cast<const gsBasis<real_t>& (MultiBasis::*)(const size_t) const>(&MultiBasis::basis),
             nb::rv_policy::reference_internal)
        .def("uniformRefine", static_cast<void (MultiBasis::*)(int, int, int)>(&MultiBasis::uniformRefine),
             nb::arg("numKnots") = 1, nb::arg("mul") = 1, nb::arg("dir") = -1)
        .def("degreeElevate", static_cast<void (MultiBasis::*)(short_t, short_t)>(&MultiBasis::degreeElevate),
             nb::arg("elevationSteps") = 1, nb::arg("dir") = -1)
        .def("degreeReduce", static_cast<void (MultiBasis::*)(short_t)>(&MultiBasis::degreeReduce),
             nb::arg("reductionSteps") = 1);

    nb::class_<gsDofMapper>(m, "gsDofMapper")
        .def(nb::init<>())
        .def("finalize", &gsDofMapper::finalize)
        .def("isFinalized", &gsDofMapper::isFinalized)
        .def("freeSize", static_cast<index_t (gsDofMapper::*)() const>(&gsDofMapper::freeSize))
        .def("boundarySize", static_cast<index_t (gsDofMapper::*)() const>(&gsDofMapper::boundarySize))
        .def("coupledSize", &gsDofMapper::coupledSize)
        .def("taggedSize", &gsDofMapper::taggedSize)
        .def("size", static_cast<index_t (gsDofMapper::*)() const>(&gsDofMapper::size))
        .def("mapSize", &gsDofMapper::mapSize)
        .def("index", static_cast<index_t (gsDofMapper::*)(index_t, index_t, index_t) const>(&gsDofMapper::index),
             nb::arg("i"), nb::arg("k") = 0, nb::arg("c") = 0)
        .def("is_free", static_cast<bool (gsDofMapper::*)(index_t, index_t, index_t) const>(&gsDofMapper::is_free),
             nb::arg("i"), nb::arg("k") = 0, nb::arg("c") = 0)
        .def("__repr__", [](const gsDofMapper& d) {
            std::ostringstream os;
            if (d.isFinalized())
                os << "gsDofMapper(size=" << d.size() << ", free=" << d.freeSize()
                   << ", boundary=" << d.boundarySize() << ")";
            else
                os << "gsDofMapper(not finalized)";
            return os.str();
        });
}

#endif // GISMO_WITH_NANOBIND
