/** @file gsParaview_.cpp

    @brief Instantiations for gsParaview class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsParaview.h>
#include <gsIO/gsParaview.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsParaview<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsParaview(py::module &m)
{
    using Class = gsParaview<real_t>;
    py::class_<Class>(m, "gsParaview")

    // Constructors
    .def(py::init<>())
    .def(py::init<const gsOptionList &>())

    // Options
    .def("options", static_cast<gsOptionList & (Class::*)()>(&Class::options),
         "Returns a reference to the options list")
    .def_static("defaultOptions", &Class::defaultOptions,
                "Returns the default options")

    // Write methods - basic types
    .def("write", py::overload_cast<const gsGeometry<real_t> &, const std::string &>
         (&Class::write, py::const_),
         py::arg("geo"), py::arg("fn") = "geometry",
         "Export a gsGeometry to Paraview file")

    .def("write", py::overload_cast<const gsMultiPatch<real_t> &, const std::string &>
         (&Class::write, py::const_),
         py::arg("mp"), py::arg("fn") = "multipatch",
         "Export a gsMultiPatch to Paraview file")

    .def("write", py::overload_cast<const gsField<real_t> &, const std::string &>
         (&Class::write, py::const_),
         py::arg("field"), py::arg("fn") = "field",
         "Export a gsField to Paraview file")

    .def("write", py::overload_cast<const gsBasis<real_t> &, const std::string &>
         (&Class::write, py::const_),
         py::arg("basis"), py::arg("fn") = "basis",
         "Export a gsBasis to Paraview file")

    .def("write", py::overload_cast<const gsMesh<real_t> &, const std::string &>
         (&Class::write, py::const_),
         py::arg("mesh"), py::arg("fn") = "mesh",
         "Export a gsMesh to Paraview file")

    // Points
    .def("writePoints", py::overload_cast<const gsMatrix<real_t> &, const std::string &>
         (&Class::writePoints, py::const_),
         py::arg("points"), py::arg("fn") = "points",
         "Export points to Paraview file")

    ;
}

#endif // GISMO_WITH_PYBIND11

} // namespace gismo
