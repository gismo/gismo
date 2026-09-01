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

// Explicit instantiation of template member functions for gsHBox (d=2,3 only)
template void gsParaview<real_t>::write<2>(const gsHBox<2,real_t> &, const std::string &) const;
template void gsParaview<real_t>::write<3>(const gsHBox<3,real_t> &, const std::string &) const;

// Explicit instantiation of template member functions for gsHBoxContainer (d=2,3 only)
template void gsParaview<real_t>::write<2>(const gsHBoxContainer<2,real_t> &, const std::string &) const;
template void gsParaview<real_t>::write<3>(const gsHBoxContainer<3,real_t> &, const std::string &) const;

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
    .def("write", static_cast<void (Class::*)(const gsGeometry<real_t> &, const std::string &) const>
         (&Class::write),
         py::arg("geo"), py::arg("fn") = "geometry",
         "Export a gsGeometry to Paraview file")

    .def("write", static_cast<void (Class::*)(const gsMultiPatch<real_t> &, const std::string &) const>
         (&Class::write),
         py::arg("mp"), py::arg("fn") = "multipatch",
         "Export a gsMultiPatch to Paraview file")

    .def("write", static_cast<void (Class::*)(const gsField<real_t> &, const std::string &) const>
         (&Class::write),
         py::arg("field"), py::arg("fn") = "field",
         "Export a gsField to Paraview file")

    .def("write", static_cast<void (Class::*)(const gsBasis<real_t> &, const std::string &) const>
         (&Class::write),
         py::arg("basis"), py::arg("fn") = "basis",
         "Export a gsBasis to Paraview file")

    .def("write", static_cast<void (Class::*)(const gsMesh<real_t> &, const std::string &) const>
         (&Class::write),
         py::arg("mesh"), py::arg("fn") = "mesh",
         "Export a gsMesh to Paraview file")

    // Points
    .def("writePoints", static_cast<void (Class::*)(const gsMatrix<real_t> &, const std::string &) const>
         (&Class::writePoints),
         py::arg("points"), py::arg("fn") = "points",
         "Export points to Paraview file")

    ;
}

#endif // GISMO_WITH_PYBIND11

} // namespace gismo
