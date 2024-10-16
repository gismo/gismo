/** @file gsSpringPatch_.cpp

    @brief Provides spring patch construction from boundary data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include<gsModeling/gsSpringPatch.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsSpringPatch<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsSpringPatch(py::module &m)
{
    using Class = gsSpringPatch<real_t>;
    py::class_<Class>(m, "gsSpringPatch")

    // Constructors
    .def(py::init<const gsMultiPatch<real_t> &>()) //default arguments

    // Member functions
    .def("compute", &Class::compute, "Computes the spring patch.", py::return_value_policy::reference_internal)
    .def("result", &Class::result, "Get the resulting spring patch.", py::return_value_policy::reference_internal)
    ;
}
#endif

}// namespace gismo
