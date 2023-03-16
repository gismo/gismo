/** @file gsParametrization_.cpp

    @brief Instatiation of the gsParametrization class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsParametrization.h>
#include <gsModeling/gsParametrization.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsParametrization<real_t>;


#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;
void pybind11_init_gsParametrization(py::module &m)
{
  using Class  = gsParametrization<real_t>;
  py::class_<Class>(m, "gsParametrization")
    // Constructors
    // gsParametrization(const gsMesh<T> &mesh, const gsOptionList & list = defaultOptions());
    .def(py::init<>())
    // .def( py::init<const gsMesh<real_t> &, const gsOptionList &>())
    // Member functions
    // .def("compute", &Class::compute, "Main function which performs the computation.")
    // .def("createUVmatrix", &Class::createUVmatrix, "Return the parametric coordinates in [0,1]^2.")
    // .def("createXYZmatrix", &Class::createXYZmatrix, "Retunr the corresponding mapped values in R^33 to parametric coordinates in [0,1]^2.")
    ;
}
#endif



} // namespace gismo
