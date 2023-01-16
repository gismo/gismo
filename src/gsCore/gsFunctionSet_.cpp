/** @file gsFunctionSet_.cpp

    @brief instantiation of gsFunctionSet

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFunctionSet.hpp>

namespace gismo {

CLASS_TEMPLATE_INST gsFunctionSet<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsFunctionSet(py::module &m)
{
  using Class = gsFunctionSet<real_t>;
  py::class_<Class>(m, "gsFunctionSet")

  //Constructor?
  .def("eval_into", &Class::eval_into, "Evaluates the function set into a matrix")
  .def("deriv_into", &Class::deriv_into, "Evaluates the first derivative into a matrix")
  .def("deriv2_into", &Class::deriv2_into, "Evaluates the second derivative into a matrix")
  .def("evalAllDers_into", &Class::deriv2_into, "Evaluates all derivatives upto certien order into a vector of matrices")
  .def("eval", &Class::eval, "Evaluates the function set and returns a matrix")
  .def("deriv", &Class::deriv, "Evaluates the first derivative and returns a matrix")
  .def("deriv2", &Class::deriv2, "Evaluates the second derivative and returns a matrix")
  .def("evalAllDers", &Class::evalAllDers, "Evaluates all derivatives upto certien order into a vector of matrices")
  .def("domainDim", &Class::domainDim, "Returns the domain dimension")
  .def("targetDim", &Class::targetDim, "Returns the target dimension")

  ;
}
#endif


}

namespace std {

TEMPLATE_INST void swap(gismo::gsFuncData<real_t> & f1, gismo::gsFuncData<real_t> & f2);


}
