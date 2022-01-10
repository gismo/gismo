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

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsFunctionSet(py::module &m)
{
  using Class = gsFunctionSet<real_t>;
  py::class_<Class>(m, "gsFunctionSet")

  //Constructor?
  ;
}
#endif


}

namespace std {

TEMPLATE_INST void swap(gismo::gsFuncData<real_t> & f1, gismo::gsFuncData<real_t> & f2);


}
