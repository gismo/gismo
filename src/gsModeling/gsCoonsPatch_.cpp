/** @file gsCoonsPatch_.cpp

    @brief Provides Coons's patch construction from boundary data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner, A. Mantzaflaris
*/

#include <gsCore/gsTemplateTools.h>
#include <gsModeling/gsCoonsPatch.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsCoonsPatch<real_t>;


#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsCoonsPatch(py::module &m)
{
    using Class = gsCoonsPatch<real_t>;
    py::class_<Class>(m, "gsCoonsPatch")

    // Constructors
    .def(py::init<const gsMultiPatch<real_t> &>()) //default arguments

    // Member functions
    .def("compute", &Class::compute, "Computes the Coons patch.")
    .def("result", &Class::result, "Gets the resulting Coons patch.")
    ;
}
#endif

} // namespace gismo
