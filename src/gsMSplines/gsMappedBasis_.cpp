/** @file gsMappedBasis.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedBasis.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsMappedBasis<1,real_t> ;
CLASS_TEMPLATE_INST gsMappedBasis<2,real_t> ;
CLASS_TEMPLATE_INST gsMappedBasis<3,real_t> ;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsMappedBasis2(py::module &m)
{
    using Base  = gsFunctionSet<real_t>;
    using Class = gsMappedBasis<2,real_t>;
    py::class_<Class,Base>(m, "gsMappedBasis2")

    // Constructors
    .def(py::init<gsMultiBasis<real_t> const &, const gsSparseMatrix<real_t> >() )

    // Member functions
    .def("eval_into"  , static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::eval_into)  , "Evaluates the mapped basis into a matrix")
    .def("deriv_into" , static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::deriv_into) , "Evaluates the first derivatives into a matrix")
    .def("deriv2_into", static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::deriv2_into), "Evaluates the second derivatives into a matrix")
    ;
}

#endif

} // end namespace gismo
