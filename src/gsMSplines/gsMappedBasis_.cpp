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

    .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &                   ) const> (&Class::evalSingle     ), "Evaluates the basis function i")
    .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &                   ) const> (&Class::evalSingle     ), "Evaluates the basis function i")

    // Member functions
    // .def("eval_into"  , static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::eval_into)  , "Evaluates the mapped basis into a matrix")
    // .def("deriv_into" , static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::deriv_into) , "Evaluates the first derivatives into a matrix")
    // .def("deriv2_into", static_cast<void (Class::*)(const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::deriv2_into), "Evaluates the second derivatives into a matrix")

    // .def("evalSingle_into"  , static_cast<void (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::evalSingle_into)  , "Evaluates the mapped basis into a matrix")
    // .def("derivSingle_into" , static_cast<void (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::derivSingle_into) , "Evaluates the first derivatives into a matrix")
    // .def("deriv2Single_into", static_cast<void (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &, gsMatrix<real_t> &) const> (&Class::deriv2Single_into), "Evaluates the second derivatives into a matrix")

    .def("piece", &Class::piece, "Returns a piece")

    .def("eval", &Class::eval, "Evaluates the function set and returns a matrix")
    .def("deriv" , &Class::deriv , "Evaluates the first derivatives into a matrix")
    .def("deriv2", &Class::deriv2, "Evaluates the second derivatives into a matrix")

    // .def("evalSingle"  , static_cast<void (Base::*)(const index_t, const index_t, const gsMatrix<real_t> &) const> (&Base::evalSingle)  , "Evaluates the mapped basis into a matrix")
    // .def("derivSingle" , static_cast<void (Base::*)(const index_t, const index_t, const gsMatrix<real_t> &) const> (&Base::derivSingle) , "Evaluates the first derivatives into a matrix")
    // .def("deriv2Single", static_cast<void (Base::*)(const index_t, const index_t, const gsMatrix<real_t> &) const> (&Base::deriv2Single), "Evaluates the second derivatives into a matrix")

    ;
}

void pybind11_init_gsMappedSingleBasis2(py::module &m)
{
    using Base  = gsBasis<real_t>;
    using Class = gsMappedSingleBasis<2,real_t>;
    py::class_<Class,Base>(m, "gsMappedSingleBasis2")

    // Constructors
    .def(py::init<gsMappedBasis<2,real_t> *, const unsigned >() )
    ;
}


#endif

} // end namespace gismo
