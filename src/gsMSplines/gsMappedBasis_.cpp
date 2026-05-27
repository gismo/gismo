/** @file gsMappedBasis.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedBasis.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST gsMappedBasis<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

template <short_t d>
void pybind11_init_gsMappedBasis(py::module &m)
{
    using Base  = gsFunctionSet<real_t>;
    using Class = gsMappedBasis<d,real_t>;
    py::class_<Class,Base>(m, ("gsMappedBasis" + std::to_string(d)).c_str())

    .def(py::init<gsMultiBasis<real_t> const &, const gsSparseMatrix<real_t> >() )
    .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(const index_t, const index_t, const gsMatrix<real_t> &) const> (&Class::evalSingle), "Evaluates the basis function i")
    .def("piece", &Class::piece, "Returns a piece")
    .def("eval", &Class::eval, "Evaluates the function set and returns a matrix")
    .def("deriv", &Class::deriv, "Evaluates the first derivatives into a matrix")
    .def("deriv2", &Class::deriv2, "Evaluates the second derivatives into a matrix")
    ;
}

#define INST(D) template void pybind11_init_gsMappedBasis<D>(py::module &m);
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

template <short_t d>
void pybind11_init_gsMappedSingleBasis(py::module &m)
{
    using Base  = gsBasis<real_t>;
    using Class = gsMappedSingleBasis<d,real_t>;
    py::class_<Class,Base>(m, ("gsMappedSingleBasis" + std::to_string(d)).c_str())
    .def(py::init<gsMappedBasis<d,real_t> *, const unsigned >() )
    ;
}

#define INST(D) template void pybind11_init_gsMappedSingleBasis<D>(py::module &m);
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

#endif

} // end namespace gismo
