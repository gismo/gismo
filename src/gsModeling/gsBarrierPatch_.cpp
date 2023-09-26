/** @file gsBarrierPatch_.cpp
 *
    @brief - A reference implementation of the following paper.
	If you make use of the code or the idea/algorithm in your work, please cite our paper
	Ji, Y., Yu, Y. Y., Wang, M. Y., & Zhu, C. G. (2021).
	Constructing high-quality planar NURBS parameterization for
	isogeometric analysis by adjustment control points and weights.
	Journal of Computational and Applied Mathematics, 396, 113615.
	(https://www.sciencedirect.com/science/article/pii/S0377042721002375)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji (jiye@mail.dlut.edu.cn), H.M. Verhelst
*/

#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsBarrierCore.h>
#include <gsModeling/gsBarrierCore.hpp>

#include <gsModeling/gsBarrierPatch.h>
#include <gsModeling/gsBarrierPatch.hpp>

//#include <gsModeling/gsBarrierPatchGenerator.h>
//#include <gsModeling/gsBarrierPatchGenerator.hpp>

namespace gismo
{
	CLASS_TEMPLATE_INST gsBarrierCore<2,real_t>;
	CLASS_TEMPLATE_INST gsBarrierCore<3,real_t>;

	CLASS_TEMPLATE_INST gsBarrierPatch<2,real_t>;
	CLASS_TEMPLATE_INST gsBarrierPatch<3,real_t>;

//	CLASS_TEMPLATE_INST gsBarrierPatchGenerator<2,real_t>;
//	CLASS_TEMPLATE_INST gsBarrierPatchGenerator<3,real_t>;

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

	void pybind11_init_gsBarrierPatch2(py::module &m)
	{
		using Class = gsBarrierPatch<2,real_t>;
		py::class_<Class>(m, "gsBarrierPatch2")

		// Constructors
		.def(py::init<const gsMultiPatch<real_t> &, const gsDofMapper &>())
		.def(py::init<const gsMultiPatch<real_t> &>(),
		      py::arg("patchWise") = true) //default arguments

		// Member functions
		.def("setMapper", &Class::setMapper, "Sets the mapper.")
		.def("compute", &Class::compute, "Computes analysis-suitable parameterizations using different methods.")
		.def("result", &Class::result, "Returns the result in a multi-patch format.")
		.def("options", &Class::options, "Returns the options list.")
		.def("defaultOptions", &Class::defaultOptions, "Sets the default options.")
		;
  }

	void pybind11_init_gsBarrierPatch3(py::module &m)
	{
		using Class = gsBarrierPatch<3,real_t>;
		py::class_<Class>(m, "gsBarrierPatch3")

	   // Constructors
		.def(py::init<const gsMultiPatch<real_t> &, const gsDofMapper &>())
		.def(py::init<const gsMultiPatch<real_t> &>(),
	         py::arg("patchWise") = true) //default arguments

	   // Member functions
	   .def("setMapper", &Class::setMapper, "Sets the mapper.")
	   .def("compute", &Class::compute, "Computes analysis-suitable parameterizations using different methods.")
	   .def("result", &Class::result, "Returns the result in a multi-patch format.")
	   .def("options", &Class::options, "Returns the options list.")
	   .def("defaultOptions", &Class::defaultOptions, "Sets the default options.")
	   ;
  }

#endif
} // namespace gismo
