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
#include <gsCore/gsDimMacro.h>

#include <gsModeling/gsBarrierCore.h>
#include <gsModeling/gsBarrierCore.hpp>

#include <gsModeling/gsBarrierPatch.h>
#include <gsModeling/gsBarrierPatch.hpp>

//#include <gsModeling/gsBarrierPatchGenerator.h>
//#include <gsModeling/gsBarrierPatchGenerator.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsBarrierCore<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsBarrierPatch<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

//	CLASS_TEMPLATE_INST gsBarrierPatchGenerator<2,real_t>;
//	CLASS_TEMPLATE_INST gsBarrierPatchGenerator<3,real_t>;

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

	template <short_t d>
	void pybind11_init_gsBarrierPatch(py::module &m)
	{
		using Class = gsBarrierPatch<d,real_t>;
		py::class_<Class>(m, ("gsBarrierPatch" + std::to_string(d)).c_str())

		.def(py::init<const gsMultiPatch<real_t> &, const gsDofMapper &>())
		.def(py::init<const gsMultiPatch<real_t> &>())
		.def(py::init<const gsMultiPatch<real_t> &, bool>())

		.def("setMapper", &Class::setMapper, "Sets the mapper.")
		.def("compute", &Class::compute, "Computes analysis-suitable parameterizations using different methods.")
		.def("result", &Class::result, "Returns the result in a multi-patch format.")
		.def("options", &Class::options, "Returns the options list.")
		.def("defaultOptions", &Class::defaultOptions, "Sets the default options.")
		;
	}

#define INST(D) template void pybind11_init_gsBarrierPatch<D>(py::module &m);
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

#endif
} // namespace gismo
