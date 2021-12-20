#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXml.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>

#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsHBSplineBasis.hpp>

#include <gsHSplines/gsHBSpline.h>
#include <gsHSplines/gsHBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsHTensorBasis <1,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <2,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <3,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <4,real_t>;

CLASS_TEMPLATE_INST gsHBSplineBasis <1,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <2,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <3,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <4,real_t>;

namespace internal
{

CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsHBSpline<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<4,real_t> >;
}

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsHBSplineBasis2(py::module &m)
{
  using Class = gsHBSplineBasis<2,real_t>;
  py::class_<Class>(m, "gsHBSplineBasis2")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())

    // Member functions
    .def("domainDim", 			&Class::domainDim,
    		"Returns the domain dimension"												)
    .def("eval_into", 			&Class::eval_into,
    		"Evaluates the values into a matrix"										)
    .def("deriv_into", 			&Class::deriv_into,
    		"Evaluates the derivatives into a matrix"									)
    .def("deriv2_into", 		&Class::deriv2_into,
    		"Evaluates the second derivatives into a matrix"							)
    .def("evalSingle_into", 	&Class::evalSingle_into,
    		"Evaluates the values of a single basis function into a matrix"				)
    .def("derivSingle_into", 	&Class::derivSingle_into,
    		"Evaluates the derivatives of a single basis functioninto a matrix" 		)
    .def("deriv2Single_into", 	&Class::deriv2Single_into,
    		"Evaluates the second derivatives of a single basis function into a matrix" )

    // Derived from gsHTensorBasis
    .def("size",&Class::size,"Returns the domain dimension")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(               ) const> (&Class::support), "Get the support of the basis")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(const index_t &) const> (&Class::support), "Get the support of the basis function with an index i")
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
    	py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::refine  ), "Refines the basis given a box")
    .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::unrefine), "Refines the basis given a box")
    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")
    // .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
    .def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
    .def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
    .def("asElements        ",&Class::asElements        ,"Returns the elements given refinement boxes")
    .def("asElementsUnrefine",&Class::asElementsUnrefine,"Returns the elements given refinement boxes")
    ;
}

void pybind11_init_gsHBSplineBasis3(py::module &m)
{
  using Class = gsHBSplineBasis<3,real_t>;
  py::class_<Class>(m, "gsHBSplineBasis3")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())

    // Member functions
    .def("domainDim", 			&Class::domainDim,
    		"Returns the domain dimension"												)
    .def("eval_into", 			&Class::eval_into,
    		"Evaluates the values into a matrix"										)
    .def("deriv_into", 			&Class::deriv_into,
    		"Evaluates the derivatives into a matrix"									)
    .def("deriv2_into", 		&Class::deriv2_into,
    		"Evaluates the second derivatives into a matrix"							)
    .def("evalSingle_into", 	&Class::evalSingle_into,
    		"Evaluates the values of a single basis function into a matrix"				)
    .def("derivSingle_into", 	&Class::derivSingle_into,
    		"Evaluates the derivatives of a single basis functioninto a matrix" 		)
    .def("deriv2Single_into", 	&Class::deriv2Single_into,
    		"Evaluates the second derivatives of a single basis function into a matrix" )

    // Derived from gsHTensorBasis
    .def("size",&Class::size,"Returns the domain dimension")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(               ) const> (&Class::support), "Get the support of the basis")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(const index_t &) const> (&Class::support), "Get the support of the basis function with an index i")
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
    	py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::refine  ), "Refines the basis given a box")
    .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::unrefine), "Refines the basis given a box")
    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")
    // .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
    .def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
    .def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
    .def("asElements        ",&Class::asElements        ,"Returns the elements given refinement boxes")
    .def("asElementsUnrefine",&Class::asElementsUnrefine,"Returns the elements given refinement boxes")
    ;
}

void pybind11_init_gsHBSplineBasis4(py::module &m)
{
  using Class = gsHBSplineBasis<4,real_t>;
  py::class_<Class>(m, "gsHBSplineBasis4")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())

    // Member functions
    .def("domainDim", 			&Class::domainDim,
    		"Returns the domain dimension"												)
    .def("eval_into", 			&Class::eval_into,
    		"Evaluates the values into a matrix"										)
    .def("deriv_into", 			&Class::deriv_into,
    		"Evaluates the derivatives into a matrix"									)
    .def("deriv2_into", 		&Class::deriv2_into,
    		"Evaluates the second derivatives into a matrix"							)
    .def("evalSingle_into", 	&Class::evalSingle_into,
    		"Evaluates the values of a single basis function into a matrix"				)
    .def("derivSingle_into", 	&Class::derivSingle_into,
    		"Evaluates the derivatives of a single basis functioninto a matrix" 		)
    .def("deriv2Single_into", 	&Class::deriv2Single_into,
    		"Evaluates the second derivatives of a single basis function into a matrix" )

    // Derived from gsHTensorBasis
    .def("size",&Class::size,"Returns the domain dimension")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(               ) const> (&Class::support), "Get the support of the basis")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(const index_t &) const> (&Class::support), "Get the support of the basis function with an index i")
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
    	py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::refine  ), "Refines the basis given a box")
    .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::unrefine), "Refines the basis given a box")
    .def("refine  ", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")
    // .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
    .def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
    .def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
    .def("asElements        ",&Class::asElements        ,"Returns the elements given refinement boxes")
    .def("asElementsUnrefine",&Class::asElementsUnrefine,"Returns the elements given refinement boxes")
    ;
}

void pybind11_init_gsHBSpline2(py::module &m)
{
	using Class = gsHBSpline<2,real_t>;
	py::class_<Class>(m, "gsHBSpline2")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<2,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<2,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<2,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

void pybind11_init_gsHBSpline3(py::module &m)
{
	using Class = gsHBSpline<3,real_t>;
	py::class_<Class>(m, "gsHBSpline3")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<3,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<3,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<3,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

void pybind11_init_gsHBSpline4(py::module &m)
{
	using Class = gsHBSpline<4,real_t>;
	py::class_<Class>(m, "gsHBSpline4")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<4,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<4,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<4,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

#endif

}
