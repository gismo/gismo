#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBSpline<real_t>;

CLASS_TEMPLATE_INST internal::gsXml<gsBSpline<real_t> >;


#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBSpline(py::module &m)
{
  using Base  = gsGeometry<real_t>;
  using Class = gsBSpline<real_t>;
  py::class_<Class,Base>(m, "gsBSpline")

    // Constructors
    .def(py::init<real_t,real_t,unsigned, int, gsMatrix<real_t>,unsigned,bool>())
    .def(py::init<gsKnotVector<real_t>, gsMatrix<real_t> >() )

    // Member functions
    .def("degree", &Class::degree, "Returns the degree of the B-Spline") //needs default argument i=0
    .def("insertKnot", &Class::insertKnot, "Insert a knot")
    .def("degreeElevate", &Class::degreeElevate, "Elevate the degree")
    .def("domainDim", &Class::domainDim, "Returns the parametric dimension of the B-Spline")
    .def("targetDim", &Class::targetDim, "Returns the target dimension of the B-Spline")
    .def("coefDim", &Class::coefDim, "Returns the number of coefficients defining this B-Spline")
    .def("coefs", (gsMatrix<real_t>& (Class::*)())&Class::coefs,
           py::return_value_policy::reference_internal,
           "Returns the coeffcient matrix (as a reference)") //there are 2 versions of coefs()
    .def("knots", static_cast<gsKnotVector<real_t>& (Class::*)(int)> (&Class::knots), "Set the pet's age")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Set the pet's age")
    .def("numCoefs", &Class::numCoefs, "Returns the number of coefficients")
    .def("sample", &Class::sample, "Returns samples on the Bspline curve")
    .def("eval", &Class::eval, "Returns the evaluation of the Bspline curve at the input point")
    .def("eval_into", &Class::eval_into, "Evaluation of the Bspline curve at the input point")
    .def("deriv", &Class::deriv, "Returns the derivative of the Bspline curve at the input point")
    .def("deriv_into", &Class::deriv_into, "Derivative of the Bspline curve at the input point")
    .def("deriv2", &Class::deriv2, "Returns the second derivative of the Bspline curve at the input point")
    ;
}

#endif

}
