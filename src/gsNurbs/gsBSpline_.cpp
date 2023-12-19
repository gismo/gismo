#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBSpline<real_t>;

CLASS_TEMPLATE_INST internal::gsXml<gsBSpline<real_t> >;


#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBSpline(py::module &m)
{
  using Base  = gsGeometry<real_t>;
  using Class = gsBSpline<real_t>;
  py::class_<Class,Base>(m, "gsBSpline")

    // Constructors
      .def(py::init<real_t,real_t,unsigned, int, gsMatrix<real_t>,unsigned,bool>(),
           py::arg("leftKnot"), py::arg("rightKnot"), py::arg("coefs"),
           py::arg("interior"), py::arg("degree"),
           py::arg("mult_interior") = 1, py::arg("periodic") = false)
      .def(py::init<gsBSplineBasis<real_t>, gsMatrix<real_t>>(),
           py::arg("basis"), py::arg("coefs") )
      .def(py::init<gsKnotVector<real_t>, gsMatrix<real_t>, bool>(),
           py::arg("knotvector"), py::arg("coefs"), py::arg("periodic")=false)

    // Member functions
    .def("degree", &Class::degree,py::arg("direction") = 0, "Returns the degree of the B-Spline")
    .def("insertKnot", &Class::insertKnot, "Insert a knot")
    .def("degreeElevate", &Class::degreeElevate, "Elevate the degree")
    .def("coefDim", &Class::coefDim, "Returns the number of coefficients defining this B-Spline")
    .def("knots", static_cast<gsKnotVector<real_t>& (Class::*)(int)> (&Class::knots),
         py::arg("direction") = 0,
         py::return_value_policy::reference_internal,
          "Gets the knots")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots),
         py::arg("direction") = 0,
         py::return_value_policy::reference_internal,
         "Gets the knots as a reference")
    .def("domainStart", &Class::domainStart, "Returns the left end of the domain")
    .def("domainEnd", &Class::domainEnd, "Returns the right end of the domain")
    .def("numCoefs", &Class::numCoefs, "Returns the number of coefficients")
    .def("sample", &Class::sample, "Returns samples on the Bspline curve")
    ;
}

#endif

}
