#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsFunctionExpr.hpp>


namespace gismo
{

CLASS_TEMPLATE_INST gsFunctionExpr<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsFunctionExpr<real_t> >;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsFunctionExpr(py::module &m)
{
  using Base  = gsFunction<real_t>;
  using Class = gsFunctionExpr<real_t>;
  py::class_<Class,Base>(m, "gsFunctionExpr")

    // Constructors
    .def(py::init<std::string,short_t>())
    .def(py::init<std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,std::string,std::string,std::string,std::string,std::string,std::string,short_t>())

    // Member functions
    .def("domainDim", &Class::domainDim, "Returns the parametric dimension of the B-Spline")
    .def("targetDim", &Class::targetDim, "Returns the target dimension of the B-Spline")
    .def("eval", &Class::eval, "Returns the evaluation of the Bspline curve on the input")
    .def("eval_into", &Class::eval_into, "Evaluation of the Bspline curve on the input")
    //define eval(..)   // This is defined in gsCore/gsGeometry.h
    .def("deriv", &Class::deriv, "Evaluates the first derivatives of the B-spline curve")
    .def("deriv_into", &Class::deriv_into, "Evaluates the first derivatives of the B-spline curve")
    .def("deriv2", &Class::deriv2, "Evaluated the second derivatives of the B-spline curve")
    .def("deriv2_into", &Class::deriv2_into, "Evaluates the second derivatives of the B-spline curve")
    ;
}

#endif

}
