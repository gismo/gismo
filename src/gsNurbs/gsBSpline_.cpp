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
  using Class = gsBSpline<real_t>;
  py::class_<Class>(m, "gsBSpline")

    // Constructors
    .def(py::init<real_t,real_t,unsigned, int, gsMatrix<real_t>,unsigned,bool>())

    // Member functions
    .def("degree", &Class::degree, "Returns the degree of the B-Spline") //needs default argument i=0
    .def("domainDim", &Class::domainDim, "Returns the parametric dimension of the B-Spline")
    .def("targetDim", &Class::targetDim, "Returns the target dimension of the B-Spline")
    .def("coefDim", &Class::coefDim, "Returns the number of coefficients defining this B-Spline")
      .def("coefs", (gsMatrix<real_t>& (Class::*)())&Class::coefs,
           "Returns the coeffcient matrix (as a reference)") //there are 2 versions of coefs()
    .def("numCoefs", &Class::numCoefs, "Returns the number of coefficients")
    .def("sample", &Class::sample, "Returns samples on the Bspline curve")
    .def("eval", &Class::eval, "Returns the evaluation of the Bspline curve on the input")
    .def("eval_into", &Class::eval_into, "Evaluation of the Bspline curve on the input")
    //define eval(..)   // This is defined in gsCore/gsGeometry.h
    ;
}

#endif

}
