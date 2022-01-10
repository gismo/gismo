
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometry.hpp>

#include <gsCore/gsCurve.h>
#include <gsCore/gsCurve.hpp>

#include <gsCore/gsSurface.h>
#include <gsCore/gsSurface.hpp>

#include <gsCore/gsVolume.h>
#include <gsCore/gsVolume.hpp>

#include <gsCore/gsBulk.h>
#include <gsCore/gsBulk.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsGeometry<real_t>;
    CLASS_TEMPLATE_INST gsCurve   <real_t> ;
    CLASS_TEMPLATE_INST gsSurface <real_t> ;
    CLASS_TEMPLATE_INST gsVolume  <real_t> ;
    CLASS_TEMPLATE_INST gsBulk    <real_t> ;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsGeometry(py::module &m)
{
  using Class = gsGeometry<real_t>;
  py::class_<Class>(m, "gsGeometry")

  // Member functions
  .def("domainDim", &Class::domainDim, "Gives the domain dimension")
  .def("targetDim", &Class::targetDim, "Gives the target dimension")
  .def("parDim", &Class::targetDim, "Gives the parameter dimension")
  .def("geoDim", &Class::targetDim, "Gives the geometry dimension")

  .def("eval", &Class::eval, "Evaluates points into a matrix")
  .def("eval_into", &Class::eval_into, "Evaluates points into a matrix")
  .def("coefs", static_cast<      gsMatrix<real_t>& (Class::*)()      > (&Class::coefs), "Get the coefficients as a reference")
  .def("coefs", static_cast<const gsMatrix<real_t>& (Class::*)() const> (&Class::coefs), "Get the coefficients as a const reference")
  ;
}

#endif

}
