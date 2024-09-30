
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

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsGeometry(py::module &m)
{
  using Base = gsFunction<real_t>;
  using Class = gsGeometry<real_t>;
  py::class_<Class, Base>(m, "gsGeometry")

  // Member functions
  .def("parDim", &Class::targetDim, "Gives the parameter dimension")
  .def("geoDim", &Class::targetDim, "Gives the geometry dimension")
  .def("coefs", static_cast<      gsMatrix<real_t>& (Class::*)()      > (&Class::coefs),
         py::return_value_policy::reference_internal,
         "Get the coefficients as a reference")
  .def("coefs", static_cast<const gsMatrix<real_t>& (Class::*)() const> (&Class::coefs),
         py::return_value_policy::reference_internal,
         "Get the coefficients as a reference")

  .def("setCoefs", &Class::setCoefs, "Sets the coefficients")
  .def("basis", static_cast<const gsBasis<real_t>& (Class::*)() const>(&Class::basis),
         py::return_value_policy::reference_internal,
         "Returns the bspline basis")
  .def("basis", static_cast<gsBasis<real_t>& (Class::*)()>(&Class::basis),
         py::return_value_policy::reference_internal,
         "Returns the bspline basis as a reference")
  .def("rotate", (void (Class::*)(real_t, const gsVector<real_t,3>&)) &Class::rotate, "Apply 3D Rotation by an angle radians around axis")
  .def("rotate", (void (Class::*)(real_t)) &Class::rotate, "Apply 2D Rotation by an angle radians")
  .def("closestPointTo", (void (Class::*)(real_t)) &Class::rotate, "Get the closest position to a given point in space")

  .def("refineElements",&Class::refineElements  ,"Refines the geometry given elements  ")
  .def("unrefineElements",&Class::unrefineElements,"Unrefines the geometry given elements")

  ;
}

#endif

}
