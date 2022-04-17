#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiPatch.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsMultiPatch<real_t> ;

#ifdef GISMO_BUILD_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsMultiPatch(py::module &m)
  {
    using Base1 = gsFunctionSet<real_t>;
    using Base2 = gsBoxTopology;
    using Class = gsMultiPatch<real_t>;
    py::class_<Class,Base1, Base2>(m, "gsMultiPatch")

      // Constructors
      .def(py::init<>())

      // Member functions
      .def("domainDim", &Class::domainDim, "Returns the domain dimension of the multipatch")
      .def("targetDim", &Class::targetDim, "Returns the target dimension of the multipatch")
      .def("nPatches", &Class::nPatches, "Returns the number of patches stored in the multipatch")
      .def("patch", static_cast<gsGeometry<real_t> & (Class::*)(size_t) const> ( &Class::patch), "Access the a patch of the multipatch",py::return_value_policy::reference)
      // Note: Bindings with unique pointers are not possible https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html
      // .def("addPatch", static_cast<void (Class::*)(typename gsGeometry<real_t>::uPtr)> (&Class::addPatch), "Adds a patch")
      .def("addPatch", static_cast<void (Class::*)(   const gsGeometry<real_t> &    )> (&Class::addPatch), "Adds a patch")

      .def("degreeElevate", &Class::degreeElevate, "Elevates the degree", py::arg("elevationSteps")=1)
      .def("uniformRefine", &Class::uniformRefine, "Refines uniformly"  , py::arg("numKnots")=1 , py::arg("mul")=1)

      .def("basis", static_cast<gsBasis<real_t> & (Class::*)(const size_t) const> ( &Class::basis), "Access the i-th basis of the multipatch",py::return_value_policy::reference)

      .def("computeTopology", static_cast<bool (Class::*)(real_t, bool, bool)> ( &Class::computeTopology), "Compute Topology")
      ;
  }

#endif

}
