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
    using Class = gsMultiPatch<real_t>;
    py::class_<Class>(m, "gsMultiPatch")

      // Constructors
      .def(py::init<>())

      // Member functions
      .def("domainDim", &Class::domainDim, "Returns the domain dimension of the multipatch")
      .def("targetDim", &Class::targetDim, "Returns the target dimension of the multipatch")
      .def("nPatches", &Class::nPatches, "Returns the number of patches stored in the multipatch")
      .def("patch", &Class::patch, "Access the a patch of the multipatch")

      ;
  }

#endif

}
