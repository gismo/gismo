#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsReadFile.h>

#ifdef GISMO_BUILD_PYBIND11
#include <gsCore/gsMultiPatch.h>
#endif


namespace gismo
{

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;
  void pybind11_init_gsReadFile(py::module &m) {

    using Class = gsReadFile<real_t>;
    py::class_<Class> fd(m, "gsReadFile");

    fd.def(py::init<const std::string&, gsMultiPatch<real_t>&>())
      // .def(py::init<const std::string&, gsTensorBSpline<1,real_t>>())
      // .def(py::init<const std::string&, gsTensorBSpline<2,real_t>>())
      // .def(py::init<const std::string&, gsTensorBSpline<3,real_t>>())
      
      // .def("read", &Class::read)
      ;    
  }
  
#endif // GISMO_BUILD_PYBIND11
  
} // end namespace gismo
