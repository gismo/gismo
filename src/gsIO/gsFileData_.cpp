#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsFileData.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsFileData<real_t>;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;
  void pybind11_init_gsFileData(py::module &m) {

    py::class_<gsFileData<real_t>> fd(m, "gsFileData");

    fd.def(py::init<>())
      .def(py::init<const std::string&>())
      
      .def("read", &gsFileData<>::read)
      .def("clear", &gsFileData<>::clear)
      .def("numData", &gsFileData<>::numData)
      .def("save", (void (gsFileData<>::*) ())                         &gsFileData<>::save)
      .def("save", (void (gsFileData<>::*) (const std::string&))       &gsFileData<>::save)
      .def("save", (void (gsFileData<>::*) (const std::string&, bool)) &gsFileData<>::save)
      .def("saveCompressed", (void (gsFileData<>::*) ())
           &gsFileData<>::saveCompressed)
      .def("saveCompressed", (void (gsFileData<>::*) (const std::string&))
           &gsFileData<>::saveCompressed)
      .def("dump", (void (gsFileData<>::*) ())                         &gsFileData<>::dump)
      .def("dump", (void (gsFileData<>::*) (const std::string&))       &gsFileData<>::dump)
      .def("addComment", &gsFileData<>::addComment)
      .def("lastPath", &gsFileData<>::lastPath)
      .def("setFloatPrecision", &gsFileData<>::setFloatPrecision)
      .def("getFloatPrecision", &gsFileData<>::getFloatPrecision)

      .def("bufferSize", &gsFileData<>::bufferSize)
      .def("print", &gsFileData<>::print)
      .def("contents", &gsFileData<>::contents)
      .def("numTags", &gsFileData<>::numTags)
      
      //.def("getId", (void (gsFileData<real_t>::*)(const int&, Object&)) &gsFileData<>::getId<Object>)

      .def("__repr__",
           [](const gsFileData<> &obj) {
             std::stringstream os;
             os << obj;
             return os.str();
           })
      ;    
  }
  
#endif // GISMO_BUILD_PYBIND11
  
} // end namespace gismo
