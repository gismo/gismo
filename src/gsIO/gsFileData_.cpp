#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsFileData.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsFileData<real_t>;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;
  void pybind11_init_gsFileData(py::module &m) {

    using Class = gsFileData<real_t>;
    py::class_<Class> fd(m, "gsFileData");

    fd.def(py::init<>())
      .def(py::init<const std::string&>())
      
      .def("read", &Class::read)
      .def("clear", &Class::clear)
      .def("numData", &Class::numData)
      .def("save", (void (Class::*) ())                         &Class::save)
      .def("save", (void (Class::*) (const std::string&))       &Class::save)
      .def("save", (void (Class::*) (const std::string&, bool)) &Class::save)
      .def("saveCompressed", (void (Class::*) ())
           &Class::saveCompressed)
      .def("saveCompressed", (void (Class::*) (const std::string&))
           &Class::saveCompressed)
      .def("dump", (void (Class::*) ())                         &Class::dump)
      .def("dump", (void (Class::*) (const std::string&))       &Class::dump)
      .def("addComment", &Class::addComment)
      .def("lastPath", &Class::lastPath)
      .def("setFloatPrecision", &Class::setFloatPrecision)
      .def("getFloatPrecision", &Class::getFloatPrecision)

      // .def("getId", static_cast<const gsBasis<real_t> & (Class::*)(const size_t) const > (&Class::getId))
      // .def("getId", static_cast<void (Class::*)(const int &, gsMultiPatch<real_t>) const > (&Class::getId<gsMultiPatch<real_t>), "Gets a const reference to basis with index i")


      .def("bufferSize", &Class::bufferSize)
      .def("print", &Class::print)
      .def("contents", &Class::contents)
      .def("numTags", &Class::numTags)
      
      // .def("getId", (void (gsFileData<real_t>::*)(const int&, Object&)) &Class::getId<Object>)

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
