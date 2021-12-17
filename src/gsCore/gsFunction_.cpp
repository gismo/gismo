#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunction.h>
#include <gsCore/gsFunction.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsFunction<real_t> ;

#ifdef GISMO_BUILD_PYBIND11  

void pybind11_init_gsFunction(py::module &m)
{
  using Class = gsFunction<real_t>;
  py::class_<Class>(m, "gsFunction")
  
  //Constructor?
  
  //Member functions
  .def("deriv_into", &Class::deriv_into, "Returns the first derivatives")
  ;
}
#endif

}
