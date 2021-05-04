/* Symbol export for G+Smo shared object */

#define gsKnotVector_EXPORT

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp> //dependency

namespace gismo
{

//CLASS_TEMPLATE_INST gsKnotVector<real_t>;
CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;

#ifdef GISMO_BUILD_PYBIND11
template<typename T>
void pybind11_init_gsKnotVector(py::module &m)
{
  typedef index_t mult_t;
  using Class = gsKnotVector<real_t>;
  py::class_<Class>(m, "gsKnotVector")

    // Constructors
    // Empty constructor
    .def(py::init<>()) // How to set the degree to -1?
    
    .def(py::init<std::vector<real_t>, short_t>()) // knot container

    // Member functions
    .def("size", &Class::size, "Returns the KnotVector number of elements")
    //There are two version of the function insert, so I need to declare its arguments.
    //.def("insert", (void (Class::*)(T, mult_t)) &Class::insert) -> this does not work, I broke pygismo
    ;
}
#endif


}
