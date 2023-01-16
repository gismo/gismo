/* Symbol export for G+Smo shared object */

#define gsKnotVector_EXPORT

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp> //dependency

namespace gismo
{

//CLASS_TEMPLATE_INST gsKnotVector<real_t>;
CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;

#ifdef GISMO_WITH_PYBIND11
void pybind11_init_gsKnotVector(py::module &m)
{
  typedef gsKnotVector<real_t>::mult_t mult_t;
  using Class = gsKnotVector<real_t>;
  py::class_<Class>(m, "gsKnotVector")

    // Constructors
    // Empty constructor
    .def(py::init<>()) // How to set the degree to -1?
    
    .def(py::init<std::vector<real_t>, short_t>()) // knot container

    // Member functions
    .def("get", &Class::get, "Returns the knot vector data")
    .def("degree", &Class::degree, "Returns the degree of the knot vector")
    .def("size", &Class::size, "Returns the KnotVector number of knots including repetitions")
    .def("uSize", &Class::uSize, "Returns the KnotVector number of knots without repetitions")
    .def("knot", &Class::operator [], "Returns the i-th knot")
    //.def("asArray", (const std::vector<real_t>& (Class::*)()) &Class::operator (), "Returns the array of knots")
    // operator const knotContainer& () const {return m_repKnots;}
    .def("uValue", &Class::uValue, "Returns the value of the i-th knot")
    .def("numElements", &Class::numElements, "Returns the number of knot intervals inside the domain")
    .def("multiplicities", &Class::multiplicities, "Returns vector of multiplicities of the knots")
    //There are two version of the function insert, so I need to declare its arguments.
    .def("insert", (void (Class::*)(real_t, mult_t)) &Class::insert)
    .def("uFind", &Class::uFind, "Returns poiter to the knot at the beginning of the _knot interval_ containing the knot")
    .def("iFind", &Class::iFind, "Returns pointer to the last occurrence of the knot at the beginning of the _knot interval_ containing the knot")
    .def("first", &Class::first, "Returns the first knot")
    .def("last", &Class::last, "Returns the last knot")
    .def("check", &Class::check, "Checks whether the knot vector is in a consistent state")
    .def("isConsistent", &Class::isConsistent, "Sanity check")
    .def("inDomain", &Class::inDomain, "Checks, whether the given value is inside the domain")
    .def("greville",static_cast<gsMatrix<real_t>
         (Class::*)(void) const>(&Class::greville), "Returns the Greville points")
    ;
}
#endif


}
