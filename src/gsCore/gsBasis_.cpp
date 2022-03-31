
#include <gsCore/gsDebug.h> // to purge warning on MSVC
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsBasis.h>
#include <gsCore/gsBasis.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBasis<real_t>;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBasis(py::module &m)
{
  using Class = gsBasis<real_t>;
  py::class_<Class>(m, "gsBasis")

  // Member functions
  .def("eval", &Class::eval, "Evaluates points into a matrix")
  .def("anchors", &Class::anchors, "Returns the anchor points of the basis")
  .def("collocationMatrix", &Class::collocationMatrix, "Computes a (sparse) collocation matrix")
  .def("numElements", static_cast<size_t (Class::*)(boxSide const & ) const> ( &Class::numElements), "Number of elements")
  .def("numElements", static_cast<size_t (Class::*)() const> ( &Class::numElements), "Number of elements")
  .def("component", static_cast<gsBasis<real_t> & (Class::*)(short_t ) > ( &Class::component), "Return the basis of component",py::return_value_policy::reference)
  ;
}

#endif

}
