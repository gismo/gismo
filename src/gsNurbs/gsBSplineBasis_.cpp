/* Symbol export for G+Smo shared object */

#define gsBSplineBasis_EXPORT

#include <gsIO/gsXml.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.hpp> //dependency (otherwise already included)

namespace gismo
{

// CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t>;
// CLASS_TEMPLATE_INST gsBSplineBasis<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsBSplineBasis<real_t> >;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBSplineBasis(py::module &m)
{
  using Class = gsBSplineBasis<real_t>;
  py::class_<Class>(m, "gsBSplineBasis")

    // Constructors
    .def(py::init<gsKnotVector<real_t>      >())
    .def(py::init<gsKnotVector<real_t>, bool>())
    .def(py::init<real_t,real_t,unsigned, int,unsigned,bool>())

    // Member functions
    .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
	.def("component", static_cast<      gsBSplineBasis<real_t>& (Class::*)(int)      > (&Class::component), "Returns the basis component as a reference")
    .def("component", static_cast<const gsBSplineBasis<real_t>& (Class::*)(int) const> (&Class::component), "Returns the basis component as a const reference")
    // Inherited from gsTensorBasis
    .def("dim", &Class::dim, "Returns the dimension")
    // Inherited from gsBasis
    .def("eval", &Class::eval, "Evaluates points into a matrix")
    // .def("eval_into", &Class::eval_into, "Evaluates points into a matrix")
    .def("numElements", static_cast<size_t (Class::*)() const> (&Class::numElements), "Returns the number of Elements")
    ;
}

#endif

}
