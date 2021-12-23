/* Symbol export for G+Smo shared object */

//#define gsTensorBSplineBasis_EXPORT

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsTensorBSplineBasis<2,real_t>;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<3,real_t>;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<4,real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t> >;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsTensorBSplineBasis2(py::module &m)
{
  using Base = gsBasis<real_t>;
  using Class = gsTensorBSplineBasis<2,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSplineBasis2")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,    gsKnotVector<real_t>    >())
    .def(py::init<gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *>())

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
    ;
}

void pybind11_init_gsTensorBSplineBasis3(py::module &m)
{
  using Base = gsBasis<real_t>;
  using Class = gsTensorBSplineBasis<3,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSplineBasis3")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,    gsKnotVector<real_t>,    gsKnotVector<real_t>    >())
    .def(py::init<gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *>())

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
    ;
}

void pybind11_init_gsTensorBSplineBasis4(py::module &m)
{
  using Base = gsBasis<real_t>;
  using Class = gsTensorBSplineBasis<4,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSplineBasis4")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,    gsKnotVector<real_t>,    gsKnotVector<real_t>,    gsKnotVector<real_t>    >())
    .def(py::init<gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *,gsBSplineBasis<real_t> *>())

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
    ;
}

#endif

}
