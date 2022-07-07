
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
    .def("dim", &Class::dim, "Returns the dimension of the basis")
    .def("eval", &Class::eval, "Evaluates points into a matrix")
    .def("anchors", &Class::anchors, "Returns the anchor points of the basis")
    .def("collocationMatrix", &Class::collocationMatrix, "Computes a (sparse) collocation matrix")

    .def("support", static_cast<gsMatrix<real_t> (Class::*)(const index_t&) const> (&Class::support), "Gives the support of basis function i")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)() const> (&Class::support), "Gives the support of the basis")

    .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &                   ) const> (&Class::evalSingle    ), "Evaluates the basis function i")
    .def("evalSingle_into", static_cast<void        (Class::*)(index_t, const gsMatrix<real_t> &, gsMatrix<real_t>&) const> (&Class::evalSingle_into), "Evaluates the basis function i")

    .def("numElements", static_cast<size_t (Class::*)(boxSide const & ) const> ( &Class::numElements), "Number of elements")
    .def("numElements", static_cast<size_t (Class::*)() const> ( &Class::numElements), "Number of elements")
    .def("component", static_cast<gsBasis<real_t> & (Class::*)(short_t ) > ( &Class::component), "Return the basis of component",py::return_value_policy::reference)
    ;
}

#endif

}
