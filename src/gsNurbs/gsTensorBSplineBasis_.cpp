/* Symbol export for G+Smo shared object */

//#define gsTensorBSplineBasis_EXPORT

#include <gsCore/gsBasisFun.h>
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

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;
// Helper binder to add dimension-specific constructors without
// producing hard template substitution errors in the function body.
template<short_t D, typename C, typename B>
struct tensor_ctor_binder { static void bind(py::class_<C,B>&) {} };

template<typename C, typename B>
struct tensor_ctor_binder<1,C,B>
{ static void bind(py::class_<C,B>& c)
    { c.def(py::init<gsKnotVector<real_t>>(), "Constructor for 1D tensor product B-spline basis"); } };

template<typename C, typename B>
struct tensor_ctor_binder<2,C,B>
{ static void bind(py::class_<C,B>& c)
    { c.def(py::init<gsKnotVector<real_t>, gsKnotVector<real_t>>(), "Constructor for 2D tensor product B-spline basis"); } };

template<typename C, typename B>
struct tensor_ctor_binder<3,C,B>
{ static void bind(py::class_<C,B>& c)
    { c.def(py::init<gsKnotVector<real_t>, gsKnotVector<real_t>, gsKnotVector<real_t>>(), "Constructor for 3D tensor product B-spline basis"); } };

template<typename C, typename B>
struct tensor_ctor_binder<4,C,B>
{ static void bind(py::class_<C,B>& c)
    { c.def(py::init<gsKnotVector<real_t>, gsKnotVector<real_t>, gsKnotVector<real_t>, gsKnotVector<real_t>>(), "Constructor for 4D tensor product B-spline basis"); } };

template <short_t d>
void pybind11_init_gsTensorBSplineBasis(py::module &m)
{
        using Base = gsBasis<real_t>;
        using Class = gsTensorBSplineBasis<d,real_t>;

        py::class_<Class,Base> cls(m, ("gsTensorBSplineBasis" + std::to_string(d)).c_str());

        cls.def(py::init<>());

        // add the appropriate constructors for this dimension
        tensor_ctor_binder<d,Class,Base>::bind(cls);

        cls
        .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
        .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
        .def("component", static_cast<      gsBSplineBasis<real_t>& (Class::*)(int)      > (&Class::component), "Returns the basis component as a reference")
        .def("component", static_cast<const gsBSplineBasis<real_t>& (Class::*)(int) const> (&Class::component), "Returns the basis component as a const reference")
        .def("size", static_cast<index_t (Class::*)() const> (&Class::size), "Returns the size")
        .def("dim", &Class::dim, "Returns the dimension")
        .def("eval", &Class::eval, "Evaluates points into a matrix")
        .def("function", &Class::function, "Returns the basis function i")
        .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &) const> (&Class::evalSingle), "Evaluates the basis function i")
        .def("degree", &Class::degree, "Returns the degree")
        ;
}

template void pybind11_init_gsTensorBSplineBasis<2>(py::module &m);
template void pybind11_init_gsTensorBSplineBasis<3>(py::module &m);
template void pybind11_init_gsTensorBSplineBasis<4>(py::module &m);

#endif

}
