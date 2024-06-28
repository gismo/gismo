
#include <gsCore/gsDebug.h> // to purge warning on MSVC
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsBasis.h>
#include <gsCore/gsBasis.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBasis<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBasis(py::module &m)
{
    using Base = gsFunctionSet<real_t>;
    using Class = gsBasis<real_t>;
    py::class_<Class,Base>(m, "gsBasis")
    // Member functions
    .def("size",&Class::size,"Returns the size of the basis")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(               ) const> (&Class::support), "Get the support of the basis")
    .def("support", static_cast<gsMatrix<real_t> (Class::*)(const index_t &) const> (&Class::support), "Get the support of the basis function with an index i")

    .def("uniformRefine", static_cast<void (Class::*)(int, int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
        py::arg("numKnots") = 1, py::arg("mul") = 1, py::arg("dir") = -1) //default arguments
    .def("uniformCoarsen", static_cast<void (Class::*)(int)> (&Class::uniformCoarsen), "Refines the basis uniformly",
        py::arg("numKnots") = 1) //default arguments
    .def("degreeElevate", static_cast<void (Class::*)(short_t const &, short_t const)> (&Class::degreeElevate), "Elevate the degree of the basis by the given amount, preserve smoothness",
        py::arg("i") = 1, py::arg("dir") = -1) //default arguments
    .def("degreeReduce" , static_cast<void (Class::*)(short_t const &, short_t const)> (&Class::degreeReduce) , "Elevate the degree of the basis by the given amount, preserve smoothness",
        py::arg("i") = 1, py::arg("dir") = -1) //default arguments
    .def("degreeIncrease", static_cast<void (Class::*)(short_t const &, short_t const)> (&Class::degreeIncrease), "Elevate the degree of the basis by the given amount, preserve smoothness",
        py::arg("i") = 1, py::arg("dir") = -1) //default arguments
    .def("degreeDecrease", static_cast<void (Class::*)(short_t const &, short_t const)> (&Class::degreeDecrease), "Elevate the degree of the basis by the given amount, preserve smoothness",
        py::arg("i") = 1, py::arg("dir") = -1) //default arguments

    .def("refine", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::refine  ), "Refines the basis given a box")
    .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &, int)> (&Class::unrefine), "Refines the basis given a box")
    // .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
    .def("refineElements",&Class::refineElements  ,"Refines the basis given elements  ")
    .def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
    .def("asElements",&Class::asElements        ,"Returns the elements given refinement boxes")
    .def("asElementsUnrefine",&Class::asElementsUnrefine,"Returns the elements given refinement boxes")


    .def("dim", &Class::dim, "Returns the dimension of the basis")
    .def("anchors", &Class::anchors, "Returns the anchor points of the basis")
    .def("collocationMatrix", &Class::collocationMatrix, "Computes a (sparse) collocation matrix")

    .def("uniformRefine", &Class::uniformRefine, "Performs uniform refinement")

    .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &                   ) const> (&Class::evalSingle    ), "Evaluates the basis function i")
    .def("evalSingle_into", static_cast<void        (Class::*)(index_t, const gsMatrix<real_t> &, gsMatrix<real_t>&) const> (&Class::evalSingle_into), "Evaluates the basis function i")

    .def("derivSingle", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &                   ) const> (&Class::derivSingle    ), "Evaluates the derivative of basis function i")
    .def("derivSingle_into", static_cast<void        (Class::*)(index_t, const gsMatrix<real_t> &, gsMatrix<real_t>&) const> (&Class::derivSingle_into), "Evaluates the derivative of basis function i")

    .def("deriv2Single", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &                   ) const> (&Class::deriv2Single    ), "Evaluates the second derivative of basis function i")
    .def("deriv2Single_into", static_cast<void        (Class::*)(index_t, const gsMatrix<real_t> &, gsMatrix<real_t>&) const> (&Class::deriv2Single_into), "Evaluates the second derivative of basis function i")

    .def("numElements", static_cast<size_t (Class::*)(boxSide const & ) const> ( &Class::numElements), "Number of elements")
    .def("numElements", static_cast<size_t (Class::*)() const> ( &Class::numElements), "Number of elements")
    .def("component", static_cast<gsBasis<real_t> & (Class::*)(short_t ) > ( &Class::component), "Return the basis of component",py::return_value_policy::reference)
    ;
}


void pybind11_init_PPN(pybind11::module &m)
{
    pybind11::module ppn = m.def_submodule("ppn");
    // .def("dim", &Class::dim, "Returns the dimension of the basis")
//    ppn.def("collocationMatrix1", &gismo::collocationMatrix1<real_t>, "returns the collocation matrix and its derivatives.");
}

#endif

}
