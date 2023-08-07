#include <gsAssembler/gsExprAssembler.h>

namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsExprAssembler(py::module &m)
    {
        using Class = gsExprAssembler<real_t>;
        py::class_<Class>(m, "gsExprAssembler")
                // Constructors
                .def(py::init<index_t, index_t>(),
                        //Member functions
                     py::arg("rBlocks") = 1, py::arg("cBlocks") = 1)
                 .def("setOption", &Class::setOptions, "Set the assembler options")
                 .def("getMap", &Class::getMap, "Registers g as an isogeometric geometry map and return a "
                                                "handle to it")
                 .def("getSpace", &Class::getSpace, "Registers mp as an isogeometric (both trial and test)"
                                                    " space and returns a handle to it")
                 .def("getCoeff",
                      static_cast<expr::gsComposition<real_t> (Class::*)(const gsFunctionSet<real_t> &, Class::geometryMap &)> (&Class::getCoeff),
                      "Registers func as a variable and returns a handle to it")
                .def("getCoeff",
                     static_cast<Class::variable (Class::*)(const gsFunctionSet<real_t> &)> (&Class::getCoeff),
                     "Registers func as a variable and returns a handle to it")
                 .def("getSolution",
                      &Class::getSolution,
                      "Registers a representation of a solution variable from space s, based on the vector cf.")
                  .def("numDofs", &Class::numDofs, "Number of degrees of freedom")
                  //.def("assemble", &Class::assemble<real_t, const expr >, "Adds the expressions args to the system matrix/rhs")
                  //.def("assembleBdr", &Class::assembleBdr, "Assemble boundary terms")
                  .def("matrix", &Class::matrix, "Returns the assembled matrix")
                .def("rhs", &Class::rhs, "Returns the assembled right-hand-side")
                ;
    }
    //.def("constructSolution", static_cast<void (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &, short_t) const> (&Class::constructSolution), "Constuct the solution function")
#endif

}