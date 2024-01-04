#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsBiharmonicExprAssembler.h>
#include <gsAssembler/gsBiharmonicExprAssembler.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsBiharmonicExprAssembler<real_t>;


#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsBiharmonicExprAssembler(py::module &m)
  {
    using Class = gsBiharmonicExprAssembler<real_t>;
    py::class_<Class>(m, "gsBiharmonicExprAssembler")

    // Constructors
    .def(py::init<gsMultiPatch<real_t>&,
                  gsMultiBasis<real_t>&,
                  gsFunctionSet<real_t> &,
                  gsBoundaryConditions<real_t> & >())



    // Member functions
    .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")
    .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")
    .def("assemble", &Class::assemble, "Assembles the linear system")

    .def("matrix", &Class::matrix, "Returns the matrix")
    .def("rhs"   , &Class::rhs,    "Returns the RHS vector")

    .def("penalty", &Class::penalty, "Returns the penalty factor for Nitsche's method")

    .def("l2error", &Class::l2error, "Returns the L2 error between the discrete solution and the exact solution")
    .def("h1error", &Class::h1error, "Returns the H1 error between the discrete solution and the exact solution")
    .def("h2error", &Class::h2error, "Returns the H2 error between the discrete solution and the exact solution")
    .def("errors", &Class::errors, "Returns the L2, H1 and H2 errors between the discrete solution and the exact solution")
    .def("interfaceError", &Class::interfaceError, "Returns the interface error between the discrete solution and the exact solution")

    .def("options", &Class::options, "Returns the options")
    .def("setOptions", &Class::setOptions, "Sets the options, ignores unknown options")

    ;
  }

#endif
}

