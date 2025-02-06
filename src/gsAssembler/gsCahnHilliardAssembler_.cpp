#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsCahnHilliardAssembler.h>
#include <gsAssembler/gsCahnHilliardAssembler.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsCahnHilliardAssembler<real_t>;


#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsCahnHilliardAssembler(py::module &m)
  {
    using Class = gsCahnHilliardAssembler<real_t>;
    py::class_<Class>(m, "gsCahnHilliardAssembler")

    // Constructors
    .def(py::init<gsMultiPatch<real_t>&,
                  gsMultiBasis<real_t>&,
                  gsBoundaryConditions<real_t> & >())

    // Member functions
    .def("initialize", &Class::initialize, "Initializes the assembler")
    .def("options", &Class::options, "Returns the options")
    .def("setOptions", &Class::setOptions, "Sets the options, ignores unknown options")
    .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")

    .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")
    .def("assembleResidual", &Class::assembleResidual, "Assembles the residual")
    .def("assembleJacobian", &Class::assembleJacobian, "Assembles the Jacobian matrix")
    .def("assembleJacobian", &Class::assembleJacobian, "Assembles the Jacobian matrix")
    .def("assembleNitscheVector", &Class::assembleNitscheVector, "Assembles the vector containing Nitsche boundary terms on weakly clamped boundaries")
    .def("assembleNitscheMatrix", &Class::assembleNitscheMatrix, "Assembles the matrix containing Nitsche boundary terms on weakly clamped boundaries")
    .def("assembleMassMatrix", &Class::assembleMassMatrix, "Assembles the mass matrix")

    .def("matrix_into", &Class::matrix_into, "Returns the matrix")
    .def("rhs_into"   , &Class::rhs_into,    "Returns the RHS vector")

    .def("constructSolution", static_cast<void (Class::*)(gsMatrix<real_t> &, gsMultiPatch<real_t> &)>(&Class::constructSolution), "Constructs a spline solution from a solution vector")
    .def("constructSolution", static_cast<void (Class::*)(gsMatrix<real_t> &, gsMappedSpline<2,real_t> &)>(&Class::constructSolution), "Constructs a mapped spline solution from a solution vector")
    .def("constructSolution", static_cast<void (Class::*)(const gsMultiPatch<real_t> &, gsMatrix<real_t> &)>(&Class::constructSolution), "Constructs a solution vector from a multi-patch solution")
    ;
  }

#endif
}

