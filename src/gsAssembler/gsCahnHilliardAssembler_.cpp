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
    .def(py::init<const gsMultiPatch<real_t>&,
                  const gsMultiBasis<real_t>&,
                  const gsBoundaryConditions<real_t> & >())

    // Member functions
    .def("initialize", &Class::initialize, "Initializes the assembler")
    .def("options", &Class::options, py::return_value_policy::reference_internal, "Returns the options (reference)")
    .def("setOptions", &Class::setOptions, "Sets the options, ignores unknown options")
    .def("numDofs", &Class::numDofs, "Returns the number of degrees of freedom of the system")

    .def("setSpaceBasis", &Class::setSpaceBasis, "Sets the basis used for discretization (but not for quadrature)")

    .def("assemble", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; self.assemble(C); }, "Assembles the full system for given coefficient vectors")
    .def("assemble", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; self.assemble(C); }, "Assembles the full system for given function sets")
    .def("assembleResidual", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; self.assembleResidual(C); }, "Assembles the spatial residual vector for given coefficient vector (caller adds M*dC)")
    .def("assembleResidual", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; self.assembleResidual(C); }, "Assembles the spatial residual vector for given function set (caller adds M*dC)")
    .def("assembleJacobian", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; self.assembleJacobian(C); }, "Assembles the Jacobian matrix for given coefficient vector", py::arg("C"))
    .def("assembleJacobian", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; self.assembleJacobian(C); }, "Assembles the Jacobian matrix for given function set", py::arg("C"))
    .def("assembleNitscheVector", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; self.assembleNitscheVector(C); }, "Assembles the Nitsche boundary vector for given coefficient vector", py::arg("C"))
    .def("assembleNitscheVector", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; self.assembleNitscheVector(C); }, "Assembles the Nitsche boundary vector for given function set", py::arg("C"))
    .def("assembleNitscheMatrix", [](Class &self) { py::gil_scoped_release release; self.assembleNitscheMatrix(); }, "Assembles the matrix containing Nitsche boundary terms on weakly clamped boundaries")
    .def("assembleMassMatrix",    [](Class &self) { py::gil_scoped_release release; self.assembleMassMatrix(); }, "Assembles the mass matrix")
    
    .def("matrix_into", &Class::matrix_into, "Returns the matrix")
    .def("rhs_into"   , &Class::rhs_into,    "Returns the RHS vector")

    .def("rhs", &Class::rhs, "Returns the RHS vector as a numpy array (convenience function)")

    .def("constructSolution", static_cast<void (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &) const>(&Class::constructSolution), "Constructs the solution (to a multi-patch) from the vector of coefficients", py::arg("C"), py::arg("solution"))
    .def("constructSolution", static_cast<void (Class::*)(const gsMatrix<real_t> &, gsMappedSpline<2,real_t> &) const>(&Class::constructSolution), "Constructs the solution (to a mapped spline) from the vector of coefficients", py::arg("C"), py::arg("solution"))
    .def("constructSolution", [](Class & self, const gsMultiPatch<real_t> &C) { gsMatrix<real_t> Cmat; self.constructSolution(C,Cmat); return Cmat; }, "Constructs the vector of coefficients from the solution (multi-patch)", py::arg("solution"))
    
    .def("computeMass", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; return self.computeMass(C); }, "Computes the mass M(c) = ∫ c dx", py::arg("C"))
    .def("computeMass", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; return self.computeMass(C); }, "Computes the mass M(c) = ∫ c dx from a coefficient vector", py::arg("C"))
    .def("computeDissipation", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; return self.computeDissipation(C); }, "Computes the dissipation D(c) = ∫ |∇mu|^2 dx", py::arg("C"))
    .def("computeDissipation", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; return self.computeDissipation(C); }, "Computes the dissipation D(c) = ∫ |∇mu|^2 dx from a coefficient vector", py::arg("C"))
    .def("computeDissipation", [](Class &self, const gsFunctionSet<real_t> &C, const gsFunctionSet<real_t> &mu) { py::gil_scoped_release release; return self.computeDissipation(C, mu); }, "Computes the dissipation D(c) = ∫ M(c) |∇mu|^2 dx where mu = -c + c^3 - lambda Δc", py::arg("C"), py::arg("mu"))
    .def("computeDissipation", [](Class &self, const gsMatrix<real_t> &C, const gsMatrix<real_t> &mu) { py::gil_scoped_release release; return self.computeDissipation(C, mu); }, "Computes the dissipation D(c) = ∫ M(c) |∇mu|^2 dx where mu = -c + c^3 - lambda Δc from coefficient vectors", py::arg("C"), py::arg("mu"))
    .def("computeEnergy", [](Class &self, const gsFunctionSet<real_t> &C) { py::gil_scoped_release release; return self.computeEnergy(C); }, "Computes the energy E(c) = ∫ [ 1/4 (c^2 - 1)^2 + lambda/2 |∇c|^2 ] dx", py::arg("C"))
    .def("computeEnergy", [](Class &self, const gsMatrix<real_t> &C) { py::gil_scoped_release release; return self.computeEnergy(C); }, "Computes the energy E(c) = ∫ [ 1/4 (c^2 - 1)^2 + lambda/2 |∇c|^2 ] dx from a coefficient vector", py::arg("C"))
    ;
  }

#endif
}

