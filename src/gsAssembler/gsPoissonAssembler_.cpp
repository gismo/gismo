
#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsPoissonAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsPoissonAssembler<real_t>;

#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsPoissonAssembler(py::module &m)
    {
        using Base  = gsAssembler<real_t>;
        using Class = gsPoissonAssembler<real_t>;
        py::class_<Class, Base>(m, "gsPoissonAssembler")
                // Constructors
                .def(py::init<gsMultiPatch<real_t> const &,
                        gsMultiBasis<real_t> const &,
                        gsBoundaryConditions<real_t> const &,
                        const gsFunction<real_t> &,
                        dirichlet::strategy,
                        iFace::strategy>(),
                //Member functions
                     py::arg("patches"), py::arg("basis"), py::arg("bconditions"),
                     py::arg("rhs"), py::arg("dirStrategy") = 11,
                     py::arg("intstrategy") = 1)
                ;
    }

#endif
}
