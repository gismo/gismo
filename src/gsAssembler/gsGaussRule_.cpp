#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsAssembler/gsGaussRule.hpp>

namespace gismo
{

    CLASS_TEMPLATE_INST gsGaussRule<real_t> ;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsGaussRule(py::module &m)
{
    using Base = gsQuadRule<real_t>;
    using Class = gsGaussRule<real_t>;
    py::class_<Class,Base>(m, "gsGaussRule")
    // Constructors
    .def(py::init<const gsBasis<real_t> &, const real_t, const index_t, short_t>(),
         "Initialize a tensor-product Gauss quadrature rule for basis using quA *deg_i + quB nodes (direction-wise)",
         py::arg("basis"), py::arg("quA"), py::arg("quB"), py::arg("fixDir") = -1)
    .def(py::init<const gsBasis<real_t> &, const gsOptionList &, short_t>(),
         "Initialize a tensor-product Gauss quadrature rule for basis using quA *deg_i + quB nodes (direction-wise). Values of quA and quB are taken from the options",
         py::arg("basis"), py::arg("options"), py::arg("fixDir") = -1)
    ;
}
#endif

}
