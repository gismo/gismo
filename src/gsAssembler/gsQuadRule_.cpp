#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsQuadRule.hpp>

namespace gismo
{

    CLASS_TEMPLATE_INST gsQuadRule<real_t> ;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsQuadRule(py::module &m)
{
    using Class = gsQuadRule<real_t>;
    py::class_<Class>(m, "gsQuadRule")
    .def("referenceNodes", &Class::referenceNodes, "Returns the reference nodes")
    .def("referenceWeights", &Class::referenceWeights, "Returns the reference weights")
    .def("numNodes", &Class::numNodes, "Returns the number of nodes")
    .def("dim", &Class::dim, "Returns the dimension of the quadrature rule")
    .def("mapTo", static_cast<void (Class::*)(const gsVector<real_t> &, const gsVector<real_t> &, gsMatrix<real_t> &, gsVector<real_t> &) const> (&Class::mapTo), "Maps quadrature rule (i.e., points and weights) from the reference domain to an element.")
    .def("mapTo", static_cast<void (Class::*)(real_t, real_t, gsMatrix<real_t> &, gsVector<real_t> &) const> (&Class::mapTo), "Maps a univariate quadrature rule (i.e., points and weights) from the reference interval to an arbitrary interval.")
    ;
}
#endif

}
