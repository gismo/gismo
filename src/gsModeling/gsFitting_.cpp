#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsFitting.h>
#include <gsModeling/gsFitting.hpp>
#include <gsNurbs/gsBSplineBasis.h>

namespace gismo
{

CLASS_TEMPLATE_INST gsFitting<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;
void pybind11_init_gsFitting(py::module &m)
{
  using Class = gsFitting<real_t>;
  py::class_<Class>(m, "gsFitting")

    // Constructors
    .def( py::init<>() ) // Empty constructor
    .def( py::init<gsMatrix<real_t> const &, gsMatrix<real_t> const &, gsBasis<real_t>&>() )

    // Member functions
    .def("compute", &Class::compute, "Computes the least square fit for a gsBasis.",
          py::arg("lambda") = 0)
    .def("result", &Class::result,
          py::return_value_policy::reference,
          "Returns the result.")
    .def("applySmoothing", &Class::applySmoothing, "apply smoothing to the input matrix.")
    .def("smoothingMatrix", &Class::smoothingMatrix, "get the amoothing matrix.")
    .def("parameterCorrection", &Class::parameterCorrection, "Apply parameter correction steps.",
          py::arg("accuracy") = 1e-8,
          py::arg("maxIter") = 10,
          py::arg("tolOrth") = 1e-6)

    ;
}
#endif


} // namespace gismo
