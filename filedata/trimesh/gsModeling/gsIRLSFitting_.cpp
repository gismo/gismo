#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsIRLSFitting.h>
#include <gsModeling/gsIRLSFitting.hpp>
#include <gsNurbs/gsBSplineBasis.h>

namespace gismo
{

CLASS_TEMPLATE_INST gsIRLSFitting<real_t>;

//#ifdef GISMO_BUILD_PYBIND11
//
//namespace py = pybind11;
//void pybind11_init_gsIRLSFitting(py::module &m)
//{
//  using Class = gsIRLSFitting<real_t>;
//  py::class_<Class>(m, "gsIRLSFitting")
//
//    // Constructors
//    .def( py::init<>() ) // Empty constructor
//    .def( py::init<gsMatrix<real_t> const &, gsMatrix<real_t> const &, gsBasis<real_t>&>() )
//
//    // Member functions
//    .def("compute", &Class::compute, "Computes the least square fit for a gsBasis.")
//    .def("applySmoothing", &Class::applySmoothing, "apply smoothing to the input matrix.")
//    .def("smoothingMatrix", &Class::smoothingMatrix, "get the amoothing matrix.")
//    .def("parameterCorrection", &Class::parameterCorrection, "Apply parameter correction steps.")
//    ;
//}
//#endif


} // namespace gismo
