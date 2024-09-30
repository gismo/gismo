#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunction.h>
#include <gsCore/gsFunction.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsFunction<real_t> ;

#ifdef GISMO_WITH_PYBIND11  

namespace py = pybind11;

void pybind11_init_gsFunction(py::module &m)
{
    using Base = gsFunctionSet<real_t>;
    using Class = gsFunction<real_t>;
    py::class_<Class, Base>(m, "gsFunction")
        .def("jacobian",  &Class::jacobian, "Returns the Jacobian")
        .def("hessian",   &Class::hessian, "Returns the Hessian")
        .def("laplacian", &Class::laplacian, "Returns the Laplacian")
        .def("argMin", &Class::argMin, "Returns the position of the minimum",
             py::arg("accuracy") = 1e-6, py::arg("max_loop") = 100,
             py::arg("init") = gsMatrix<real_t>(),
             py::arg("damping_factor") = 1 )
  ;
}
#endif

}
