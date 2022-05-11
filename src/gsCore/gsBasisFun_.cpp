#include <gsCore/gsBasis.h>
#include <gsCore/gsBasisFun.h>

namespace gismo
{
#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsBasisFun(py::module &m)
{
    using Class = gsBasisFun<real_t>;
    py::class_<Class>(m, "gsBasisFun")
        // Member functions
    .def("eval", &Class::eval, "Evaluates points into a matrix")
        ;
}

#endif

}
