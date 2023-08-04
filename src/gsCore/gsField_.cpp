
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsField.h>
#include <gsCore/gsField.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsField<real_t>;

#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsField(py::module &m)
    {
        using Class = gsField<real_t>;
        py::class_<Class>(m, "gsField")
                // Constructors
                //Member functions
                .def("distanceL2", static_cast< real_t (Class::*)(gsFunctionSet<real_t> const &, bool, int) const> (&Class::distanceL2), "Compute L2 distance")
                .def("distanceH1", static_cast< real_t (Class::*)(gsFunctionSet<real_t> const &, bool, int) const> (&Class::distanceH1), "Compute H1 semi-norm distance")
                ;
    }


#endif
}
