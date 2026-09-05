#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsConstantFunction.hpp>
#include <gsIO/gsXmlRegistry.h>


namespace gismo
{

CLASS_TEMPLATE_INST gsConstantFunction<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsConstantFunction<real_t> >;


// XML dispatch registration (priorities: see gsCoreXmlRegistration.h)
GISMO_XML_REGISTER(gsFunction<real_t>, gsConstantFunction<real_t>, 100)
#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsConstantFunction(py::module &m)
{
    using Base  = gsFunction<real_t>;
    using Class = gsConstantFunction<real_t>;
    py::class_<Class,Base>(m, "gsConstantFunction")

    // Constructors
    .def(py::init<const gsVector<real_t>&, short_t>(), py::arg("val"), py::arg("domainDim"))
    .def(py::init<real_t, short_t>(), py::arg("x"), py::arg("domainDim"))
    .def(py::init<real_t, real_t, short_t>(), py::arg("x"), py::arg("y"), py::arg("domainDim"))
    .def(py::init<real_t, real_t, real_t, short_t>(), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("domainDim"))
    .def(py::init<real_t, real_t, real_t, real_t, short_t>(), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("w"), py::arg("domainDim"))
    ;
}

#endif

}
