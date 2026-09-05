#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsFunctionExpr.hpp>
#include <gsIO/gsXmlRegistry.h>


namespace gismo
{

CLASS_TEMPLATE_INST gsFunctionExpr<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsFunctionExpr<real_t> >;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsFunctionExpr(py::module &m)
{
  using Base  = gsFunction<real_t>;
  using Class = gsFunctionExpr<real_t>;
  py::class_<Class,Base>(m, "gsFunctionExpr")

    // Constructors
    .def(py::init<std::string,short_t>())
    .def(py::init<std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,std::string,short_t>())
    .def(py::init<std::string,std::string,std::string,std::string,std::string,std::string,std::string,std::string,std::string,short_t>())

    // set_u binding (Phase 0b: update the named variable "u" in expression for each load step)
    .def("set_u", &Class::set_u, py::arg("value"), "Update the named variable 'u' in the expression")
    ;
}

#endif


// XML dispatch registration (priorities: see gsCoreXmlRegistration.h)
GISMO_XML_REGISTER(gsFunction<real_t>, gsFunctionExpr<real_t>, 110)

}
