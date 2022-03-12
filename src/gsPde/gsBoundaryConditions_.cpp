
#include <gsCore/gsTemplateTools.h>

#include <gsPde/gsBoundaryConditions.h>
#include <gsPde/gsBoundaryConditions.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBoundaryConditions<real_t>;
    CLASS_TEMPLATE_INST internal::gsXml< gsBoundaryConditions<real_t> >;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_enum_gsBoundaryConditions(py::module &m)
{
    py::enum_<condition_type::type>(m, "bctype")
        .value("unknownType"   , condition_type::unknownType )
        .value("dirichlet"     , condition_type::dirichlet )
        .value("weak_dirichlet", condition_type::weak_dirichlet)
        .value("neumann"       , condition_type::neumann)
        .value("robin"         , condition_type::robin)
        .value("clamped"       , condition_type::clamped )
        .value("weak_clamped"  , condition_type::weak_clamped)
        .value("collapsed"     , condition_type::collapsed)
        .value("laplace"       , condition_type::laplace)
        .export_values();
}

void pybind11_init_gsBoundaryConditions(py::module &m)
{
  using Class = gsBoundaryConditions<real_t>;
  py::class_<Class>(m, "gsBoundaryConditions")

    // Constructors
    .def(py::init<>())
    // Member functions
    .def("clear", &Class::clear, "Clears the gsBoundaryConditions object")
    .def("size", &Class::size, "Number of boundary conditions assigned")

    .def("add", static_cast<void (Class::*)(int, boxSide, const std::string &, const gsFunction<real_t>::Ptr &, short_t, int, bool)> (&Class::add), "Adds a boundary condition")
    .def("add", static_cast<void (Class::*)(int, boxSide, const std::string &,       gsFunction<real_t>      *, short_t, int, bool)> (&Class::add), "Adds a boundary condition")
    .def("add", static_cast<void (Class::*)(int, boxSide, const std::string &, const gsFunction<real_t>      &, short_t, int, bool)> (&Class::add), "Adds a boundary condition")


    .def("addCondition", static_cast<void (Class::*)(int, boundary::side, condition_type::type, const gsFunction<real_t> &     , short_t, bool, int)> (&Class::addCondition),
                            "Adds a boundary condition"//,
                            // py::arg("unknown") = 0,
                            // py::arg("parametric") = false,
                            // py::arg("comp") = -1
                            )
    .def("addCondition", static_cast<void (Class::*)(int, boundary::side, condition_type::type,       gsFunction<real_t> *     , short_t, bool, int)> (&Class::addCondition),
                            "Adds a boundary condition"//,
                            // py::arg("unknown") = 0,
                            // py::arg("parametric") = false,
                            // py::arg("comp") = -1
                            )
    .def("addCondition", static_cast<void (Class::*)(const patchSide&, condition_type::type,       gsFunction<real_t> *     , short_t, bool, int)> (&Class::addCondition),
        "Adds a boundary condition"//,
        // py::arg("unknown") = 0,
        // py::arg("parametric") = false,
        // py::arg("comp") = -1
        )
    .def("addCornerValue", static_cast<void (Class::*)(boundary::corner, real_t, int, short_t, short_t)> (&Class::addCornerValue),
                            "Adds a boundary condition"//,
                            // py::arg("p") = 0,
                            // py::arg("unknown") = 0
                            )

    .def("setGeoMap", &Class::setGeoMap, "Sets the geometry map for the boundary computations")
    ;
}

#endif

}

