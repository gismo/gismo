
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
        .def(py::init<>())
        .def(py::init<const gsFunctionSet<real_t> &, typename gsFunctionSet<real_t>::Ptr, const bool>(), py::arg("mp"), py::arg("fs"), py::arg("isparam"))
        .def(py::init<const gsGeometry<real_t> &, const gsFunctionSet<real_t> &, const bool>(), py::arg("sp"), py::arg("pf"), py::arg("isparam") = false)
        .def(py::init<const gsMultiPatch<real_t> &, const gsFunctionSet<real_t> &, const bool>(), py::arg("mp"), py::arg("f"), py::arg("isparam") = false)

        // Member functions
        .def("point", &Class::point, "Maps a point from the parametric to the physical domain", py::arg("u"), py::arg("i") = 0)
        .def("value", &Class::value, "Evaluates the field at given points", py::arg("u"), py::arg("i") = 0)
        .def("pvalue", &Class::pvalue, "Evaluates the field at given physical points", py::arg("u"), py::arg("i"))
        .def("distanceL2", py::overload_cast<const gsField<real_t> &>(&Class::distanceL2, py::const_), "Computes the L2-distance between two fields", py::arg("field"))
        .def("distanceL2", py::overload_cast<const gsFunctionSet<real_t> &, bool>(&Class::distanceL2, py::const_), "Computes the L2-distance between a field and a function", py::arg("func"), py::arg("isFunc_param")=false)
        .def("distanceL2", py::overload_cast<const gsFunctionSet<real_t> &, const gsMultiBasis<real_t> &, bool>(&Class::distanceL2, py::const_), "Computes the L2-distance between a field and a function using a given mesh", py::arg("func"), py::arg("mesh"), py::arg("isFunc_param")=false)
        .def("distanceH1", py::overload_cast<const gsFunctionSet<real_t> &, bool>(&Class::distanceH1, py::const_), "Computes the H1-seminorm of the difference between a field and a function", py::arg("func"), py::arg("isFunc_param")=false)
        .def("distanceH1", py::overload_cast<const gsFunctionSet<real_t> &, const gsMultiBasis<real_t> &, bool>(&Class::distanceH1, py::const_), "Computes the H1-seminorm of the difference between a field and a function using a given mesh", py::arg("func"), py::arg("mesh"), py::arg("isFunc_param")=false)
        .def("distanceH2", &Class::distanceH2, "Computes the H2-seminorm of the difference between a field and a function", py::arg("func"), py::arg("isFunc_param")=false)
        .def("distanceDG", &Class::distanceDG, "Computes the DG-distance between a field and a function", py::arg("func"), py::arg("isFunc_param")=false)
        .def("__str__",
            [](const Class& self) {
                std::ostringstream os;
                self.print(os);
                return os.str();
            },
            "Returns the string representation of the field")
        .def("parDim", &Class::parDim, "Returns the parametric dimension of the field")
        .def("geoDim", &Class::geoDim, "Returns the geometric dimension of the field")
        .def("dim", &Class::dim, "Returns the dimension of the field")

        .def("geometry", &Class::geometry, "Returns the geometry of the field", py::return_value_policy::reference_internal)
        .def("patches", &Class::patches, "Returns the multipatch domain of the field", py::return_value_policy::reference_internal)
        .def("fields", &Class::fields, "Returns the fields defined on the patches", py::return_value_policy::reference_internal)
        .def("patch", &Class::patch, "Returns the geometry of the i-th patch", py::arg("i")=0, py::return_value_policy::reference_internal)
        .def("function", &Class::function, "Returns the field defined on the i-th patch", py::arg("i")=0, py::return_value_policy::reference_internal)
        ;
    }

#endif

}
