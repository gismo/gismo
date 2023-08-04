#include <gsAssembler/gsAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsAssembler<real_t>;

#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_enum_gsAssemblerOptions(py::module &m)
    {
        py::enum_<dirichlet::strategy>(m, "dirichletStrategy")
                .value("elimination" , dirichlet::elimination )
                .value("penalize" , dirichlet::penalize )
                .value("nitsche", dirichlet::nitsche)
                .value("none", dirichlet::none)
                .export_values();

        py::enum_<dirichlet::values>(m, "values")
                .value("homogeneous", dirichlet::homogeneous)
                .value("interpolation", dirichlet::interpolation)
                .value("l2Projection", dirichlet::l2Projection)
                .value("user", dirichlet::user)
                .export_values();

        py::enum_<iFace::strategy>(m, "iFaceStrategy")
                .value("conforming", iFace::conforming)
                .value("glue", iFace::glue)
                .value("dg", iFace::dg)
                .value("smooth", iFace::smooth)
                .value("none", iFace::none)
                .export_values();
    }

    void pybind11_init_gsAssembler(py::module &m)
    {
        //using Base  = gsGeometry<real_t>;
        using Class = gsAssembler<real_t>;
        py::class_<Class>(m, "gsAssembler")
                // Constructors
                //Member functions
                .def("options", &Class::options, "Returns the options")
                .def("matrix", &Class::matrix, "Returns the assembled matrix")
                .def("rhs", static_cast<const gsMatrix<real_t> & (Class::*) () const> (&Class::rhs), "Returhs the assembled right-hand-side")
                .def("rhs", static_cast<gsMatrix<real_t> & (Class::*) ()> (&Class::rhs), "Returhs the assembled right-hand-side")
                .def("constructSolution", static_cast<void (Class::*)(const gsMatrix<real_t> &, gsMultiPatch<real_t> &, short_t) const> (&Class::constructSolution), "Constuct the solution function")
                ;
    }


#endif
}
