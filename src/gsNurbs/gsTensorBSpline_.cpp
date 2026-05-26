
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSpline.hpp>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{
TEMPLATE_INST
void constructCoefsForSlice<1, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 1> & sizes,
                                       gsMatrix<real_t> & result
                                      );

TEMPLATE_INST
void constructCoefsForSlice<2, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 2> & sizes,
                                       gsMatrix<real_t> & result
                                      );

TEMPLATE_INST
void constructCoefsForSlice<3, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 3> & sizes,
                                       gsMatrix<real_t> & result
                                      );
TEMPLATE_INST
void constructCoefsForSlice<4, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 4> & sizes,
                                       gsMatrix<real_t> & result
                                      );


CLASS_TEMPLATE_INST gsTensorBSpline<1,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<2,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<3,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<4,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<5,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<6,real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<1,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<2,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<3,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<4,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<5,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<6,real_t> >;


#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

short_t tensor_spline_dimension_from_args(const py::args & args)
{
    if (args.size() == 0)
        throw py::value_error("Cannot infer tensor spline dimension without constructor arguments.");

    const py::handle first = args[0];
    if (!py::hasattr(first, "dim"))
        throw py::value_error("Cannot infer tensor spline dimension: first argument has no dim() method.");

    const short_t d = py::cast<short_t>(first.attr("dim")());
    if (d < 2 || d > 6)
        throw py::value_error("Expected inferred dimension in [2, 6] for gsTensorBSpline factory.");

    return d;
}

void pybind11_init_gsTensorBSpline_factory(py::module &m)
{
    m.def("gsTensorBSpline", [module = py::module_(m)](py::args args) -> py::object
    {
        const short_t d = tensor_spline_dimension_from_args(args);
        const std::string className = "gsTensorBSpline" + std::to_string(d);
        py::object ctor = module.attr(className.c_str());

        PyObject * result = PyObject_CallObject(ctor.ptr(), args.ptr());
        if (!result)
            throw py::error_already_set();

        return py::reinterpret_steal<py::object>(result);
    },
    "Factory constructor that dispatches to gsTensorBSpline2..6 based on inferred dimension");
}

template <short_t d>
void pybind11_init_gsTensorBSpline(py::module &m)
{
    using Base  = gsGeometry<real_t>;
    using Class = gsTensorBSpline<d,real_t>;
    py::class_<Class,Base>(m, ("gsTensorBSpline" + std::to_string(d)).c_str())

    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<d,real_t>, gsMatrix<real_t>>())

    .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
    .def("insertKnot", &Class::insertKnot, "Insert a knot in the knot vector")
    .def("degree", &Class::degree, "Returns the degree")
    .def("degreeElevate", &Class::degreeElevate, "Elevate the degree of the spline")
    ;
}

template void pybind11_init_gsTensorBSpline<2>(py::module &m);
template void pybind11_init_gsTensorBSpline<3>(py::module &m);
template void pybind11_init_gsTensorBSpline<4>(py::module &m);
template void pybind11_init_gsTensorBSpline<5>(py::module &m);
template void pybind11_init_gsTensorBSpline<6>(py::module &m);

#endif

}
