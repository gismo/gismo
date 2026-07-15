#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsIO/gsXml.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.hpp>

#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSpline.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST gsHTensorBasis<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsTHBSplineBasis<D,real_t,true>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsTHBSplineBasis<D,real_t,false>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsTHBSpline<D,real_t,true>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsTHBSpline<D,real_t,false>;
GISMO_DIM_FOREACH(INST)
#undef INST

namespace internal
{

#define INST(D) CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<D,real_t> >;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<D,real_t,true> >;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<D,real_t,false> >;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsXml< gsTHBSpline<D,real_t,true> >;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST gsXml< gsTHBSpline<D,real_t,false> >;
GISMO_DIM_FOREACH(INST)
#undef INST

}

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

short_t h_spline_dimension_from_args(const py::args & args, const char * factoryName)
{
    if (args.size() == 0)
        throw py::value_error(std::string("Cannot infer dimension for ") + factoryName + " without constructor arguments.");

    const py::handle first = args[0];

    short_t d = 0;
    if (py::hasattr(first, "dim"))
    {
        d = py::cast<short_t>(first.attr("dim")());
    }
    else if (py::hasattr(first, "domainDim"))
    {
        d = py::cast<short_t>(first.attr("domainDim")());
    }
    else
    {
        throw py::value_error(std::string("Cannot infer dimension for ") + factoryName + ": first argument has neither dim() nor domainDim().");
    }



    return d;
}

py::object call_dimensioned_constructor(py::module module, const std::string & classPrefix, const short_t d, const py::args & args)
{
    const std::string className = classPrefix + std::to_string(d);
    py::object ctor = module.attr(className.c_str());

    PyObject * result = PyObject_CallObject(ctor.ptr(), args.ptr());
    if (!result)
        throw py::error_already_set();

    return py::reinterpret_steal<py::object>(result);
}

void pybind11_init_gsTHBSplineBasis_factory(py::module &m)
{
    m.def("gsTHBSplineBasis", [module = py::module_(m)](py::args args) -> py::object
    {
        const short_t d = h_spline_dimension_from_args(args, "gsTHBSplineBasis");
        return call_dimensioned_constructor(module, "gsTHBSplineBasis", d, args);
    },
    "Factory constructor that dispatches to gsTHBSplineBasis1..N based on inferred dimension");

    m.def("gsHBSplineBasis", [module = py::module_(m)](py::args args) -> py::object
    {
        const short_t d = h_spline_dimension_from_args(args, "gsHBSplineBasis");
        return call_dimensioned_constructor(module, "gsHBSplineBasis", d, args);
    },
    "Factory constructor that dispatches to gsHBSplineBasis1..N based on inferred dimension");
}

void pybind11_init_gsTHBSpline_factory(py::module &m)
{
    m.def("gsTHBSpline", [module = py::module_(m)](py::args args) -> py::object
    {
        const short_t d = h_spline_dimension_from_args(args, "gsTHBSpline");
        return call_dimensioned_constructor(module, "gsTHBSpline", d, args);
    },
    "Factory constructor that dispatches to gsTHBSpline1..N based on inferred dimension");

    m.def("gsHBSpline", [module = py::module_(m)](py::args args) -> py::object
    {
        const short_t d = h_spline_dimension_from_args(args, "gsHBSpline");
        return call_dimensioned_constructor(module, "gsHBSpline", d, args);
    },
    "Factory constructor that dispatches to gsHBSpline1..N based on inferred dimension");
}
template <short_t d>
void pybind11_init_gsHTensorBasis(py::module &m)
{
    using Base  = gsBasis<real_t>;
    using Class = gsHTensorBasis<d,real_t>;
    py::class_<Class,Base>(m, ("gsHTensorBasis" + std::to_string(d)).c_str())

    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("refine", static_cast<void (Class::*)(gsMatrix<real_t> const &)> (&Class::refine), "Refines the basis given a box")
    ;
}

#define INST(D) template void pybind11_init_gsHTensorBasis<D>(py::module &m);
GISMO_DIM_FOREACH(INST)
#undef INST

template <short_t d, bool Trunc>
void pybind11_init_gsTHBSplineBasis(py::module &m)
{
    using Base  = gsHTensorBasis<d,real_t>;
    using Class = gsTHBSplineBasis<d,real_t,Trunc>;
    py::class_<Class,Base>(m, (std::string(Trunc ? "gsTHBSplineBasis" : "gsHBSplineBasis") + std::to_string(d)).c_str())

    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<d,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<d,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())
    ;
}

#define INST(D) template void pybind11_init_gsTHBSplineBasis<D, true>(py::module &m);
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) template void pybind11_init_gsTHBSplineBasis<D, false>(py::module &m);
GISMO_DIM_FOREACH(INST)
#undef INST

template <short_t d, bool Trunc>
void pybind11_init_gsTHBSpline(py::module &m)
{
    using Base  = gsGeometry<real_t>;
    using Class = gsTHBSpline<d,real_t,Trunc>;
    py::class_<Class,Base>(m, (std::string(Trunc ? "gsTHBSpline" : "gsHBSpline") + std::to_string(d)).c_str())

    .def(py::init<>())
    .def(py::init<const gsTHBSplineBasis<d,real_t,Trunc> &, const gsMatrix<real_t> & >())
    .def(py::init<const gsTensorBSpline<d,real_t> &>())
    ;
}

#define INST(D) template void pybind11_init_gsTHBSpline<D, true>(py::module &m);
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) template void pybind11_init_gsTHBSpline<D, false>(py::module &m);
GISMO_DIM_FOREACH(INST)
#undef INST

#endif

}
