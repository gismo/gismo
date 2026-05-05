
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

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<1,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<2,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<3,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<4,real_t> >;


#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

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

#endif

}
