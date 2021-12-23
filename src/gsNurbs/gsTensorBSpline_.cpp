
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


#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsTensorBSpline2(py::module &m)
{
  using Base  = gsGeometry<real_t>;
  using Class = gsTensorBSpline<2,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSpline2")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,    gsKnotVector<real_t>,   gsMatrix<real_t> >())
    .def(py::init<gsTensorBSplineBasis<2,real_t>, gsMatrix<real_t>>())

    // Member functions
    .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
    .def("degree", &Class::degree, "Returns the degree")
    ;
}

void pybind11_init_gsTensorBSpline3(py::module &m)
{
  using Base  = gsGeometry<real_t>;
  using Class = gsTensorBSpline<3,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSpline3")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,
                  gsKnotVector<real_t>,
                  gsKnotVector<real_t>,
                  gsMatrix<real_t>
                  >())
    .def(py::init<gsTensorBSplineBasis<3,real_t>, gsMatrix<real_t>>())

    // Member functions
    .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
    .def("degree", &Class::degree, "Returns the degree")
    ;
}

void pybind11_init_gsTensorBSpline4(py::module &m)
{
  using Base  = gsGeometry<real_t>;
  using Class = gsTensorBSpline<4,real_t>;
  py::class_<Class,Base>(m, "gsTensorBSpline4")

    // Constructors
    .def(py::init<gsKnotVector<real_t>,    gsKnotVector<real_t>,    gsKnotVector<real_t>,    gsKnotVector<real_t>    ,   gsMatrix<real_t> >())
    .def(py::init<gsTensorBSplineBasis<4,real_t>, gsMatrix<real_t>>())

    // Member functions
    .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
    .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
    .def("degree", &Class::degree, "Returns the degree")
    ;
}

#endif

}
