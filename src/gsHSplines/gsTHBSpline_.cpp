#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXml.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.hpp>

#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsHTensorBasis <1,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <2,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <3,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <4,real_t>;
  
CLASS_TEMPLATE_INST gsTHBSplineBasis <1,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <2,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <3,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <4,real_t,true>;

CLASS_TEMPLATE_INST gsTHBSplineBasis <1,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <2,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <3,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <4,real_t,false>;

CLASS_TEMPLATE_INST gsTHBSpline      <1,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSpline      <2,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSpline      <3,real_t,true>;
CLASS_TEMPLATE_INST gsTHBSpline      <4,real_t,true>;

CLASS_TEMPLATE_INST gsTHBSpline      <1,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSpline      <2,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSpline      <3,real_t,false>;
CLASS_TEMPLATE_INST gsTHBSpline      <4,real_t,false>;

namespace internal
{

CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<4,real_t> >;
  
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<1,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<2,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<3,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<4,real_t,true> >;

CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<1,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<2,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<3,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<4,real_t,false> >;

CLASS_TEMPLATE_INST gsXml< gsTHBSpline<1,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<2,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<3,real_t,true> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<4,real_t,true> >;

CLASS_TEMPLATE_INST gsXml< gsTHBSpline<1,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<2,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<3,real_t,false> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<4,real_t,false> >;

}

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

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

template void pybind11_init_gsHTensorBasis<1>(py::module &m);
template void pybind11_init_gsHTensorBasis<2>(py::module &m);
template void pybind11_init_gsHTensorBasis<3>(py::module &m);
template void pybind11_init_gsHTensorBasis<4>(py::module &m);

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

template void pybind11_init_gsTHBSplineBasis<1, true>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<2, true>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<3, true>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<4, true>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<1, false>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<2, false>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<3, false>(py::module &m);
template void pybind11_init_gsTHBSplineBasis<4, false>(py::module &m);

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

template void pybind11_init_gsTHBSpline<1, true>(py::module &m);
template void pybind11_init_gsTHBSpline<2, true>(py::module &m);
template void pybind11_init_gsTHBSpline<3, true>(py::module &m);
template void pybind11_init_gsTHBSpline<4, true>(py::module &m);
template void pybind11_init_gsTHBSpline<1, false>(py::module &m);
template void pybind11_init_gsTHBSpline<2, false>(py::module &m);
template void pybind11_init_gsTHBSpline<3, false>(py::module &m);
template void pybind11_init_gsTHBSpline<4, false>(py::module &m);

#endif

}
