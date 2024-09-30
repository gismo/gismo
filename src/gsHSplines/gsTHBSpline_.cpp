#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.hpp>

#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsTHBSplineBasis <1, real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <2,real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <3,real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <4,real_t>;

CLASS_TEMPLATE_INST gsTHBSpline      <1,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <2,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <3,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <4,real_t>;

namespace internal
{
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsTHBSpline<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<4,real_t> >;
}

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsTHBSplineBasis2(py::module &m)
{
  using Base  = gsHTensorBasis<2,real_t>;
  using Class = gsTHBSplineBasis<2,real_t>;
  py::class_<Class,Base>(m, "gsTHBSplineBasis2")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())
    ;
}

void pybind11_init_gsTHBSplineBasis3(py::module &m)
{
  using Base  = gsHTensorBasis<3,real_t>;
  using Class = gsTHBSplineBasis<3,real_t>;
  py::class_<Class,Base>(m, "gsTHBSplineBasis3")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())
    ;
}

void pybind11_init_gsTHBSplineBasis4(py::module &m)
{
  using Base  = gsHTensorBasis<4,real_t>;
  using Class = gsTHBSplineBasis<4,real_t>;
  py::class_<Class,Base>(m, "gsTHBSplineBasis4")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())
    ;
}

void pybind11_init_gsTHBSpline2(py::module &m)
{
  using Base  = gsGeometry<real_t>;
	using Class = gsTHBSpline<2,real_t>;
	py::class_<Class,Base>(m, "gsTHBSpline2")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsTHBSplineBasis<2,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsTHBSplineBasis<2,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<2,real_t> &                      >())
	;
}

void pybind11_init_gsTHBSpline3(py::module &m)
{
  using Base  = gsGeometry<real_t>;
	using Class = gsTHBSpline<3,real_t>;
	py::class_<Class,Base>(m, "gsTHBSpline3")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsTHBSplineBasis<3,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsTHBSplineBasis<3,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<3,real_t> &                      >())
	;
}

void pybind11_init_gsTHBSpline4(py::module &m)
{
  using Base  = gsGeometry<real_t>;
	using Class = gsTHBSpline<4,real_t>;
	py::class_<Class,Base>(m, "gsTHBSpline4")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsTHBSplineBasis<4,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsTHBSplineBasis<4,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<4,real_t> &                      >())
	;
}

#endif

}
