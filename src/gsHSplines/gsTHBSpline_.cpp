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

void pybind11_init_gsHTensorBasis2(py::module &m)
{
  using Base  = gsBasis<real_t>;
  using Class = gsHTensorBasis<2,real_t>;
  using Mat = gsSparseMatrix<real_t, RowMajor>;
  py::class_<Class,Base>(m, "gsHTensorBasis2")

    // Accessor methods
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("refine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")

    // Level-based refinement methods
    .def("refineToLevel", &Class::refineToLevel, py::arg("minLevel"),
         "Refine the basis uniformly up to minLevel")
    .def("refineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.refineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Refine uniformly to minLevel and return transfer matrix")
    .def("refineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.refineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Refine uniformly to minLevel and update coefficients")

    .def("refineCoarsestLevel", &Class::refineCoarsestLevel,
         "Refine only the coarsest level by one")
    .def("refineCoarsestLevel_withTransfer",
         [](Class &self, Mat &T) { self.refineCoarsestLevel_withTransfer(T); },
         py::arg("transfer"), "Refine coarsest level with transfer matrix")
    .def("refineCoarsestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.refineCoarsestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Refine coarsest level and update coefficients")

    // Level-based unrefinement methods
    .def("unrefineToLevel", &Class::unrefineToLevel, py::arg("minLevel"),
         "Unrefine the basis uniformly down to minLevel")
    .def("unrefineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.unrefineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Unrefine uniformly to minLevel and return transfer matrix")
    .def("unrefineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.unrefineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Unrefine uniformly to minLevel and update coefficients")

    .def("unrefineFinestLevel", &Class::unrefineFinestLevel,
         "Unrefine only the finest level by one")
    .def("unrefineFinestLevel_withTransfer",
         [](Class &self, Mat &T) { self.unrefineFinestLevel_withTransfer(T); },
         py::arg("transfer"), "Unrefine finest level with transfer matrix")
    .def("unrefineFinestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.unrefineFinestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Unrefine finest level and update coefficients")

    // withCoefs and withTransfer variants
    .def("uniformRefine_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk, int m, short_t d) {
           self.uniformRefine_withCoefs(coefs, nk, m, d);
         },
         py::arg("coefs"), py::arg("numKnots")=1, py::arg("mul")=1, py::arg("dir")=-1,
         "Uniformly refine and update coefficients")
    .def("uniformCoarsen_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk) {
           self.uniformCoarsen_withCoefs(coefs, nk);
         },
         py::arg("coefs"), py::arg("numKnots")=1,
         "Uniformly coarsen and update coefficients")

    .def("refineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements and update coefficients")
    .def("refineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements and return transfer matrix")
    .def("refineElements_withTransfer2",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer2(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements (variant 2) and return transfer matrix")
    .def("refineElements_withCoefs2",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs2(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements (variant 2) and update coefficients")

    .def("unrefineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.unrefineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Unrefine elements and update coefficients")
    .def("unrefineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.unrefineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Unrefine elements and return transfer matrix")
    ;
}

void pybind11_init_gsHTensorBasis3(py::module &m)
{
  using Base  = gsBasis<real_t>;
  using Class = gsHTensorBasis<3,real_t>;
  using Mat = gsSparseMatrix<real_t, RowMajor>;
  py::class_<Class,Base>(m, "gsHTensorBasis3")

    // Accessor methods
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("refine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")

    // Level-based refinement methods
    .def("refineToLevel", &Class::refineToLevel, py::arg("minLevel"),
         "Refine the basis uniformly up to minLevel")
    .def("refineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.refineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Refine uniformly to minLevel and return transfer matrix")
    .def("refineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.refineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Refine uniformly to minLevel and update coefficients")

    .def("refineCoarsestLevel", &Class::refineCoarsestLevel,
         "Refine only the coarsest level by one")
    .def("refineCoarsestLevel_withTransfer",
         [](Class &self, Mat &T) { self.refineCoarsestLevel_withTransfer(T); },
         py::arg("transfer"), "Refine coarsest level with transfer matrix")
    .def("refineCoarsestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.refineCoarsestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Refine coarsest level and update coefficients")

    // Level-based unrefinement methods
    .def("unrefineToLevel", &Class::unrefineToLevel, py::arg("minLevel"),
         "Unrefine the basis uniformly down to minLevel")
    .def("unrefineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.unrefineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Unrefine uniformly to minLevel and return transfer matrix")
    .def("unrefineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.unrefineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Unrefine uniformly to minLevel and update coefficients")

    .def("unrefineFinestLevel", &Class::unrefineFinestLevel,
         "Unrefine only the finest level by one")
    .def("unrefineFinestLevel_withTransfer",
         [](Class &self, Mat &T) { self.unrefineFinestLevel_withTransfer(T); },
         py::arg("transfer"), "Unrefine finest level with transfer matrix")
    .def("unrefineFinestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.unrefineFinestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Unrefine finest level and update coefficients")

    // withCoefs and withTransfer variants
    .def("uniformRefine_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk, int m, short_t d) {
           self.uniformRefine_withCoefs(coefs, nk, m, d);
         },
         py::arg("coefs"), py::arg("numKnots")=1, py::arg("mul")=1, py::arg("dir")=-1,
         "Uniformly refine and update coefficients")
    .def("uniformCoarsen_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk) {
           self.uniformCoarsen_withCoefs(coefs, nk);
         },
         py::arg("coefs"), py::arg("numKnots")=1,
         "Uniformly coarsen and update coefficients")

    .def("refineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements and update coefficients")
    .def("refineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements and return transfer matrix")
    .def("refineElements_withTransfer2",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer2(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements (variant 2) and return transfer matrix")
    .def("refineElements_withCoefs2",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs2(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements (variant 2) and update coefficients")

    .def("unrefineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.unrefineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Unrefine elements and update coefficients")
    .def("unrefineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.unrefineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Unrefine elements and return transfer matrix")
    ;
}

void pybind11_init_gsHTensorBasis4(py::module &m)
{
  using Base  = gsBasis<real_t>;
  using Class = gsHTensorBasis<4,real_t>;
  using Mat = gsSparseMatrix<real_t, RowMajor>;
  py::class_<Class,Base>(m, "gsHTensorBasis4")

    // Accessor methods
    .def("tensorLevel",&Class::tensorLevel,"Returns the tensor basis on level i")
    .def("refine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::refine  ), "Refines the basis given a box")

    // Level-based refinement methods
    .def("refineToLevel", &Class::refineToLevel, py::arg("minLevel"),
         "Refine the basis uniformly up to minLevel")
    .def("refineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.refineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Refine uniformly to minLevel and return transfer matrix")
    .def("refineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.refineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Refine uniformly to minLevel and update coefficients")

    .def("refineCoarsestLevel", &Class::refineCoarsestLevel,
         "Refine only the coarsest level by one")
    .def("refineCoarsestLevel_withTransfer",
         [](Class &self, Mat &T) { self.refineCoarsestLevel_withTransfer(T); },
         py::arg("transfer"), "Refine coarsest level with transfer matrix")
    .def("refineCoarsestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.refineCoarsestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Refine coarsest level and update coefficients")

    // Level-based unrefinement methods
    .def("unrefineToLevel", &Class::unrefineToLevel, py::arg("minLevel"),
         "Unrefine the basis uniformly down to minLevel")
    .def("unrefineToLevel_withTransfer",
         [](Class &self, index_t lvl, Mat &T) { self.unrefineToLevel_withTransfer(lvl, T); },
         py::arg("minLevel"), py::arg("transfer"),
         "Unrefine uniformly to minLevel and return transfer matrix")
    .def("unrefineToLevel_withCoefs",
         [](Class &self, index_t lvl, gsMatrix<real_t> &coefs) { self.unrefineToLevel_withCoefs(lvl, coefs); },
         py::arg("minLevel"), py::arg("coefs"),
         "Unrefine uniformly to minLevel and update coefficients")

    .def("unrefineFinestLevel", &Class::unrefineFinestLevel,
         "Unrefine only the finest level by one")
    .def("unrefineFinestLevel_withTransfer",
         [](Class &self, Mat &T) { self.unrefineFinestLevel_withTransfer(T); },
         py::arg("transfer"), "Unrefine finest level with transfer matrix")
    .def("unrefineFinestLevel_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs) { self.unrefineFinestLevel_withCoefs(coefs); },
         py::arg("coefs"), "Unrefine finest level and update coefficients")

    // withCoefs and withTransfer variants
    .def("uniformRefine_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk, int m, short_t d) {
           self.uniformRefine_withCoefs(coefs, nk, m, d);
         },
         py::arg("coefs"), py::arg("numKnots")=1, py::arg("mul")=1, py::arg("dir")=-1,
         "Uniformly refine and update coefficients")
    .def("uniformCoarsen_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, int nk) {
           self.uniformCoarsen_withCoefs(coefs, nk);
         },
         py::arg("coefs"), py::arg("numKnots")=1,
         "Uniformly coarsen and update coefficients")

    .def("refineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements and update coefficients")
    .def("refineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements and return transfer matrix")
    .def("refineElements_withTransfer2",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.refineElements_withTransfer2(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Refine elements (variant 2) and return transfer matrix")
    .def("refineElements_withCoefs2",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.refineElements_withCoefs2(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Refine elements (variant 2) and update coefficients")

    .def("unrefineElements_withCoefs",
         [](Class &self, gsMatrix<real_t> &coefs, const std::vector<index_t> &boxes) {
           self.unrefineElements_withCoefs(coefs, boxes);
         },
         py::arg("coefs"), py::arg("boxes"),
         "Unrefine elements and update coefficients")
    .def("unrefineElements_withTransfer",
         [](Class &self, const std::vector<index_t> &boxes, Mat &T) {
           self.unrefineElements_withTransfer(boxes, T);
         },
         py::arg("boxes"), py::arg("transfer"),
         "Unrefine elements and return transfer matrix")
    ;
}

void pybind11_init_gsHBSplineBasis2(py::module &m)
{
  using Base  = gsHTensorBasis<2,real_t>;
  using Class = gsHBSplineBasis<2,real_t>;
  py::class_<Class,Base>(m, "gsHBSplineBasis2")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<2,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())
    ;
}

void pybind11_init_gsHBSplineBasis3(py::module &m)
{
	using Base  = gsHTensorBasis<3,real_t>;
	using Class = gsHBSplineBasis<3,real_t>;
	py::class_<Class,Base>(m, "gsHBSplineBasis3")


    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<3,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())

    ;
}

void pybind11_init_gsHBSplineBasis4(py::module &m)
{
	using Base  = gsHTensorBasis<4,real_t>;
	using Class = gsHBSplineBasis<4,real_t>;
	py::class_<Class,Base>(m, "gsHBSplineBasis4")

    // Constructors
    .def(py::init<>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, std::vector<index_t>   &>())
    .def(py::init<gsTensorBSplineBasis<4,real_t> const&, gsMatrix<real_t> const &>())
    .def(py::init<gsBasis<real_t> const&                                         >())

    ;
}

void pybind11_init_gsHBSpline2(py::module &m)
{
	using Base  = gsGeometry<real_t>;
	using Class = gsHBSpline<2,real_t>;
	py::class_<Class,Base>(m, "gsHBSpline2")


	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<2,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<2,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<2,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1, py::arg("dir") = -1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

void pybind11_init_gsHBSpline3(py::module &m)
{
	using Base  = gsGeometry<real_t>;
	using Class = gsHBSpline<3,real_t>;
	py::class_<Class,Base>(m, "gsHBSpline3")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<3,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<3,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<3,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1, py::arg("dir") = -1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

void pybind11_init_gsHBSpline4(py::module &m)
{
	using Base  = gsGeometry<real_t>;
	using Class = gsHBSpline<4,real_t>;
	py::class_<Class,Base>(m, "gsHBSpline4")

	// Constructors
	.def(py::init<>())
	// this one does not work:
	// .def(py::init<const gsHBSplineBasis<4,real_t> *, const gsMatrix<real_t> * >())
	.def(py::init<const gsHBSplineBasis<4,real_t> &, const gsMatrix<real_t> & >())
	.def(py::init<const gsTensorBSpline<4,real_t> &                           >())

	// Member functions
	.def("domainDim", 			&Class::domainDim,
			"Returns the domain dimension"					)
	.def("eval_into", 			&Class::eval_into,
			"Evaluates the values into a matrix"			)
	.def("deriv_into", 			&Class::deriv_into,
			"Evaluates the derivatives into a matrix"		)
	.def("deriv2_into", 		&Class::deriv2_into,
			"Evaluates the second derivatives into a matrix")

	// Derived from gsHTensorBasis
	.def("size",&Class::size,"Returns the domain dimension")
	.def("uniformRefine", static_cast<void (Class::*)(int, int, int)> (&Class::uniformRefine), "Refines the basis uniformly",
		py::arg("numKnots") = 1, py::arg("mul") = 1, py::arg("dir") = -1) //default arguments

	// .def("unrefine", static_cast<void (Class::*)(gsMatrix<real_t> const &     )> (&Class::unrefine), "Refines the basis given a box")
	.def("refineElements  ",&Class::refineElements  ,"Refines the basis given elements  ")
	.def("unrefineElements",&Class::unrefineElements,"Unrefines the basis given elements")
	;
}

#endif

}
