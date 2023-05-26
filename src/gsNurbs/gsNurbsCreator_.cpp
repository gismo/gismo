#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsNurbsCreator.h>
#include <gsNurbs/gsNurbsCreator.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsNurbsCreator<real_t>;

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

void pybind11_init_gsNurbsCreator(py::module &m)
{
  	using Class = gsNurbsCreator<real_t>;
	py::class_<Class>(m, "gsNurbsCreator")
        .def(py::init<>())
	    .def_static("lift3D", py::detail::overload_cast_impl<gsTensorBSpline<2,real_t> const &, real_t>()(&Class::lift3D),
	    			"Lifts a 2D geometry into 3D space")
	    .def_static("lift4D", py::detail::overload_cast_impl<gsTensorBSpline<3,real_t> const &, real_t>()(&Class::lift4D),
	    			"Lifts a 3D geometry into 4D space")
	    .def_static("lift3D", py::detail::overload_cast_impl<gsTensorNurbs<2,real_t> const &, real_t>()(&Class::lift3D),
	    			"Lifts a 2D geometry into 3D space")
	    .def_static("lift4D", py::detail::overload_cast_impl<gsTensorNurbs<3,real_t> const &, real_t>()(&Class::lift4D),
	    			"Lifts a 3D geometry into 4D space")

        .def_static("BSplineUnitInterval", &Class::BSplineUnitInterval,
        			"Creates a B-Spline unit interval")
	    .def_static("BSplineRectangle", &Class::BSplineRectangle,
	    			"2d-rectange [low_x,upp_x] x [low_y,upp_y], rotated by turndeg degrees.",
	    			py::arg("low_x")=0, py::arg("low_y")=0, py::arg("upp_x")=1, py::arg("upp_y")=1, py::arg("turn_def")=0)
	    .def_static("BSplineRectangleWithPara", &Class::BSplineRectangleWithPara,
					"Rectangle described by the identity mapping over the given parameter domain, using tensor product B-splines.",
	    			py::arg("low_x")=0, py::arg("low_y")=0, py::arg("upp_x")=1, py::arg("upp_y")=1)
	    .def_static("BSplineSquareGrid", &Class::BSplineSquareGrid,
					"Creates a \em n X \em m rectangle multipatch consisting of B-splines squares with lower left corner at at (lx,ly)",
	    			py::arg("n"),py::arg("m"),py::arg("r")=1, py::arg("lx")=0, py::arg("ly")=0)
	    .def_static("BSplineSquareDeg", &Class::BSplineSquareDeg,
					"The unit square represented as a tensor B-spline of degree deg",
	    			py::arg("deg"),py::arg("scale")=1)
	    .def_static("BSplineSquare", py::detail::overload_cast_impl<gsMatrix<real_t> const &>()(&Class::BSplineSquare),
	    			"Square of side defined by two corners (columsn)")
	    .def_static("BSplineSquare", py::detail::overload_cast_impl<real_t const &, real_t const &, real_t const &>()(&Class::BSplineSquare),
	    			"Square of side r, with lower left corner at (x,y)")

	    .def_static("BSplineCube", py::detail::overload_cast_impl<real_t const &, real_t const &, real_t const &, real_t const &>()(&Class::BSplineCube),
	    			"Cube of side \a r, with lower left corner at (x,y,z)")
	    .def_static("BSplineCube", py::detail::overload_cast_impl<short_t>()(&Class::BSplineCube),
	    			"The unit cube represented as a tensor B-spline of degree \a deg")

	    .def_static("BSplineCubeGrid", &Class::BSplineCubeGrid,
					"Creates a \em n X \em m X p rectangle multipatch consisting of B-splines squares with lower left corner at at (lx,ly,lz)",
	    			py::arg("n"),py::arg("m"),py::arg("p"),py::arg("r")=1, py::arg("lx")=0, py::arg("ly")=0, py::arg("lz")=0)
	    .def_static("BSplineHalfCube", &Class::BSplineHalfCube,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
	    .def_static("NurbsCube", &Class::NurbsCube,
					"Cube of side \a r, with lower left corner at (x,y,z) using NURBS",
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
	    .def_static("NurbsQuarterAnnulus", &Class::NurbsQuarterAnnulus,
					"Exact annulus using NURBS",
	    			py::arg("r0")=0, py::arg("r1")=1)
        .def_static("NurbsAnnulus", &Class::NurbsAnnulus,
    				"Exact full annulus using NURBS",
        			py::arg("r0")=0, py::arg("r1")=1)
	    .def_static("BSplineSaddle", &Class::BSplineSaddle,
					"Saddle using B-splines")
	    .def_static("BSplineQuarterAnnulus", &Class::BSplineQuarterAnnulus,
					"Inexact annulus using B-splines",
	    			py::arg("deg")=2)
	    .def_static("BSplineFatQuarterAnnulus", &Class::BSplineFatQuarterAnnulus,
					"Fat annulus using B-splines, discarding the weights of the exact NURBS",
	    			py::arg("r0")=0, py::arg("r1")=1)
	    .def_static("NurbsSphere", &Class::NurbsSphere,
					"Sphere using NURBS",
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
	    .def_static("NurbsCircle", &Class::NurbsCircle,
					"Circle using NURBS",
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineFatCircle", &Class::BSplineFatCircle,
					" Inexact circle using B-splines",
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineFatDisk", &Class::BSplineFatDisk,
					"Inexact disk using B-splines",
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsCurve1", &Class::NurbsCurve1,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsCurve2", &Class::NurbsCurve2,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsBean", &Class::NurbsBean,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineE", &Class::BSplineE,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsAmoebaFull", &Class::NurbsAmoebaFull,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineLineSegment", &Class::BSplineLineSegment,
					"")//?
	    .def_static("BSplineSegment", &Class::BSplineSegment,
					"",
	    			py::arg("u0")=0, py::arg("u1")=1)
	    .def_static("BSplineLShape_p1", &Class::BSplineLShape_p1,
					"L-Shaped domain represented as a tensor B-spline of degree 1",
	    			py::arg("r")=1)
	    .def_static("BSplineLShape_p2C0", &Class::BSplineLShape_p2C0,
					"L-Shaped domain represented as a tensor B-spline of degree 2, with C0-continuity across the diagonal.")
	    .def_static("BSplineLShape_p2C1", &Class::BSplineLShape_p2C1,
					"L-Shaped domain represented as a tensor B-spline of degree 2, with C1-continuity and double control points at the corners.")
	    .def_static("BSplineLShapeMultiPatch_p2", &Class::BSplineLShapeMultiPatch_p2,
					"L-Shaped domain represented as a multipatch (3 patches) tensor B-spline of degree 2. 1. Patch is the middel part, 2. Patch is the upper part, 3 Patch is the right part.")
	    .def_static("BSplineAmoeba", &Class::BSplineAmoeba,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineAmoebaBig", &Class::BSplineAmoebaBig,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineAustria", &Class::BSplineAustria,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineFish", &Class::BSplineFish,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("BSplineAmoeba3degree", &Class::BSplineAmoeba3degree,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsDisk", &Class::NurbsDisk,
					"",//?
	    			py::arg("r")=1, py::arg("x")=0, py::arg("y")=0)
	    .def_static("NurbsQrtPlateWHoleC0", &Class::NurbsQrtPlateWHoleC0,
					"")//?
	    .def_static("BSplineTriangle", &Class::BSplineTriangle,
					"Makes a Isosceles triangle with height H and width W ",
    			py::arg("H")=1, py::arg("W")=0)
        .def_static("BSplineTriangle", &Class::BSplineTriangle,
    				"Makes a star with N patches, outer radius R0 and inner radius R1",
        			py::arg("H")=1, py::arg("W")=0)
    ;
}

#endif

}
