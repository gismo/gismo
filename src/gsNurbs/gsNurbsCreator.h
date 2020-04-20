/** @file gsNurbsCreator.h

    @brief Provides declaration of the NurbsCreator struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


// Note:
// will provide functions to make bsplines, nurbs etc
// like the identity, ruled, l spapes, rings, etc
// linear extrude .....
// sweep ( eg. half-cylinder by sweeping half-circle
// revolve operation

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/**
   @brief Class gsNurbsCreator provides some simple examples of Nurbs Geometries
   
   \ingroup Nurbs
*/

template<class T>
struct gsNurbsCreator
{
    typedef memory::unique_ptr<gsGeometry<T> >        GeometryPtr;
    typedef typename gsBSpline<T>::uPtr         BSplinePtr;
    typedef typename gsNurbs<T>::uPtr           NurbsPtr;
    typedef typename gsTensorBSpline<2,T>::uPtr TensorBSpline2Ptr;
    typedef typename gsTensorBSpline<3,T>::uPtr TensorBSpline3Ptr;
    typedef typename gsTensorBSpline<4,T>::uPtr TensorBSpline4Ptr;
    typedef typename gsTensorNurbs<2,T>::uPtr   TensorNurbs2Ptr;
    typedef typename gsTensorNurbs<3,T>::uPtr   TensorNurbs3Ptr;
    typedef typename gsTensorNurbs<4,T>::uPtr   TensorNurbs4Ptr;

public:
    
    static TensorBSpline3Ptr lift3D( gsTensorBSpline<2,T> const & geo, T z = 1);
    
    static TensorBSpline4Ptr lift4D( gsTensorBSpline<3,T> const & geo, T z = 1);

    static TensorNurbs3Ptr lift3D( gsTensorNurbs<2,T> const & geo, T z = 1);

    static TensorNurbs4Ptr lift4D( gsTensorNurbs<3,T> const & geo, T z = 1);

    /* Computes a set of control points, weights, and knots that define an order-3 circular arc centered at the origin
    \param X Defines the X axis of the plane containing the arc
    \param Y Defines the Y axis of the plane containing the arc
    \param StartAngle Start angle of the arc in radians
    \param EndAngle End angle of the arc in radians
    \param Segments The number of NURBS segments in the resulting arc
    \param Knots Output container for the resulting arc knot vector
    \param Weights Output container for the resulting arc control point weights
    \param ControlPoints Output container for the resulting arc control point positions

    static gsTensorNurbs<3,T> * circularArc(const vector3& X, const vector3& Y,
    const T StartAngle, const T EndAngle,
    const unsinged Segments = 1);
*/

    static BSplinePtr BSplineUnitInterval(short_t deg);

/// 2d-rectange [low_x,upp_x] x [low_y,upp_y], rotated by \a turndeg degrees.
    static TensorBSpline2Ptr BSplineRectangle( T const & low_x = 0,
                                                    T const & low_y = 0,
                                                    T const & upp_x = 1,
                                                    T const & upp_y = 1, T const & turndeg = 0);

    // Rectangle described by the identity mapping over the given parameter domain, using tensor product B-splines.
    static TensorBSpline2Ptr BSplineRectangleWithPara( T low_x = 0, T low_y = 0, T upp_x = 1, T upp_y = 1);

/// Creates a \em n X \em m rectangle multipatch consisting of B-splines squares
/// with lower left corner at at (lx,ly).
/// The numbering of the patches are (for \em n = 4 and \em m = 3):
/// --|---|---|--
/// --|---|---|--
/// 2 | 5 | 8 |11
/// 1 | 4 | 7 |10
/// 0 | 3 | 6 |9
/// \param n number of squares in x-direction.
/// \param m number of squares in y-direction.
/// \param r with length of the side of the squares.
/// \param lx x-coordinate for lower left corner of the rectangle.
/// \param ly y-coordinate for lower left corner of the rectangle.
    static gsMultiPatch<T> BSplineSquareGrid(int n, int m, T const & r = 1,
                                             T const & lx = 0, T const & ly = 0);

    static TensorBSpline2Ptr BSplineSquare( gsMatrix<T> const & Box);

    // Note: this can probably be removed once we have degree elevation for tensor B-splines.
    //
    /// The unit square represented as a tensor B-spline of degree \a deg
    static TensorBSpline2Ptr BSplineSquareDeg(short_t deg, T scale = 1.0);

    /// Square of side \a r, with lower left corner at (x,y)
    static TensorBSpline2Ptr BSplineSquare( T const & r = 1, T const & x = 0, T const & y = 0  );

    static TensorBSpline3Ptr BSplineCube( T const & r = 1, T const & x = 0,
                                               T const & y = 0, T const & z = 0  );

    // Note: this can probably be removed once we have degree elevation for tensor B-splines.
    //
    /// The unit cube represented as a tensor B-spline of degree \a deg
    static TensorBSpline3Ptr BSplineCube(short_t deg);

    static gsMultiPatch<T> BSplineCubeGrid(int n, int m,int p, T const & r = 1,
                                           T const & lx = 0, T const & ly = 0, T const & lz = 0);

    static TensorBSpline3Ptr BSplineHalfCube( T const & r = 1, T const & x = 0,
                                                   T const & y = 0, T const & z = 0  );
    
    static TensorNurbs3Ptr NurbsCube( T const & r = 1, T const & x = 0,
                                           T const & y = 0, T const & z = 0 );

    static TensorNurbs2Ptr NurbsQuarterAnnulus( T const & r0 =1, T const & r1 =2);
    static TensorNurbs3Ptr BSplineSaddle();
    /// Inexact annulus using B-splines
    static GeometryPtr BSplineQuarterAnnulus(const short_t & deg = 2);

    //static TensorNurbs2Ptr NurbsQuarterAnnulusMixedWithLShape();
    //static GeometryPtr BSplineQuarterAnnulusMixedWithLShape(const short_t & deg = 2);

    /// Fat annulus using B-splines, discarding the weights of the exact NURBS
    ///  Analytical formulation (when r0 = 1 and r1 = 2):
    /// (x, y) = (1 + s - s*t*t - t*t, 2*s*t -s*t*t + 2*t - t*t)
    static TensorBSpline2Ptr BSplineFatQuarterAnnulus( T const & r0 =1, T const & r1 =2);

    static TensorNurbs2Ptr NurbsSphere( T const & r =1, T const & x = 0, T const & y = 0, T const & z = 0);

    static NurbsPtr  NurbsCircle( T const & r =T(1), T const & x = 0, T const & y = 0);

    static BSplinePtr BSplineFatCircle( T const & r =T(1), T const & x = 0, T const & y = 0);

    static TensorBSpline2Ptr BSplineFatDisk (T const & r=1, T const & x=0, T const & y = 0);

    static NurbsPtr NurbsCurve1 (T const & r=1, T const & x=0, T const & y = 0);

    static NurbsPtr NurbsCurve2 (T const & r=1, T const & x=0, T const & y = 0);

    static NurbsPtr NurbsBean(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineE (T const & r=1, T const & x=0, T const & y = 0);

    static NurbsPtr NurbsAmoebaFull(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineLineSegment(gsMatrix<T> const & p0, gsMatrix<T> const & p1 );

    static BSplinePtr BSplineSegment(T const u0 = 0, T const u1 = 1);

    /// L-Shaped domain represented as a tensor B-spline of degree 1
    static TensorBSpline2Ptr BSplineLShape_p1(T r = 1.0);

    /// L-Shaped domain represented as a tensor B-spline of degree 2
    /// with C0-continuity across the diagonal.
    static TensorBSpline2Ptr BSplineLShape_p2C0();

    /// L-Shaped domain represented as a tensor B-spline of degree 2
    /// with C1-continuity and double control points at the corners.
    static TensorBSpline2Ptr BSplineLShape_p2C1();

    /// L-Shaped domain represented as a multipatch (3 patches) tensor B-spline of degree 2
    /// 1. Patch is the middel part, 2. Patch is the upper part, 3 Patch is the right part.
    static gsMultiPatch<T> BSplineLShapeMultiPatch_p2();

    static BSplinePtr BSplineAmoeba(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineAmoebaBig(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineAustria(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineFish(T const & r=1, T const & x=0, T const & y = 0);

    static BSplinePtr BSplineAmoeba3degree(T const & r=1, T const & x=0, T const & y = 0);

    static TensorNurbs2Ptr NurbsDisk(T const & r=1, T const & x=0, T const & y = 0);

    static TensorBSpline2Ptr NurbsQrtPlateWHoleC0();
}; // struct

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNurbsCreator.hpp)
#endif
