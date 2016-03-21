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
    typedef typename gsBSpline<T>::uPtr         BSplinePtr;
    typedef typename gsNurbs<T>::uPtr           NurbsPtr;
    typedef typename gsTensorBSpline<2,T>::uPtr TensorBSpline2Ptr;
    typedef typename gsTensorNurbs<2,T>::uPtr   TensorNurbs2Ptr;

public:
    
    static gsTensorBSpline<3,T> * lift3D( gsTensorBSpline<2,T> const & geo, T z = 1);
    
    static gsTensorBSpline<4,T> * lift4D( gsTensorBSpline<3,T> const & geo, T z = 1);

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

    static gsBSpline<T> * BSplineUnitInterval(int deg);

/// 2d-rectange [low_x,upp_x] x [low_y,upp_y], rotated by \a turndeg degrees.
    static gsTensorBSpline<2,T> * BSplineRectangle( T const & low_x = 0,
                                                    T const & low_y = 0,
                                                    T const & upp_x = 1,
                                                    T const & upp_y = 1, T const & turndeg = 0);

    // Rectangle described by the identity mapping over the given parameter domain, using tensor product B-splines.
    static gsTensorBSpline<2,T> * BSplineRectangleWithPara( T low_x = 0, T low_y = 0, T upp_x = 1, T upp_y = 1);

/// Square of side \a r, with lower left corner at (x,y)
    static gsTensorBSpline<2,T> * BSplineSquare( T const & r = 1, T const & x = 0, T const & y = 0  );

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
    static gsMultiPatch<T> * BSplineSquareGrid(int n, int m, T const & r = 1,
                                               T const & lx = 0, T const & ly = 0);

    static gsTensorBSpline<2,T> * BSplineSquare( gsMatrix<T> const & Box);

    // Note: this can probably be removed once we have degree elevation for tensor B-splines.
    //
    /// The unit square represented as a tensor B-spline of degree \a deg
    static gsTensorBSpline<2,T> * BSplineSquare(int deg, T scale = 1.0);
    
    static gsTensorBSpline<3,T> * BSplineCube( T const & r = 1, T const & x = 0,
                                               T const & y = 0, T const & z = 0  );

    // Note: this can probably be removed once we have degree elevation for tensor B-splines.
    //
    /// The unit cube represented as a tensor B-spline of degree \a deg
    static gsTensorBSpline<3,T> * BSplineCube(int deg);

    static gsTensorBSpline<3,T> * BSplineHalfCube( T const & r = 1, T const & x = 0,
                                                   T const & y = 0, T const & z = 0  );
    
    static gsTensorNurbs<3,T> * NurbsCube( T const & r = 1, T const & x = 0,
                                           T const & y = 0, T const & z = 0 );

    static gsTensorNurbs<2,T> * NurbsQuarterAnnulus( T const & r0 =1, T const & r1 =2);

    static gsTensorNurbs<3,T> * BSplineSaddle();

    /// Inexact annulus using B-splines
    static gsGeometry<T> * BSplineQuarterAnnulus(int const & deg = 2);

    /// Fat annulus using B-splines, discarding the weights of the exact NURBS
    ///  Analytical formulation (when r0 = 1 and r1 = 2):
    /// (x, y) = (1 + s - s*t*t - t*t, 2*s*t -s*t*t + 2*t - t*t)
    static gsTensorBSpline<2,T> * BSplineFatQuarterAnnulus( T const & r0 =1, T const & r1 =2);

    static gsTensorNurbs<2,T> * NurbsSphere( T const & r =1, T const & x = 0, T const & y = 0, T const & z = 0);

    static gsNurbs<T> * NurbsCircle( T const & r =T(1), T const & x = 0, T const & y = 0);

    static gsBSpline<T> * BSplineFatCircle( T const & r =T(1), T const & x = 0, T const & y = 0);

    static gsTensorBSpline<2,T> *BSplineFatDisk (T const & r=1, T const & x=0, T const & y = 0);

    static gsNurbs<T> *NurbsCurve1 (T const & r=1, T const & x=0, T const & y = 0);

    static gsNurbs<T> *NurbsCurve2 (T const & r=1, T const & x=0, T const & y = 0);

    static gsNurbs<T> *NurbsBean(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineE (T const & r=1, T const & x=0, T const & y = 0);

    static gsNurbs<T> *NurbsAmoebaFull(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineLineSegment(gsMatrix<T> const & p0, gsMatrix<T> const & p1 );

    static BSplinePtr BSplineSegment(T const u0 = 0, T const u1 = 1);

    /// L-Shaped domain represented as a tensor B-spline of degree 1
    static gsTensorBSpline<2,T> * BSplineLShape_p1(T r = 1.0);

    /// L-Shaped domain represented as a tensor B-spline of degree 2
    /// with C0-continuity across the diagonal.
    static gsTensorBSpline<2,T> * BSplineLShape_p2C0();

    /// L-Shaped domain represented as a tensor B-spline of degree 2
    /// with C1-continuity and double control points at the corners.
    static gsTensorBSpline<2,T> * BSplineLShape_p2C1();

    /// L-Shaped domain represented as a multipatch (3 patches) tensor B-spline of degree 2
    /// 1. Patch is the middel part, 2. Patch is the upper part, 3 Patch is the right part.
    static gsMultiPatch<T> * BSplineLShapeMultiPatch_p2();

    static gsBSpline<T> *BSplineAmoeba(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineAmoebaBig(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineAustria(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineFish(T const & r=1, T const & x=0, T const & y = 0);

    static gsBSpline<T> *BSplineAmoeba3degree(T const & r=1, T const & x=0, T const & y = 0);

    static gsTensorNurbs<2,T> *NurbsDisk(T const & r=1, T const & x=0, T const & y = 0);

    static gsTensorBSpline<2,T> * NurbsQrtPlateWHoleC0();
}; // struct

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNurbsCreator.hpp)
#endif
