/** @file testEmbedding.hpp

    @brief Provides implementation of embed curve on surface
 Returns the spatial positions of the points on 3D space belonging to the curve

 Input:
 p_curve: the polynomial order of the curve
 Xi_curve: the knot vector of the curve in its parameter space
 */


#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbs.h>
#include "testEmbedding.h"


namespace gismo{
template<typename T>
T testEmbedding<T>::EmbedCurvesOnSurface(   const gsGeometry<T> & curve,
                                            const gsGeometry<T> & surface,
                                            gsGeometry<T> & result,
                                            int numEval)
{
    result.clear();
    //Step 0: Read Input
    gsKnotVector<T> curve_knot = curve.knots();
    gsBasis<T> curve_basis = curve.basis();
    index_t m_xi = curve_knot.size(); //The number of knots of the curve
    index_t n_xi = curve_knot.size(); //The number of basis vectors
    index_t p_xi = curve.degree();       // Degree of polynomial

//    gsMatrix<T> Dom


    //Compute an increment
    index_t dXi = (curve_knot.last()- curve_knot.at(0))/numEval;
    //Initialize output array
    gsMatrix<> out_xyz(numEval, 3);


    //Points in u-direction
    gsVector<T> x_u = curve_knot.at(p_xi+1);

    gsMatrix<T> curveBasis, curveVal, curveDeriv, curveDeriv2, surfVal, surfDeriv, surfDeriv2;

    for(int j = 0; j < numEval;j++){
        //Find the span index where xi is located
        curve_knot.iFind(x_u) - curve_knot.begin();

        // Compute the basis functions at xi
        curve_basis.eval_into(x_u, curveBasis);
        curve.jacobian_into(x_u,curveDeriv);
        curve.deriv2_into(x_u, curveDeriv2);

        // Compute the Cartesian coordinates of xi in the parametric space of the surface
        curve.eval_into(x_u, curveVal);

        // Compute the Cartesian coordinates of xi in the physical space of the surface
        surface.eval_into(curveVal,surfVal);

        // Update the parametric location of the point
        x_u += dXi;
    }

}
}