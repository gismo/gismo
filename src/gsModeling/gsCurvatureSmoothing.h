/** @file gsCurvatureSmoothing.h

    @brief Computes a closed B-spline curve with a smaller number of curvature
    extrema compared to a given closed B-spline curve  i.e. some kind of
    smoothing the curvature of the curve. This smoothting can be done with the
    help of two methods - total variation and Hadenfeld's algorithm (see Jan
    Hadenfeld, Iteratives Glätten von B-Spline Kurven und B-Spline Flächen,
    Shaker Verlag, PhD Thesis)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl
*/

#pragma once

#include <iostream>
#include <stdexcept>

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsBSpline.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{
/**
    @brief Class for computing a closed B-spline curve with a smaller number of curvature extrema compared to a given closed B-spline curve

    i.e. some
    kind of smoothing the curvature of the curve with the help of two different methods

    \ingroup Modeling
*/

template<class T>
class gsCurvatureSmoothing
{
public:
    /// default constructor
    gsCurvatureSmoothing()
    {
        m_curve_original = NULL;
        m_curve_smooth = NULL;
    }

    /// constructor
    gsCurvatureSmoothing(const gsBSpline<T> & init_curve, const gsMatrix<T>& param_values, const gsMatrix<T>& points)
    {
        m_curve_original = &init_curve;
        m_curve_smooth = init_curve.clone().release();
        m_param_values = param_values;
        m_points = points;
    }

    /// Destructor
    ~gsCurvatureSmoothing()
    {
        delete m_curve_smooth; // deleting the B-spline curve
    }


    /// gives back the original B-spline curve
    const gsBSpline<T> & curveOriginal() const { return *m_curve_original; }
    /// gives back the smoother B-spline curve
    const gsBSpline<T> & curveSmooth() const { return *m_curve_smooth; }


    /// smooth the curve by total variation -- computes the stepsize by itsself (with the help of a backtracking line search method)
    /// this method should be used instead of the two methods void smoothTotalVariationSelectLamda
    void smoothTotalVariation(const T omega1, const T omega2, const T lamda, const T tau, const unsigned iter=50);

    /// smooth the curve by total variation -- uses different stepsizes (in listlamdas) in the gradient descent method (look for the best from the list!)
    /// if possible use the method void smoothTotalVariation
    void smoothTotalVariationSelectLamda(const T omega1, const T omega2, const gsMatrix<T> listlamdas, const unsigned iter=50);

    /// smooth the curve by total variation -- uses always the same stepsize lamda - but which has to be chosen!!
    /// if possible use the method void smoothTotalVariation
    void smoothTotalVariationSelectLamda(const T omega1, const T omega2, const T lamda, const unsigned iter=50);

    /// smooth the curve by smoothing only one cofficient in each step using the Hadenfeld algorithm --- the usual Hadenfeld algorithm --this method should be used
    void smoothHadenfeld(const unsigned smooth_degree, const T delta, const index_t iter_step, const index_t iter_total, gsVector<index_t> &iterated, const bool original=true);


    /// smooth the curve in one step for all coefficients using the Hadenfeld algorithm. Be aware of the fact
    /// that it is not ensured that we get a nice result --- can work but do not have to work (not sure that it will converge!!)
    /// if possible use method void smoothHadenfeld
    void smoothAllHadenfeld(const unsigned smooth_degree=4, const unsigned iter=500);

    /// writes the smooth curve to a file, which can be visualized in Mathematica (Name of Mathematica File VisualizationCurvatureSmoothing)
    void write(std::ostream &os);

    /// computes approximation error of the smoother curve to the original point cloud
    void computeApproxError(T & error);

    /// computes the L^{2}-norm approximation error of the smoother curve
    void computeApproxErrorL2(T & error);

    /// computes the L-max-norm approximation error of the smoother curve
    void computeApproxErrorLMax(T & error);

    /// computes the max approximation error between the coeffeicientsof the original and the smoother curve
    void computeApproxErrorCoef(T & error);

    /// computes the curvature error of the smoother curve
    void computeCurvatureError(T & error);

private:

    /// the original B-spline curve
    const gsBSpline<T> * m_curve_original;
    /// the smoother B-spline curve
    gsBSpline<T> * m_curve_smooth;
    /// the parameter values of the original point cloud
    gsMatrix<T> m_param_values;
    /// the points of the original point cloud
    gsMatrix<T> m_points;

    /// computes all values and derivatives (up to three) at the parameter values u for the given coefs
    void compute_AllValues(gsBSplineBasis<T> * basis, gsMatrix<T> u, gsMatrix<T> *coefs, gsMatrix<T> & values0, gsMatrix<T> & values1, gsMatrix<T> & values2, gsMatrix<T> & values3);

    /// computes the objective function for given coefs and omega1 and omega2 -- objective function = omega1*ApproximationFunction + omega2*CurvatureFunction
    void compute_ObjectiveFunction(gsBSplineBasis<T> * basis, gsMatrix<T> *coefs, const T omega1, const T omega2, T &value);

    /// set the smooth curve to the the original curve
    // TODO: What is exactly the purpose of this function; why do we want to change output?
    void reset(gsBSpline<T> * newCurve)
    {
        if (m_curve_smooth)
            delete m_curve_smooth;
        m_curve_smooth = newCurve;
    }

}; // class gsCurvatureSmoothing

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCurvatureSmoothing.hpp)
#endif
