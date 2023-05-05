/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>

#include <gismo.h>

#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsGradientDescent.h>

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif


using namespace gismo;

// To define an optimization problem we inherit from gsOptProblem class
// and implement the default constructor and few inherited virtual functions

//! [OptProblemExample Class]
template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    //! [OptProblemExample Constructor]
    // The constructor defines all properties of our optimization problem
    gsOptProblemExample(const gsMultiPatch<T>& mp, const gsMatrix<T> & params, const gsMatrix<T> & X)
    :
    m_mp(mp),
    m_params(params),
    m_X(X)
    {
        // Initial guess: coefficients + parameters
        // Number of design variables
        m_numDesignVars  = m_mp[0].coefs().size() + m_params.size(); // dim * Dofs + numPts

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);


        for(index_t i = 0; i < m_mp[0].coefs().size(); i++)
        {
          m_desLowerBounds[i] = -1.0e19; // lower bound on the coefficients
          m_desUpperBounds[i] =  1.0e19; // upper bound on the coefficients
        }

        for(index_t i = 0; i < m_params.size(); i++)
        {
          m_desLowerBounds[m_mp[0].coefs().size() + i] = 0.; // lower bound on the parameters, take care or the finate differences
          m_desUpperBounds[m_mp[0].coefs().size() + i] = 1.; // upper bound on the parameters
          // m_desLowerBounds[m_mp[0].coefs().size() + i] = 0. + 3e-5; // lower bound on the parameters, take care or the finate differences
          // m_desUpperBounds[m_mp[0].coefs().size() + i] = 1. - 1e-5; // upper bound on the parameters
        }


        m_curDesign.resize(m_numDesignVars);
        m_curDesign << m_mp[0].coefs().reshape(m_mp[0].coefs().size(),1), m_params.transpose();

    }
    //! [OptProblemExample Constructor]

public:

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    T evalObj( const gsAsConstVector<T> & u ) const
    {
        gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within bounds
        u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
        // return the value of the objective function
        //gsInfo << "Starting guess:\n" << u << "\n";
        gsMatrix<T> currentcoefs( m_mp[0].coefs().rows(), m_mp[0].coefs().cols() );
        for(index_t i=0; i<m_mp[0].coefs().rows(); i++){
          currentcoefs(i,0) = u(i);
          currentcoefs(i,1) = u(i + m_mp[0].coefs().rows());
        }

        // gsInfo << "coefs 4 evaluation:\n" << currentcoefs << "\n";

        gsMatrix<T> currentparams(1, m_X.cols());
        for(index_t i = m_mp[0].coefs().size(); i < u.size(); i++){
          currentparams(0, i - m_mp[0].coefs().size() ) = u(i);
        }

        gsBSpline<T> current = dynamic_cast<gsBSpline<T>&>(*m_mp[0].basis().makeGeometry(give(currentcoefs)));
        //( bs, give(currentcoefs)); // make a multipatch object?

        gsMatrix<T> pred;
        current.eval_into(currentparams, pred);

        // gsInfo << "current curve:\n" << current << "\n";
        // gsInfo << "control points:\n" << current.coefs() << "\n";
        //gsInfo << "current prediction: " << pred.rows() << "x" << pred.cols() << "\n";
        //gsInfo << "       point cloud: " << m_X.rows() << "x" << m_X.cols() << "\n";

        T mse = 0.;
        // gsInfo << "Computing the mean square error:" << "\n";
        // gsInfo << "start: " << mse << "\n";
        for (index_t i = 0; i < m_X.cols(); i++)
        {
            T err = std::pow( (m_X.col(i) - pred.col(i)).norm(), 2) ;
            mse += err;
            // gsInfo << "diff (" << i << "):\n" << m_X.col(i) - pred.col(i) << "\n";
            // gsInfo << "err: " << err << "\n";
        }

        return 0.5 * mse;
    }
    //! [OptProblemExample evalObj]

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)
    // void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    // {
    //     gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within bounds
    //     u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
    //
    //     result.resize(m_numDesignVars,1);
    //     result.setZero();
    //
    //     gsMatrix<T> currentcoefs( m_mp[0].coefs().rows(), m_mp[0].coefs().cols() );
    //     for(index_t i=0; i<m_mp[0].coefs().rows(); i++)
    //     {
    //       currentcoefs(i,0) = u(i);
    //       currentcoefs(i,1) = u(i + m_mp[0].coefs().rows());
    //     }
    //
    //     // gsInfo << "coefs 4 evaluation:\n" << currentcoefs << "\n";
    //
    //     gsMatrix<T> currentparams(1, m_X.cols());
    //     for(index_t i = m_mp[0].coefs().size(); i < u.size(); i++){
    //       currentparams(0, i - m_mp[0].coefs().size() ) = u(i);
    //     }
    //
    //     // gsInfo << "params 4 evaluation: " << currentparams.rows() << " x " << currentparams.cols() << "\n";
    //     // gsInfo << currentparams << "\n";
    //
    //     gsBSpline<T> current = dynamic_cast<gsBSpline<T>&>(*m_mp[0].basis().makeGeometry(give(currentcoefs)));
    //     gsMatrix<T> deriv;
    //     current.deriv_into(currentparams, deriv);
    //     gsInfo << "first derivative: " << deriv.rows() << " x " << deriv.cols() << "\n";
    //
    //     gsMatrix<T> pred;
    //     current.eval_into(currentparams, pred);
    //
    //     gsMatrix<T> diff(pred.rows(), pred.cols());
    //     diff = pred - m_X;
    //
    //     gsInfo << "First indeces:\n";
    //     for(index_t el=0; el<m_mp[0].coefs().rows(); el++)
    //     {
    //       gsInfo << "indices: " << el << ", " << el + m_mp[0].coefs().rows() << "\n";
    //       for(index_t i=0; i<m_X.cols(); i++)
    //       {
    //         gsMatrix<T> basis_fun_val = m_mp[0].basis().evalSingle(el, currentparams.col(i));;
    //         result(el) += diff(0, i) * basis_fun_val(0,0);
    //         result(el + m_mp[0].coefs().rows()) += diff(1, i) * basis_fun_val(0,0);
    //       }
    //     }
    //
    //     gsInfo << "last indeces:\n";
    //     for(index_t i = m_mp[0].coefs().size(); i < u.size(); i++)
    //     {
    //       gsInfo << i << "\n";
    //       result(i) = diff.col(i - m_mp[0].coefs().size()).transpose() * deriv.col(i - m_mp[0].coefs().size());
    //     }
    //
    // }


private:
    const gsMultiPatch<T> m_mp; // the spline object
    const gsMatrix<T> m_params;
    const gsMatrix<T> m_X;
    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;


    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_curDesign;
};
//! [OptProblem]




int main(int argc, char *argv[])
{
    index_t numPts  = 100;
    index_t maxIter = 100;
    index_t mupdate = 5;
    index_t verbosity = 0;
    real_t minGradLen = 1e-6;


    gsCmdLine cmd("Planar B-spline curve fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");
    cmd.addInt("p", "psize", "number of planar points to sample.", numPts);
    cmd.addInt("i", "iter", "number of maximum iterations.", maxIter);
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addInt("v", "verbosity", "print verbose.", verbosity);
    cmd.addReal("t", "tolerance", "tolerance for the HLBFGS optimization.", minGradLen);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Data generation.
    // Make a BSpline curve
    int deg = 3;
    gsKnotVector<> kv(0, 1, 3, deg+1);//start,end,interior knots, start/end multiplicites of knots1
    gsBSplineBasis<> basis(kv);
    gsMatrix<> coefs(7, 2);
    coefs << 0, 0,
             0., 1,
             0.25, 1,
             0.5, 0,
             0.75, 1,
             1, 1,
             1, 0;

    gsBSpline<> original( kv, give(coefs));

    gsWriteParaview( original, "originalPlanarBsplineCurve", 1000);
    gsWriteParaview( original, "originalcnet", 1000, false, true);

    gsMatrix<> X;

    real_t umin = 0.;
    real_t umax = 1.;

    real_t urange = umax - umin; // umax (=1) - umin (=0)
    //gsInfo << "u-range sampling: " << urange << "\n";
    gsMatrix<> mu = gsMatrix<>::Random(1,numPts); // 1xnumPts Matrix filled with random numbers between (-1,1)
    mu = (mu + gsMatrix<>::Constant(1,numPts,1.))*urange/2.; // add umax to the matrix to have values between 0 and 2; multiply with range/2
    mu = (mu + gsMatrix<>::Constant(1,numPts,umin)); //set LO (=0) as the lower bound (offset)


    gsInfo << "sampling dimension: " << mu.rows() << " x " << mu.cols() << "\n";
    for(auto row : mu.rowwise())
      std::sort(row.begin(), row.end());
    //gsInfo << "Here is the sorted parameters:\n" << mu << "\n";



    // gsMatrix<> uline = gsVector<>::LinSpaced(numPts, 0, 1);
    // original.eval_into(uline.transpose(), X);

    original.eval_into(mu, X);
    gsMatrix<> xk = X.row(0);
    gsMatrix<> yk = X.row(1);

    //gsInfo << "Sampled points:\n" << X.rows() << " x " << X.cols() << "\n";
    gsWriteParaviewPoints(xk, yk, "inputPlanarPts");
    // gsInfo << "scattered parameters:\n" << mu << "\n";
    // gsInfo << "point cloud:\n" << X << "\n";

    // TODO: add noise as in http://dx.doi.org/10.1016/j.cagd.2012.03.004
    // gsMatrix<> X_noisy = X + 1e-3 * gsMatrix<>::Random(X.rows(), X.cols());
    // gsMatrix<> xn = X_noisy.row(0);
    // gsMatrix<> yn = X_noisy.row(1);
    //
    // //gsInfo << "Sampled points:\n" << X_noisy.rows() << " x " << X_noisy.cols() << "\n";
    // gsWriteParaviewPoints(xn, yn, "inputPlanarPts_noisy");

    // step 1. Specify an initial curve P(t), personal choice : (penalized) least square approximation
    // Assumption: knots are given and fixed (same as the one of the original B-spline curve).
    // gsMatrix<> P(basis_dim, 2); // try to use same notation of the paper. P = control points
    gsMatrix<> param_values(1,X.cols());
    param_values.row(0) = gsVector<>::LinSpaced(X.cols(), 0, 1);
    gsFitting<> fitted_curve(param_values, X, basis);
    // fitted_curve.compute(1e-7);
    fitted_curve.compute(0); // standard least square, to avoid the functional derivative computation
    gsBSpline<>curve(kv, give(fitted_curve.result()->coefs()));

    gsInfo << "Least square fitting, with smoothing:\n " << curve << "\n";
    gsWriteParaview( curve, "curve_step0", 1000);
    gsWriteParaview( curve, "cnet_step0", 1000, false, true);




    // step 2. Find the foot point P(t_k) on P(t) for evry data point X_k.
    // already as initial parameters, the foot point projection are provided.
    gsMatrix<> params(1, X.cols());
    params(0,0) = 0;
    params(0,X.cols()-1) = 1;
    for (index_t i = 1; i < X.cols()-1; ++i)
    {
        gsVector<> newParam;
        // curve.closestPointTo(X.col(i).transpose(), newParam, 1e-10, false);
        curve.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
        params.col(i) = newParam;
    }

    gsWriteParaviewPoints( params, "params_step0");

    gsInfo << "Fixed spline space:\n" << curve.basis() << "\n";
    gsInfo << "Initial coefficients:\n" << curve.coefs().rows() << " x " << curve.coefs().cols() << "\n";
    gsInfo << "Initial parameters:\n" << params.rows() << " x " << params.cols() << "\n";
    gsInfo << "x0 size: " << curve.coefs().size() + params.size()  << "\n";



    gsOptProblemExample<real_t> problem(curve, params, X); // FIRST THE COEFFICIENTS, THEN THE PARAMETERS AND THEN THE POINTS.
    gsOptimizer<real_t> * optimizer;

    optimizer = new gsHLBFGS<real_t>(&problem);
    //
    //! [Optimizer options]
    // Set number of iterations as stop criterion.
    optimizer->options().setInt("MaxIterations",maxIter);
    //
    // Set verbosity
    optimizer->options().setInt("Verbose", verbosity);

    // see the HLBFGS method for the options
    optimizer->options().setReal("MinGradientLength", minGradLen);
    // optimizer->options().setReal("MinStepLength",1);
    // optimizer->options().setInt("LBFGSUpdates", mupdate);


    //! [Solve]
    // Start the optimization
    gsVector<> in(curve.coefs().size() + params.size());
    in << curve.coefs().reshape( curve.coefs().size() ,1), params.transpose();
    optimizer->solve(in);
    // gsVector<> result = optimizer->currentDesign();
    // gsInfo << "Initial settings:\n" << result << "\n";
    // optimizer->solve(result);
    // //
    //! [Output]
    // Print final design info
    gsInfo << "\nNumber of iterations : " << optimizer->iterations() <<"\n";
    gsInfo << "Final objective value: " << optimizer->objective() <<"\n";
    //gsInfo << "Final design:\n" << optimizer->currentDesign() <<"\n";

    gsMatrix<> finaldesign = optimizer->currentDesign();

    gsMatrix<> currentcoefs( curve.coefs().rows(), curve.coefs().cols() );
    for(index_t i=0; i< curve.coefs().rows(); i++){
      currentcoefs(i,0) = finaldesign(i,0);
      currentcoefs(i,1) = finaldesign(i + curve.coefs().rows(),0);
    }

    //gsInfo << "coefs 4 evaluation:\n" << currentcoefs << "\n";

    gsMatrix<> currentparams(1, X.cols());
    for(index_t i = 2*curve.coefs().rows(); i<finaldesign.size(); i++){
      currentparams(0, i - 2*curve.coefs().rows() ) = finaldesign(i,0);
    }



    gsBSpline<> final( kv, give(currentcoefs));
    gsWriteParaview( final, "curve_final", 1000);


    delete optimizer;
    return EXIT_SUCCESS;
}
