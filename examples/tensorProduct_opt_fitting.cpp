/** @file tensorProduct_opt_fitting.cpp

    @brief Tensor product BSpline surface fitting with HLBFGS optimization.

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
#include <gsIpOpt/gsIpOpt.h>
#endif


using namespace gismo;

// This optimization problem addresses the optimization of coefficients and paramters
// for a tensor product bspline fitting model.
// Assumption: the spline space is fixed and does not change withing the optimization loop
//             the basis is encoded in the gsOptProblem member m_mp, i.e. m_mp->result()->basis()

// Input: gsFitting<T>& mp fitting object;
//        const gsMatrix<T> & params parameters in [0,1]^2;
//        const gsMatrix<T> & X points in [0,1]^3;

// Output: gsVector<T> u containing the optimized coefficients (cx, cy, cz) in R^3,  and parameters (u, v) in [0,1]^2.
//         dimension of the tensor-product spline space: n = m_mp->result()->basis().size(),
//         number of points m = X.cols()
//         then, u = [cx_0, cx_1, ... , cx_{n-1}, cy_1, cy_2, ... , cy_{n-1}, cz_1, cz_2, ... , cz_{n-1}, u_0, ..., u_{m-1}, v_0, ..., v_{m-1}]

//! [OptProblemExample Class]
template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    //! [OptProblemExample Constructor]
    // The constructor defines all properties of our optimization problem
    gsOptProblemExample(gsFitting<T>& mp, const gsMatrix<T> & params, const gsMatrix<T> & X)
    :
    m_mp(&mp),
    m_params(params),
    m_X(X)
    {
        // Number of design variables: how many variables we optimize, i.e. coefficiets + parametric values
        // m_numDesignVars  = m_mp[0].coefs().size() + m_params.size(); // dim * spline-dofs + numPts
        m_numDesignVars  = m_mp->result()->coefs().size() + m_params.size();

        // design bounds
        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // coefficiets in R^3, namely no bounds
        for(index_t i = 0; i < m_mp->result()->coefs().size(); i++)
        {
          m_desLowerBounds[i] = -1.0e19; // lower bound on the coefficients
          m_desUpperBounds[i] =  1.0e19; // upper bound on the coefficients
        }

        // parameters in [0,1]^2
        currentparams = m_params.transpose(); // m_mp.returnParamValues().transpose();
        for(index_t i = 0; i < currentparams.size(); i++)
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the parameters
        }

        // Initialization of the smoothing matrix that we need to define th objective function.
        m_G.resize(m_mp->result()->coefs().rows(), m_mp->result()->coefs().rows());
        m_G.reservePerColumn( cast<T,index_t>( (2 * m_mp->result()->basis().maxDegree() + 1) * 1.333 ) );
        m_mp->applySmoothing(m_mp->lambda(), m_G);

        // design variables: whant we do optimize.
        // c_x, c_y, c_z, u, v
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign << m_mp->result()->coefs().reshape(m_mp->result()->coefs().size(),1), currentparams.reshape(currentparams.size(),1);
    }
    //! [OptProblemExample Constructor]

public:

    //temporary storage
    mutable gsMatrix<T> tmp, currentparams;
    mutable gsSparseMatrix<T> c;
    mutable std::vector<gsSparseMatrix<T>> c_matrices;

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    // look at the manuscript for its rigurous definition
    // idea: 1/2 * ( (spline_model - points)^2 + lambda * smoothing_term )
    T evalObj( const gsAsConstVector<T> & u ) const
    {

        gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within design bounds
        u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);

        gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
        currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();

        c = m_mp->result()->basis().collocationMatrix( currentparams ) ;
        tmp.noalias() = c * currentcoefs - m_X.transpose();
        return 0.5 * ( (tmp * tmp.transpose()).trace() + (currentcoefs.transpose() * m_G * currentcoefs).trace() * m_mp->lambda() );
    }

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)

    //   compute the derivative with respect to the coefficients (cx, cy, cz)
    //   compute the partial derivative with respect to u-parameter
    //   compute the partial derivatice with respect to v-parameter
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
      // make sure that everything stays in bounds.
      gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() );
      u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);

      // finite differences, as if is left unimplemented; use to make the check if the implementation is correct.
      // gsOptProblem<T>::gradObj_into(u, check_result);

      gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
      currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();

      c_matrices = collocationMatrix1(m_mp->result()->basis(), currentparams);
      // c_matrices[0] : collocation matrix
      // c_matrices[1] : basis partial first derivatives
      tmp.noalias() = c_matrices[0] * currentcoefs - m_X.transpose();

      result.head(currentcoefs.size()).noalias() = // d_coefs.asVector();
      ( ( c_matrices[0].transpose()*c_matrices[0] + m_mp->lambda() * m_G ) * currentcoefs - c_matrices[0].transpose() * m_X.transpose() )
          .reshaped(currentcoefs.size(),1);

      result.middleRows(currentcoefs.size(),m_X.cols()).noalias() = (tmp * (c_matrices[1] * currentcoefs).transpose()).diagonal();
      result.tail(m_X.cols()).noalias() = (tmp * (c_matrices[2] * currentcoefs).transpose()).diagonal();
    }



private:

    gsFitting<T> *m_mp;
    const gsMatrix<T> m_params;
    const gsMatrix<T> m_X;
    gsSparseMatrix<T> m_G;

    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_curDesign;

};
//! [OptProblem]

int main(int argc, char *argv[])
{

    index_t maxIter = 1;
    index_t mupdate = 20;
    real_t deg = 2;
    index_t numKnots = 2;
    real_t lambda = 1e-6;
    real_t gtoll = 1e-7; // to decrees to push trogh the iterations
    std::string fn = "../filedata/fitting/floaterPts_out.xml";
    index_t verbosity = 0; // 0 (no videoprint), 1 (some videoprint), 2 (a lot of videoprint)
    bool ptype = false; // keep it false.
    bool plotInParaview = false;
    index_t plotIt = -1;
    bool apdm = false;

    gsCmdLine cmd("Tensor product B-spline surface fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");

    cmd.addInt("i", "iter", "number of maximum iterations of the optimization algorithm(s).", maxIter);
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addReal("d", "degree", "bi-degree (q,q).", deg);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addInt("b", "print", "set printing verbosity", verbosity);
    cmd.addReal("g", "gtoll", "stopping criteria on ||g||", gtoll);
    cmd.addSwitch("p", "parameters", "input parameters: (0) from .xml file; (1) for foot-point projection;", ptype);
    cmd.addSwitch("z", "plot", "(0): no paraview plot generated.", plotInParaview);
    cmd.addInt("k", "kplot", "iteration of the optimization procedure to be plotted", plotIt);
    cmd.addSwitch("a", "apdm", "run the A-PDM algorithm.", apdm);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsStopwatch time;

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, X;
    fd_in.getId<gsMatrix<> >(1, uv );
    fd_in.getId<gsMatrix<> >(0, X);
    //! [Read data]

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, numKnots, deg+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, numKnots, deg+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> basis( u_knots, v_knots );

    gsFitting<real_t> fitting_object(uv, X, basis);
    fitting_object.compute(lambda);

    gsMatrix<> coefs(basis.size(), 3);
    coefs = fitting_object.result()->coefs();

    gsTensorBSpline<2, real_t> original(basis, coefs);


    // step 2. Find the foot point P(t_k) on P(t) for evry data point X_k.
    // already as initial parameters, the foot point projection are provided.
    // REMARK: we can start with any different parameterization, not necessary with the foot-point projection.
    gsMatrix<> params(2, X.cols());
    for (index_t i = 0; i < X.cols(); ++i)
    {
        if (ptype)
        {
          gsVector<> newParam;
          original.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
          params.col(i) = newParam;
        }
        else
        {
          params.col(i) = uv.col(i);
        }
    }

    gsInfo << "Fixed spline space:\n" << original.basis() << "\n";
    gsInfo << "Initial coefficients:\n" << original.coefs().rows() << " x " << original.coefs().cols() << "\n";
    gsInfo << "Initial parameters:\n" << params.rows() << " x " << params.cols() << "\n";
    gsInfo << "x0 size: " << original.coefs().size() + params.size()  << "\n";

    gsInfo << "Starting geometry:\n " << original << "\n";

    if (plotInParaview)
    {
      gsWriteParaview( original, "geo_it0", 1000, true, false);
      gsWriteParaview( original, "cnet_it0", 1000, false, true);
      gsMatrix<> originalCoefsToPlot(original.coefs().cols(), original.coefs().rows());
      originalCoefsToPlot = original.coefs().transpose();
      gsWriteParaviewPoints( originalCoefsToPlot, "coefs_it0");
      gsWriteParaviewPoints( params, "params_it0");
    }

    real_t rmse0 = 0.;
    // compute the fitting error.
    gsMatrix<> tmp0 = original.eval(params) - X;
    real_t pred0_eval = (tmp0 * tmp0.transpose()).trace();
    rmse0 += math::pow(pred0_eval/X.cols(), 0.5);

    std::cout << "Initial fitting error: " << rmse0 << std::scientific << "\n";

    // Initialization of the optimization problem.
    gsOptProblemExample<real_t> problem(fitting_object, params, X); // FIRST THE fitting object, THEN THE PARAMETERS AND THEN THE POINTS.
    gsOptimizer<real_t> * optimizer;

    // to store the information in .csv files.
    std::ofstream file_opt, file_pc;
    file_opt.open("results_CPDM.csv");
    file_opt << "it, time, rmse\n";


    gsInfo << "Fast fitting with HLBFGS:\n";
    // TODO: store time and fitting error in a proper way.
    // Uncomment the following line to avoid the foor loop on the maximum number of iterations.
    // index_t it_opt = maxIter;
    for(index_t it_opt = 1; it_opt <= maxIter; it_opt++)
    {
      optimizer = new gsHLBFGS<real_t>(&problem);
      optimizer->options().setInt("Verbose",verbosity);
      optimizer->options().setReal("MinGradientLength", gtoll);
      optimizer->options().setInt("LBFGSUpdates", mupdate);

      //! [Solve]
      // Start the optimization
      gsVector<> in(original.coefs().size() + params.size());
      gsMatrix<> t_params(params.cols(), params.rows());
      t_params = params.transpose();
      in << original.coefs().reshape( original.coefs().size() ,1), t_params.reshape( t_params.size() ,1);
      optimizer->options().setInt("MaxIterations",it_opt); // set maximum number of iterations
      time.restart(); // start optimization algorithm
      optimizer->solve(in);
      real_t finaltime = time.stop(); // end optimization algorithm

      // assemble the new geometry with optimized coefficiets and parameters
      gsMatrix<> finaldesign = optimizer->currentDesign();
      gsMatrix<> currentcoefs( original.coefs().rows(), original.coefs().cols() );
      for(index_t i=0; i < original.coefs().rows(); i++)
      {
        currentcoefs(i,0) = finaldesign(i,0);
        currentcoefs(i,1) = finaldesign(i + original.coefs().rows(),0);
        currentcoefs(i,2) = finaldesign(i + 2*original.coefs().rows(),0);
      }

      gsMatrix<> currentparams(2, X.cols());

      for(index_t i = original.coefs().size(); i < finaldesign.size(); i++)
      {
        index_t poff = X.cols() + original.coefs().size();
        if ( i < X.cols() + original.coefs().size())
        {
          currentparams(0, i - original.coefs().size() ) = finaldesign(i);
        }
        else
        {
          currentparams(1, i - poff ) = finaldesign(i);
        }

      }

      if(plotInParaview && it_opt == plotIt) // plot certain output
      {
        gsWriteParaviewPoints( currentparams, "params_CPDM_it" + internal::to_string(it_opt));
        gsMatrix<> currentCoefsToPlot(currentcoefs.cols(), currentcoefs.rows());
        currentCoefsToPlot = currentcoefs.transpose();
        gsWriteParaviewPoints( currentCoefsToPlot, "coefs_CPDM_it" + internal::to_string(it_opt));

        // gsTensorBSpline<2, real_t> final( basis, give(currentcoefs));
        gsTensorBSpline<2, real_t> final( basis, currentcoefs);
        gsWriteParaview( final, "geo_CPDM_it" + internal::to_string(it_opt), 1000, true, false);
        gsWriteParaview( final, "cnet_CPDM_it" + internal::to_string(it_opt), 1000, false, true);
      }

      gsTensorBSpline<2, real_t> currentGeo(basis, currentcoefs);
      real_t rmse = 0.;
      // compute the fitting error.
      gsMatrix<> tmp = currentGeo.eval(currentparams) - X;
      real_t pred_eval = (tmp * tmp.transpose()).trace();
      rmse += math::pow(pred_eval/X.cols(), 0.5);
      //gsInfo << it_opt << ", " << time << ", " << rmse << "\n";
      file_opt <<std::setprecision(3)<< std::to_string(it_opt)<<std::setprecision(12) << "," << std::to_string(finaltime) << "," << std::to_string(rmse) << "\n";

      gsInfo << it_opt << ", " << time << ", " << rmse <<"\n";
      if (verbosity > 1)
      {
        gsInfo << "\nNumber of iterations : " << optimizer->iterations() <<"\n";
        gsInfo << "Final objective value: " << optimizer->objective() <<"\n";
        gsInfo<<"Fitting time: "<< time <<"\n";
        // gsInfo << "Final design:\n" << optimizer->currentDesign() <<"\n"; // this plot the whole vector given output from the optimizer.
        gsInfo << "params are moved from the originals by: " << (currentparams - params).norm() << "\n";
        gsInfo << "coefficients are moved from the originals by: " << (currentcoefs - original.coefs()).norm() << "\n";
      }
    } // maxIter
    file_opt.close();



    if(apdm)
    {
      gsInfo << "Running the A-PDM algorithm for comparion.\n";
      file_pc.open("results_APDM.csv");
      file_pc << "it, time, rmse\n";
      gsFitting<real_t> ref(uv, X, basis); // original geometry, same starting point for C-PDM;
      ref.compute(lambda);
      if(verbosity > 1)
        gsInfo << "it     time     rmse\n";
      // index_t step = maxIter; // uncomment to avoid foor loop on maximum number of iterations.
      for (index_t step=1; step <= maxIter; step ++)
      {
        time.restart(); // start optimization procedure: 1 step = points projection + refit to update the control points.
        ref.parameterCorrection(1e-7, step, 1e-4); //closestPoint accuracy, orthogonality tolerance
        real_t finaltime = time.stop(); // end of the optimization algorithm

      if(plotInParaview && step == plotIt) // plot certain output
      {
        gsWriteParaviewPoints( ref.returnParamValues(), "params_APDM_it" + internal::to_string(step));
        gsMatrix<> currentCoefsToPlot(ref.result()->coefs().cols(), ref.result()->coefs().rows());
        currentCoefsToPlot = ref.result()->coefs().transpose();
        gsWriteParaviewPoints( currentCoefsToPlot, "coefs_APDM_it" + internal::to_string(step));

        gsWriteParaview( *ref.result(), "geo_APDM_it" + internal::to_string(step), 1000, true, false);
        gsWriteParaview( *ref.result(), "cnet_APDM_it" + internal::to_string(step), 1000, false, true);
      }

      // fitting error
      real_t rmse = 0.;
      gsMatrix<> tmp = ref.result()->eval(ref.returnParamValues()) - X;
      real_t pred_eval = (tmp * tmp.transpose()).trace();
      rmse += math::pow(pred_eval/X.cols(), 0.5);

      gsInfo << step << ", " << time << ", " << rmse <<"\n";

      // store data in .csv file
      file_pc << std::to_string(step) << "," << std::to_string(finaltime) << "," << std::to_string(rmse) << "\n";

      } // maxIter
      file_pc.close();
    } // fi apdm

    return EXIT_SUCCESS;

}
