/** @file adaptive_opt_fitting_example.cpp

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
    gsOptProblemExample(gsHFitting<2, T>& mp, const gsMatrix<T> & params, const gsMatrix<T> & X)
    :
    m_mp(&mp),
    m_params(params),
    m_X(X)
    {
        // Initial guess: coefficients + foot points
        // Number of design variables
        // m_numDesignVars  = m_mp[0].coefs().size() + m_params.size(); // dim * Dofs + numPts
        m_numDesignVars  = m_mp->result()->coefs().size() + m_params.size(); // dim * Dofs + numPts
        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        for(index_t i = 0; i < m_mp->result()->coefs().size(); i++)
        {
          m_desLowerBounds[i] = -1.0e19; // lower bound on the coefficients
          m_desUpperBounds[i] =  1.0e19; // upper bound on the coefficients
        }

        gsMatrix<> t_params(m_params.cols(), m_params.rows());
        t_params = m_params.transpose();
        // gsInfo << "parameters order:\n" << t_params << "\n";
        for(index_t i = 0; i < t_params.size(); i++)
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the parameters, take care or the finate differences
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the parameters
        }

        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign << m_mp->result()->coefs().reshape(m_mp->result()->coefs().size(),1), t_params.reshape(t_params.size(),1); // coefs_x, coefs_y, coefs_z, u, v

        m_G.resize(m_mp->result()->coefs().rows(), m_mp->result()->coefs().rows());
        gsInfo << "Smoothing matrix dimensions: " << m_G.rows() << " x " << m_G.cols() << "\n";
        m_mp->applySmoothing(m_mp->lambda(), m_G);

    }
    //! [OptProblemExample Constructor]

public:

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    T evalObj( const gsAsConstVector<T> & u ) const
    {

        gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within bounds
        u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
        gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
        gsAsConstMatrix<T> currentparams(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2);

        gsMatrix<T> tmp, val;
        gsMatrix<index_t> act;
        T f_eval = 0.;
        gsSparseMatrix<T> c = m_mp->result()->basis().collocationMatrix(gsMatrix<T>(currentparams.transpose() ) ) ;
        tmp = c * currentcoefs - m_X.transpose();
        f_eval = (tmp * tmp.transpose()).trace();
        f_eval += (currentcoefs.transpose() * m_G * currentcoefs).trace() * m_mp->lambda();

        return 0.5 * f_eval ;
    }

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
    // make sure that everything stays in bounds.
      gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() );
      u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);

      gsAsConstMatrix<T> currentparams(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2);
      gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());

      std::vector<gsSparseMatrix<T>> c_matrices = collocationMatrix1(m_mp->result()->basis(), gsMatrix<T>(currentparams.transpose()));

      gsMatrix<T> d_coefs = ( c_matrices[0].transpose()*c_matrices[0] + m_mp->lambda() * m_G ) * currentcoefs - c_matrices[0].transpose() * m_X.transpose();
      gsMatrix<T> tmp = c_matrices[0] * currentcoefs - m_X.transpose();
      gsMatrix<T> d_u = (tmp * (c_matrices[1] * currentcoefs).transpose()).diagonal();
      gsMatrix<T> d_v = (tmp * (c_matrices[2] * currentcoefs).transpose()).diagonal();

      result.head(currentcoefs.size()) = d_coefs.asVector();
      result.middleRows(currentcoefs.size(),d_u.size()) = d_u;
      result.tail(d_v.size()) = d_v;
    }



private:

    gsHFitting<2, T> *m_mp;
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

    index_t maxIter = 1; // for the optimization algorithm
    index_t refIter = 1; // for the adaptive loop
    index_t numURef = 0;
    index_t mupdate = 20;
    real_t deg = 2;
    index_t numKnots = 2;
    real_t lambda = 1e-6;
    real_t gtoll = 1e-9;
    std::string fn = "../filedata/fitting/floaterPts_out.xml";
    index_t verbosity = 1;
    index_t extension = 2;
    real_t threshold = 1e-02;
    real_t tolerance = 1e-02;

    gsCmdLine cmd("Tensor product B-spline surface fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");

    cmd.addInt("i", "iter", "number of maximum optmizing iterations.", maxIter);
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addReal("d", "degree", "bi-degree (q,q).", deg);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addInt("b", "print", "set printing verbosity", verbosity);
    cmd.addReal("g", "gtoll", "stopping criteria on ||g||", gtoll);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsStopwatch time;

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, X;
    fd_in.getId<gsMatrix<> >(1, uv );
    fd_in.getId<gsMatrix<> >(0, X);
    //! [Read data]

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");

    GISMO_ASSERT( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, numKnots, deg+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, numKnots, deg+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create Initial hierarchical basis
    gsTHBSplineBasis<2>  basis ( T_tbasis ) ;

    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    gsHFitting<2, real_t> opt_f(uv, X, basis, 0, ext, lambda);
    const std::vector<real_t> & errors = opt_f.pointWiseErrors();
    std::vector<real_t> errors2;
    real_t sum_of_errors2;


    for(int i = 0; i <= refIter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo << "C-PDM for the adaptive loop.\n";
        gsInfo<<"Adaptive loop iteration .. "<<i<<".."<<"\n";

        time.restart();

        opt_f.nextIteration(tolerance, threshold, 0); // no parameter correction here, but the combined approach.

        // after the fitting, we optimize the control points and the parametrica values.

        gsMatrix<> coefs(basis.size(), 3);
        coefs = opt_f.result()->coefs();

        gsMatrix<> params(2, X.cols());
        params = opt_f.returnParamValues();


        // Idea: you could initialize the gsOptProblem in a different way.
        gsOptProblemExample<real_t> problem(opt_f, params, X); // FIRST THE COEFFICIENTS, THEN THE PARAMETERS AND THEN THE POINTS.
        gsOptimizer<real_t> * optimizer;

        gsGeometry<>::uPtr tmp = (opt_f.result()->basis().makeGeometry(opt_f.result()->coefs()));
        gsTHBSpline<2,real_t> & original = dynamic_cast<gsTHBSpline<2,real_t> &>(*tmp);
        // gsTHBSpline<2, real_t> original(opt_f.result()->basis(), opt_f.result()->coefs());

        gsMatrix<> currentparams(2, X.cols());
        gsMatrix<> currentcoefs( original.coefs().rows(), original.coefs().cols() );



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
        optimizer->options().setInt("MaxIterations",maxIter);
        optimizer->solve(in);

        gsMatrix<> finaldesign = optimizer->currentDesign();

        for(index_t i=0; i < original.coefs().rows(); i++)
        {
          currentcoefs(i,0) = finaldesign(i,0);
          currentcoefs(i,1) = finaldesign(i + original.coefs().rows(),0);
          currentcoefs(i,2) = finaldesign(i + 2*original.coefs().rows(),0);
        }

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

        opt_f.updateGeometry(currentcoefs, currentparams);

        // compute mean squared error
        opt_f.get_Error(errors2, 0);
        sum_of_errors2 = std::accumulate(errors2.begin(), errors2.end(), 0.0);

        gsInfo<<"Fitting time: "<< time <<"\n";

        index_t dofs = opt_f.result()->basis().size();
        real_t minPointError = opt_f.minPointError();
        real_t maxPointError = opt_f.maxPointError();
        real_t mseError = sum_of_errors2/errors2.size();
        real_t percentagePoint = 100.0 * opt_f.numPointsBelow(tolerance)/errors.size();

        gsInfo<<"Fitted with "<< opt_f.result()->basis() <<"\n";
        gsInfo    << "DOFs         : "<< dofs <<"\n";
        std::cout << "Min distance : "<< minPointError << std::scientific <<"\n";
        std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
        std::cout << "         MSE : "<< mseError << std::scientific <<"\n";
        gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

        if ( opt_f.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<< i <<" iterations.\n";
            break;
        }
    } // iterative loop

    return EXIT_SUCCESS;

}
