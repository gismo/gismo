/** @file point_projection.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
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

    gsOptProblemExample(const gsGeometry<T> & g, const gsVector<T> & pt, const gsVector<T> & u)
    : m_g(&g), m_pt(pt), m_u(u)
    {
        // Initial guess: it is the input parameters u.
        // Number of design variables

        m_numDesignVars  = m_u.size(); // univariate or bivariate

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);


        gsMatrix<T> gsupp = m_g->support(); // d*2 matrix : [lower, upper]
        for(index_t i = 0; i < m_u.size(); i++)
        {
          m_desLowerBounds[i] = gsupp(i,0); // lower bound on the parameters
          m_desUpperBounds[i] = gsupp(i,1); // upper bound on the parameters
        }

        m_curDesign.resize(m_numDesignVars);
        m_curDesign << m_u; // the initial parameter

    }
    //! [OptProblemExample Constructor]

public:

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    T evalObj( const gsAsConstVector<T> & u ) const
    {
        // return the value of the objective function
        //gsInfo << "Starting guess:\n" << u << "\n";
        gsMatrix<T> proj;
        m_g->eval_into(u, proj);

        return 0.5 * (proj-m_pt).squaredNorm();
    }
    //! [OptProblemExample evalObj]

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        result.resize(m_numDesignVars,1);

        gsMatrix<T> tempParam = u;
        tempParam.reshape(m_g->support().rows(), m_g->support().cols());
        gsMatrix<T> tempEval;
        gsMatrix<T> tempGrad;
        gsMatrix<T> tempResult;

        m_g->eval_into(tempParam,tempEval);
        m_g->jacobian_into(tempParam,tempGrad);
        tempResult = tempGrad.transpose() * (tempEval - m_pt);
        result = tempResult.reshape(m_numDesignVars,1);

    }


private:
    const gsGeometry<T> * m_g;
    const gsVector<T> m_pt;
    const gsVector<T> m_u;
    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;
    //using gsOptProblem<T>::m_numConstraints;
    //using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    // using gsOptProblem<T>::m_conLowerBounds;
    // using gsOptProblem<T>::m_conUpperBounds;

    // using gsOptProblem<T>::m_conJacRows;
    // using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};
//! [OptProblem]




int main(int argc, char *argv[])
{
    std::string fn = "../filedata/surfaces/simple.xml";
    int numPatch = 0;
    int maxIter = 100;
    int mupdate = 1;


    gsCmdLine cmd("Point projection by HLBFGS method: https://xueyuhanlang.github.io/software/HLBFGS/");
    cmd.addString("g", "geometry", "Input geometry.", fn);
    cmd.addInt("p", "patch", "number of the patch to evaluate.", numPatch);
    cmd.addInt("i", "iter", "number of maximum iterations.", maxIter);
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd_in(fn);

    std::vector<gsGeometry<>::uPtr> geo = fd_in.getAll< gsGeometry<> >();
    if ( ! geo.empty() ){
        gsInfo<< "Got "<< geo.size() <<" patch"<<(geo.size() == 1 ? "." : "es.") <<"\n";
    }
    else{
        gsInfo<< "----------------------------------------------" <<"\n";
        gsInfo<< "No geometry has been found in input, quitting." <<"\n";
        gsInfo<< "----------------------------------------------" <<"\n";
        return 0;
    }

    if ( numPatch > geo.size()-1 ){
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        gsInfo<< "Patch number "<< numPatch <<" has not been found in input, quitting." <<"\n";
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        return 0;
    }

    gsMatrix<> u = geo[numPatch]->parameterCenter();
    gsMatrix<> tmp;
    geo[numPatch]->eval_into(u, tmp);
    gsInfo << "Parametric value:\n" << u << "\n";
    gsInfo << "Point on geometry:\n" << tmp << "\n";

    gsWriteParaview(*geo[numPatch], "ingeometry", 10000);
    gsWriteParaviewPoints(tmp, "pointOnGeo");

    gsMatrix<> mu = gsMatrix<>::Random(tmp.rows(), tmp.cols());
    mu = (mu + gsMatrix<>::Constant(tmp.rows(), tmp.cols(),1.))*0.1/2.;
    gsMatrix<> pt0(tmp.rows(), tmp.cols());
    pt0 = tmp + mu;
    gsWriteParaviewPoints(pt0, "pointToProject");

    gsMatrix<> u0(u.rows(), u.cols());
    u0 << 0.4, 0.4;


    gsOptProblemExample<real_t> problem(*geo[numPatch], pt0, u0);
    gsOptimizer<real_t> * optimizer;

    optimizer = new gsHLBFGS<real_t>(&problem);
    //
    //! [Optimizer options]
    // Set number of iterations as stop criterion.
    optimizer->options().setInt("MaxIterations",maxIter);
    //
    // Set verbosity
    optimizer->options().setInt("Verbose",2);
    //
    //
    optimizer->options().setReal("MinGradientLength", 1e-10);
    optimizer->options().setReal("MinStepLength",1e-10);
    optimizer->options().setInt("LBFGSUpdates", mupdate);
    //
    //
    //
    //! [Solve]
    // Start the optimization
    //gsVector<> in(curve.coefs().size() + params.size());
    //in << curve.coefs().reshape( curve.coefs().size() ,1), params.transpose();
    //optimizer->solve(in);
    optimizer->solve(problem.currentDesign());

    gsInfo << "\nNumber of iterations : " << optimizer->iterations() <<"\n";
    gsInfo << "Final objective value: " << optimizer->objective() <<"\n";
    gsInfo << "Final design:\n" << optimizer->currentDesign() <<"\n\n";

    gsMatrix<> ufinal = optimizer->currentDesign();
    gsMatrix<> ptfinal;
    geo[numPatch]->eval_into(ufinal, ptfinal);
    gsWriteParaviewPoints(ptfinal, "finalHLBFGpoint");

    gsVector<> u0vec(u0);
    gsInfo << "p0:\n" << pt0 << "\nu0:\n" << u0vec << "\n";

    geo[numPatch]->closestPointTo(pt0, u0vec, 1e-10, true);
    gsInfo << "gsGeometry::closestPointTo\n" << u0vec << "\n";
    gsMatrix<> ptfinalclosestPointTo;
    geo[numPatch]->eval_into(u0vec, ptfinalclosestPointTo);
    gsWriteParaviewPoints(ptfinalclosestPointTo, "finalCPTpoint");

    std::cout << "----------------------------------------\n" << std::scientific;
    gsInfo << "Difference in the parametric domain:"<< "\n";
    gsInfo << "      HLBFGS: " << (ufinal - u0).norm() << "\n";
    gsInfo << "         CPT: " << (u0vec - u0).norm() << "\n";
    gsInfo << "HLBFGS - CPT: " << (ufinal - u0vec).norm() << "\n";
    gsInfo << "Difference in the physical domain:" << "\n";
    gsInfo << "      HLBFGS: " << (ptfinal - pt0).norm() << "\n";
    gsInfo << "         CPT: " << (ptfinalclosestPointTo - pt0).norm() << "\n";
    gsInfo << "HLBFGS - CPT: " << (ptfinal - ptfinalclosestPointTo).norm() << "\n";



    return 0;
}
