/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>
#endif

#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsTimeIntegrator.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif


using namespace gismo;

int main(int argc, char *argv[])
{
#ifdef gsElasticity_ENABLED
#ifdef gsStructuralAnalysis_ENABLED

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    std::string precice_config;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe

    index_t Nsteps = 50;

    real_t dt = 1e-1;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe",method);
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "N", "Nsteps", "Number of steps",  Nsteps );
    cmd.addReal("t","dt","time step",dt);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::string methodName;
    if (method==1)
      methodName = "ExplEuler";
    else if (method==2)
      methodName = "ImplEuler";
    else if (method==3)
      methodName = "Newmark";
    else if (method==4)
      methodName = "Bathe";
    else if (method==5)
      methodName = "CentralDiff";
    else if (method==6)
      methodName = "RK4";
  
    //! [Read input file]
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-0.05,0.0,0.05,1.0));

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    for (index_t i = 0; i < numElevate; ++i)
        bases.degreeElevate();
    for (index_t i = 0; i < numRefine; ++i)
        bases.uniformRefine();

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";
    
    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    real_t mu = E / (2.0 * (1.0 + nu));
    real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Bottom side (fixed)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    // Neumann side
    gsVector<> neuData(2);
    neuData<<-1e4,0;
    gsConstantFunction<> g_N(neuData,patches.geoDim());
    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &g_N);
    //
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,2);

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(patches,bases,bcInfo,gravity);
    massAssembler.options().setReal("Density",rho);
    massAssembler.assemble();

    // creating stiffness assembler.
    gsElasticityAssembler<real_t> assembler(patches,bases,bcInfo,g);
    assembler.options().setReal("YoungsModulus",E);
    assembler.options().setReal("PoissonsRatio",nu);
    assembler.assemble();

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix().transpose();
    gsSparseMatrix<> K_T;

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) 
    {
        // to do: add time dependency of forcing
        assembler.assemble(x, fixedDofs);
        // For the shell
	//assembler.constructSolution(x, solution);
	//status = assembler.assembleMatrix(solution);
        m = assembler.matrix();
        return true;
    };


    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        assembler.assemble(x,fixedDofs);
	//Add assemble vector JL
	//assembler.constructSolution(x,solution);
        //status = assembler.assembleVector(solution);
    	result = assembler.rhs();
        return true;
    };

    //gsTimeIntegrator<real_t> timeIntegrator(M,K,Forcing,dt);
    gsTimeIntegrator<real_t> timeIntegrator(M,Jacobian,Residual,dt);
    timeIntegrator.verbose();
    timeIntegrator.setTolerance(1e-5);
    timeIntegrator.setMethod(methodName);

    //gsDebugVar(M.toDense());
    //gsDebugVar(K.toDense());


    size_t N = assembler.numDofs();
    gsMatrix<> uOld(N,1);
    gsMatrix<> vOld(N,1);
    gsMatrix<> aOld(N,1);
    uOld.setZero();
    vOld.setZero();
    aOld.setZero();

    gsMatrix<> uNew,vNew,aNew;
    uNew = uOld;
    vNew = vOld;
    aNew = aOld;
    timeIntegrator.setDisplacement(uNew);
    timeIntegrator.setVelocity(vNew);
    timeIntegrator.setAcceleration(aNew);

    timeIntegrator.constructSolution();
    gsParaviewCollection collection("./output/solution");
    real_t time = 0;

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

    gsStructuralAnalysisOutput<real_t> writer("pointData.csv",points);
    writer.init({"x","y"},{"time"});

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);
    for (index_t i=0; i<Nsteps; i++)
    {	        
    	gsStatus status = timeIntegrator.step();
        GISMO_ASSERT(status == gsStatus::Success,"Time integrator did not succeed");
    	time = timeIntegrator.currentTime();
        timeIntegrator.constructSolution();
        gsMatrix<> displacements = timeIntegrator.displacements();

        gsMultiPatch<> solution;
        assembler.constructSolution(displacements,fixedDofs,solution);

        // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);
        std::string fileName = "./output/solution" + util::to_string(i);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(i) + "0";
        collection.addTimestep(fileName,time,".vts");

        solution.patch(0).eval_into(points,pointDataMatrix);
        otherDataMatrix<<time;
        writer.add(pointDataMatrix,otherDataMatrix);

        // if (write)
        // {
        //   gsMatrix<> v(2,1);
        //   v<<  0.0,0.0;
        //   gsMatrix<> res2;
        //   solution.patch(0).eval_into(v,res2);
        //   time = timeIntegrator.currentTime();
        //   std::ofstream file;
        //   file.open(wn,std::ofstream::out | std::ofstream::app);
        //   file  << std::setprecision(10)
        //         << time << ","
        //         << res2(2,0) <<"\n";
        //   file.close();
        // }
        // Update solution with multipatch coefficients to generate geometry again
    }

    if (plot)
    {
      collection.save();
    }

    return  EXIT_SUCCESS;

#else
    GISMO_ERROR("This file requires gsElasticity and gsStructuralAnalysis");
    return EXIT_FAILURE;
#endif
#else
    GISMO_ERROR("This file requires gsElasticity");
    return EXIT_FAILURE;
#endif

}
