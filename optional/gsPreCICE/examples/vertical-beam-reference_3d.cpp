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
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif


using namespace gismo;


int main(int argc, char *argv[])
{
#ifdef gsElasticity_ENABLED
#ifdef gsStructuralAnalysis_ENABLED

    //! [Parse command line]
    bool plot = false;
    bool write = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe

    index_t Nsteps = 50;

    real_t dt = 1e-1;

    std::string dirname = "output";

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe",method);
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "N", "Nsteps", "Number of steps",  Nsteps );
    cmd.addReal("t","dt","time step",dt);
    cmd.addString("o", "outputDir","Output directory",dirname);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Create a file with point data", write);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineCube());

    real_t L, B, H;
    L = 0.1;
    B = 0.5;
    H = 1.0;
    patches.patch(0).coefs().col(0).array() *= B;
    patches.patch(0).coefs().col(1).array() *= H;
    patches.patch(0).coefs().col(2).array() *= L;
    patches.patch(0).coefs().col(2).array() -= L/2;

    // Make basis
    gsMultiBasis<> bases(patches);
    gsDebugVar(bases);

    for (index_t i = 0; i < numElevate; ++i)
        bases.degreeElevate();
    for (index_t i = 0; i < numRefine; ++i)
        bases.uniformRefine();

    gsInfo << "Patches: "<< patches.nPatches() <<", max degree: "<< bases.minCwiseDegree() <<"\n";
    for (size_t p=0; p!=patches.nPatches(); p++)
        gsInfo<<"Basis "<<p<<":\n"<<patches.basis(p)<<"\n";

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Bottom side (fixed)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 2);

    // bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr,2);
    // Neumann side
    gsVector<> neuData(3);
    neuData<<0,0,-1e4;
    gsConstantFunction<> g_N(neuData,patches.geoDim());

    bcInfo.addCondition(0, boundary::front, condition_type::neumann, &g_N);
    //
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    // source function, rhs
    gsFunctionExpr<> g("0","0","0",3);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,0.,3);

    gsMultiPatch<> mp_def = patches;




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

    gsSparseMatrix<> C = gsSparseMatrix<>(assembler.numDofs(),assembler.numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };

    gsDynamicBase<real_t> * timeIntegrator;
    if (method==1)
        timeIntegrator = new gsDynamicExplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==2)
        timeIntegrator = new gsDynamicImplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==3)
        timeIntegrator = new gsDynamicNewmark<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==4)
        timeIntegrator = new gsDynamicBathe<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==5)
    {
        timeIntegrator = new gsDynamicWilson<real_t,true>(Mass,Damping,Jacobian,Residual);
        timeIntegrator->options().setReal("gamma",1.4);
    }
    else if (method==6)
        timeIntegrator = new gsDynamicRK4<real_t,true>(Mass,Damping,Jacobian,Residual);
    else
        GISMO_ERROR("Method "<<method<<" not known");

    timeIntegrator->options().setReal("DT",dt);
    timeIntegrator->options().setReal("TolU",1e-3);
    timeIntegrator->options().setSwitch("Verbose",true);

    size_t N = assembler.numDofs();
    gsVector<> U(N), V(N), A(N);
    U.setZero();
    V.setZero();
    A.setZero();

    gsFileManager::mkdir(dirname);

    gsParaviewCollection collection(dirname + "/solution");
    real_t time = 0;

    gsMatrix<> points(3,1);
    points.col(0)<<0.5,1,0.0;

    gsStructuralAnalysisOutput<real_t> writer(dirname + "/pointData.csv",points);
    writer.init({"x","y","z"},{"time"});

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);
    for (index_t i=0; i<Nsteps; i++)
    {           
        gsStatus status = timeIntegrator->step(time,dt,U,V,A);
        GISMO_ASSERT(status == gsStatus::Success,"Time integrator did not succeed");
        time += dt;

        if (plot||write)
        {
            // Update the displacement vector
            gsMultiPatch<> solution;
            assembler.constructSolution(U,fixedDofs,solution);
            if (plot)
            {
                gsField<> solField(patches,solution);
                std::string fileName = dirname + "/solution" + util::to_string(i);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = "solution" + util::to_string(i) + "0";
                collection.addTimestep(fileName,time,".vts");
            }
            if (write)
            {
                solution.patch(0).eval_into(points,pointDataMatrix);
                otherDataMatrix<<time;
                writer.add(pointDataMatrix,otherDataMatrix);
            }
        }
    }

    if (plot)
      collection.save();

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
