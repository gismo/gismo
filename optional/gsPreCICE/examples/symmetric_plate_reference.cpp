/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
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
#ifdef gsKLShell_ENABLED
#ifdef gsStructuralAnalysis_ENABLED

    //! [Parse command line]
    bool plot = false;
    bool write = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson

    index_t Nsteps = 100;

    real_t dt = 1e-1;

    std::string dirname = "./output_ball";

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson",method);
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "N", "Nsteps", "Number of steps",  Nsteps );
    cmd.addReal("t","dt","time step",dt);
    cmd.addString("o", "outputDir","Output directory",dirname);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Create a file with point data", write);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    gsMultiPatch<> patches;
    gsMultiPatch<> solutions;
    // real_t L, B, H;
    patches.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,0));
    patches.addAutoBoundaries();
    patches.embed(3);

    // Create bases
    //true: poly-splines (not NURBS)

    for (index_t i = 0; i < numElevate; ++i)
        patches.degreeElevate();
    for (index_t i = 0; i < numRefine; ++i)
        patches.uniformRefine();

    gsMultiBasis<> bases(patches);

    gsInfo << "Patches: "<< patches.nPatches() <<", max degree: "<< bases.minCwiseDegree() <<"\n";
    for (size_t p=0; p!=patches.nPatches(); p++)
        gsInfo<<"Basis "<<p<<":\n"<<patches.basis(p)<<"\n";

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
                                                 
    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Bottom side (fixed)
    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 0);
    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 2);

    // bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr, 2);
    // bcInfo.addCondition(0, boundary::north, condition_type::clamped, nullptr, 2);
    // bcInfo.addCondition(0, boundary::east, condition_type::clamped, nullptr, 2);
    // bcInfo.addCondition(0, boundary::west, condition_type::clamped, nullptr, 2);

    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    // bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, nullptr, 0);
    // bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, nullptr, 0);
    // bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 1);

    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet,0,0,false, 0);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet,0,0, false, 1);
    bcInfo.addCondition(0, boundary::north, condition_type::clamped,0,0, false, 2);
    bcInfo.addCondition(0, boundary::west, condition_type::clamped,0,0, false, 2);    // bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);



    
    gsVector<> point(2);
    point << 0, 0;

    gsVector<> loadVec(3);
    loadVec<<0,0, -5e3;

    // gsVector<> loadVec(3);
    // loadVec<<0,1,-5e2;
    // gsConstantFunction<> surfForce(loadVec,3);
    pLoads.addLoad(point, loadVec, 0 );
    gsFunctionExpr<> force("0","0","-5e3",3);

    // Assign geometry map
    bcInfo.setGeoMap(patches);


    gsDebugVar(bcInfo);

// ----------------------------------------------------------------------------------------------

    // source function, rhs
    // gsConstantFunction<> g(0.,0.,0.,3);

    // Set up the material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E),3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu),3);
    gsFunctionExpr<> Density(std::to_string(rho),3);

    // Define thickness
    real_t thickness = 0.1;

    gsFunctionExpr<> t(std::to_string(thickness),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;

    gsMaterialMatrixBase<real_t>* materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);

    materialMatrix = getMaterialMatrix<3,real_t>(patches,t,parameters,Density,options);

    gsFunctionExpr<> foundation("0","0","1e3",3);

    gsThinShellAssembler<3, real_t, true> assembler(patches,bases,bcInfo,force,materialMatrix);
    assembler.setFoundation(foundation);
    assembler.setPointLoads(pLoads);

    // Compute mass matrix (since it is constant over time)
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.massMatrix();
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> Fext = assembler.rhs();
    loadVec.setZero();
    // surfForce.setValue(loadVec,3);
    assembler.assemble(); // Reset RHS with homogeneous surface force. HV: Bit inefficient, maybe replace with assembleLinearVector in the future.

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) 
    {
        // For the shell
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x, solutions);
        status = assembler.assembleMatrix(solutions);
        m = assembler.matrix();
        return status==ThinShellAssemblerStatus::Success;
    };


    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&dt,&Fext,&assembler,&solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        //Add assemble vector JL
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x,solutions);
        status = assembler.assembleVector(solutions);
        result = assembler.rhs();
        if (t<=dt)
        {
            gsDebug<<"IMPACT!!!\n";
            result += Fext;
        }
        return status==ThinShellAssemblerStatus::Success;
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

    //------------------------------------------------------------------------------
    // Initial Conditions
    //------------------------------------------------------------------------------

    size_t N = assembler.numDofs();
    gsVector<> U(N), V(N), A(N);
    U.setZero();
    V.setZero();
    A.setZero();

    gsFileManager::mkdir(dirname);

    gsParaviewCollection collection(dirname + "/solution");
    real_t time = 0;

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

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
            gsMultiPatch<> solution = assembler.constructDisplacement(U);
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
    GISMO_ERROR("This file requires gsKLShell and gsStructuralAnalysis");
    return EXIT_FAILURE;
#endif
#else
    GISMO_ERROR("This file requires gsKLShell");
    return EXIT_FAILURE;
#endif

}
