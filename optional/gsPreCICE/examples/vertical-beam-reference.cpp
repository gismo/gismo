/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
// #include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsPreCICEVectorFunction.h>

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>
#endif

#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/gsTimeIntegrator.h>
#endif


using namespace gismo;

int main(int argc, char *argv[])
{
#ifdef gsElasticity_ENABLED
#ifdef gsStructuralAnalysis_ENABLED

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    std::string precice_config;
    int method = 2; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe

    index_t Nsteps = 100;

    real_t dt = 1e-2;

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
    neuData<<-1e2,0;
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

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> K_T;

    gsDebugVar(assembler.rhs());

    // Function for the Residual
    gsVector<> F = assembler.rhs();
    std::function<gsMatrix<real_t> (real_t) > Forcing;
    Forcing = [&F](real_t time)
    {
      return math::sin(3.1415923565 * time ) * F;
    };

    gsTimeIntegrator<real_t> timeIntegrator(M,K,Forcing,dt);
    timeIntegrator.verbose();
    timeIntegrator.setTolerance(1e-6);
    timeIntegrator.setMethod(methodName);

    gsDebugVar(M.toDense());
    gsDebugVar(K.toDense());


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

    real_t time = 0;
    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    gsParaviewCollection collection("./solution");
    for (index_t i=0; i<Nsteps; i++)
    {
        time += dt;
        timeIntegrator.step();
        timeIntegrator.constructSolution();
        gsMatrix<> displacements = timeIntegrator.displacements();

        gsDebugVar(displacements.transpose());

        gsMultiPatch<> solution;
        assembler.constructSolution(displacements,fixedDofs,solution);

        // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);
        std::string fileName = "./solution" + util::to_string(i);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(i) + "0";
        collection.addTimestep(fileName,time,".vts");

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
#endif

}
