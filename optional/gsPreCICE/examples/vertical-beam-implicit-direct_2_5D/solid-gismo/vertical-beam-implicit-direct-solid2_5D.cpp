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
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>


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
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    std::string precice_config;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    std::string dirname = "./output";


    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    //cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson",method);
    // cmd.addSwitch("write", "Create a file with point data", write);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.5,1.0));

    // Embed the 2D geometry to 3D
    gsMultiPatch<> solutions;
    patches.addAutoBoundaries();
    patches.embed(3);


    // Create bases
    // p-refine
    if (numElevate!=0)
        patches.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    // real_t nu = 0.3;
    // real_t mu = E / (2.0 * (1.0 + nu));
    // real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    /*
     * Initialize the preCICE participant
     *
     *
     */
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    /*
     * Data initialization
     *
     * This participant manages the geometry. The follow meshes and data are made available:
     *
     * - Meshes:
     *   + KnotMesh             This mesh contains the knots as mesh vertices
     *   + ControlPointMesh:    This mesh contains the control points as mesh vertices
     *   + ForceMesh:           This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + ControlPointData:    This data is defined on the ControlPointMesh and stores the displacement of the control points
     *   + ForceData:           This data is defined on the ForceMesh and stores pressure/forces
     */
    std::string KnotMesh        = "KnotMesh";
    std::string ControlPointMesh= "ControlPointMesh";
    std::string ForceMesh       = "ForceMesh";

    std::string ControlPointData= "ControlPointData";
    std::string ForceData       = "ForceData";

    // Step 1: write the meshes to PreCICE
    // Step 1a: KnotMesh
    // get the knots in a matrix, using the utility function knotsToMatrix
    gsMatrix<> knots = knotsToMatrix(bases.basis(0)); //
    participant.addMesh(KnotMesh,knots);

    // Step 1b: ControlPointMesh
    // get the control points, in the format where every column is a control point
    gsVector<index_t> controlPointIDs; // needed for writing
    gsMatrix<> controlPoints = patches.patch(0).coefs().transpose();
    participant.addMesh(ControlPointMesh,controlPoints,controlPointIDs);

    // Step 1c: ForceMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();

    // Get the quadrature points
    gsMatrix<> quadPoints = gsQuadrature::getAllNodes(bases.basis(0),quadOptions);
    participant.addMesh(ForceMesh,quadPoints);

    // Step 2 (not needed???)
    // Needed for direct mesh coupling
    // gsMatrix<> bbox = bases.basis(0).support();
    // bbox.transposeInPlace();
    // bbox.resize(1,bbox.rows()*bbox.cols());
    // participant.setMeshAccessRegion(ForceMesh,bbox);

    // gsMatrix<> bbox(3,2);
    // bbox << 0,1,
    //         0,1,
    //         0,1;
    // bbox.transposeInPlace();
    // bbox.resize(1,bbox.rows()*bbox.cols());
    // participant.setMeshAccessRegion(ForceMesh,bbox);
    // Step 3: initialize the participant
    real_t precice_dt = participant.initialize();

//----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;

    // Bottom side
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, -1);
    bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr, 2);

    // Assign geometry map
    bcInfo.setGeoMap(patches);

    // Surface force
    gsPreCICEFunction<real_t> surfForce(&participant,ForceMesh,ForceData,patches,patches.parDim(),patches.geoDim(),false);


//------------------------------------------IF gsElasticity is Used------------------------------
    // // creating mass assembler
    // gsMassAssembler<real_t> massAssembler(patches,bases,bcInfo,gravity);
    // massAssembler.options().setReal("Density",rho);
    // massAssembler.assemble();

    // // creating stiffness assembler.
    // gsElasticityAssembler<real_t> assembler(patches,bases,bcInfo,g);
    // assembler.options().setReal("YoungsModulus",E);
    // assembler.options().setReal("PoissonsRatio",nu);
    // assembler.options().setInt("MaterialLaw",material_law::hooke);
    // assembler.assemble();


    // gsMatrix<real_t> Minv;
    // gsSparseMatrix<> M = massAssembler.matrix();
    // gsSparseMatrix<> K = assembler.matrix();
    // gsSparseMatrix<> K_T;
//------------------------------------------IF gsElasticity is Used------------------------------


    // Set up the material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E),3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu),3);
    gsFunctionExpr<> Density(std::to_string(rho),3);

    //Define thickness
    real_t thickness = 0.1;

    gsFunctionExpr<> t(std::to_string(thickness),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;

    gsMaterialMatrixBase<real_t>* materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3, real_t>(patches, t, parameters, Density, options);

    gsThinShellAssembler<3, real_t, true> assembler(patches, bases, bcInfo, surfForce, materialMatrix);
    gsOptionList assemblerOptions = quadOptions.wrapIntoGroup("ExprAssembler");
    assembler.setOptions(assemblerOptions);

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;

    // Compute the mass matrix (since it is constant over time)
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.massMatrix();
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();

    
    // Define the solution collection for Paraview
    gsFileManager::mkdir(dirname); 
    gsParaviewCollection collection(dirname + "/solution");

        // Time step
    real_t dt = precice_dt;

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) 
    {
        // to do: add time dependency of forcing
        // For the shell
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x, solutions);
        status = assembler.assembleMatrix(solutions);
        m = assembler.matrix();
        return true;
    };


    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        //Add assemble vector JL
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x,solutions);
        status = assembler.assembleVector(solutions);
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



    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    // Assemble the RHS
    gsVector<> F = assembler.rhs();

    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(),1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(),1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(),1);


    real_t time = 0;
    // Plot initial solution
    if (plot)
    {
        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements);

        // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);

        std::string fileName = dirname + "/solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileName,time,".vts");
    }

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

    gsStructuralAnalysisOutput<real_t> writer("./output/pointData.csv",points);
    writer.init({"x","y","z"},{"time"}); // point1 - x, point1 - y, point1 - z, time

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);

    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;

            gsInfo<<"Checkpoint written:\n";
            gsInfo<<"\t ||U|| = "<<U.norm()<<"\n";
            gsInfo<<"\t ||V|| = "<<V.norm()<<"\n";
            gsInfo<<"\t ||A|| = "<<A.norm()<<"\n";


            timestep_checkpoint = timestep;
        }

        assembler.assemble();
        F = assembler.rhs();
        gsDebugVar(F);

        // solve gismo timestep
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time,dt,U,V,A);
        solVector = U;

        gsInfo<<"Finished\n";

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements);
        controlPoints = solution.patch(0).coefs().transpose();
        participant.writeData(ControlPointMesh,ControlPointData,controlPointIDs,controlPoints);

        // do the coupling
        precice_dt =participant.advance(dt);

        if (participant.requiresReadingCheckpoint())
        {
            U = U_checkpoint;
            V = V_checkpoint;
            A = A_checkpoint;
            timestep = timestep_checkpoint;
        }
        else
        {
            // gsTimeIntegrator advances
            // advance variables
            time += dt;
            timestep++;

            gsField<> solField(patches,solution);
            if (timestep % plotmod==0 && plot)
            {
                // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                std::string fileName = dirname + "/solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = dirname + "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName,time,".vts");
            }
            solution.patch(0).eval_into(points,pointDataMatrix);
            otherDataMatrix<<time;
            writer.add(pointDataMatrix,otherDataMatrix);
        }
    }

    if (plot)
    {
        collection.save();
    }

    delete timeIntegrator;
    return  EXIT_SUCCESS;
}
