/** @file test-communication-IGA-IGA-fluid.cpp
 * 
 * 	@brief test communicate spline through preCICE. To see if the spline force function can be reconstructed
 * 	
 * 	This file is a part of G+Smo library.
 * 
 * 	This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 * 
 * 	Author(s): J. Li
 * 
 **/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>


#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#endif


using namespace gismo;

// template<typename T>
// T ProjectForceMesh(      const gsMatrix<>       & forceControlPoints, 
//                          const gsMultiBasis<>   & forceBasis,
//                          const gsMultiBasis<>   & geometryBasis,
//                          const index_t          & targetDim,
//                          gsMatrix<T>            & results)
// {
//     gsExprAssembler<> A(1,1);

//     gsMatrix<T> solVector;

//     // Set the discretization space
//     space u = A.getSpace(forceBasis, targetDim);

//     solution sol = A.getSolution(u, solVector);
//     auto f = A.getCoordinates(u);

//     u.setup(0);
//     A.initSystem();



//     gsMatrix<> matrixForce, matrixGeometry;
// }


int main(int argc, char *argv[])
{
    bool plot = false;
    bool get_readTime = false;
    bool get_writeTime = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;

    //! [Parse command line]
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson",method);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");
    //! [Read input file]
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    //Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.5,1.0));

    //Embed the 2D geometry to 3D
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


    /*
     * Data initialization
     *
     * This participant receives mesh information (knot vector, control points) from the solid,
     * and it creates its own mesh based on the information. 
     * And writes pressure, reads displacements from the solid participant.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + GeometryControlPointMesh:        This mesh contains the control points as mesh vertices.
     *   + GeometryKnotMesh:                This mesh contains the knots as mesh vertices.
     *   + ForceKnotMesh:           		This mesh contains the knots as mesh vertices for force mesh.
     *   + ForceControlPointMesh:   		This mesh contains the control points as mesh vertices for force mesh.
     *
     *
     * - Data:
     *   + GeometryControlPointData:        This data is defined on the GeometryControlPointMesh and stores the displacement of the control points
     *	 + GeometryKnotData
     *   + ForceControlPointData:   		This data is defined on the ForceControlPointMesh and stores the change of pressure/stress
     *	 + ForceKnotData
     */

    std::string GeometryKnotMesh               	= "GeometryKnotMesh";
    std::string GeometryControlPointMesh		= "GeometryControlPointMesh";
    std::string ForceKnotMesh           		= "ForceKnotMesh";
    std::string ForceControlPointMesh   		= "ForceControlPointMesh";

    // std::string StressData       = "StressData";
    std::string GeometryControlPointData        = "GeometryControlPointData";
    // std::string GeometryKnotData				= "GeometryKnotData";
    // std::string ForceKnotData           		= "ForceKnotMeshData";
    std::string ForceControlPointData   		= "ForceControlPointData";

    gsMatrix<> bbox(3,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(ForceControlPointMesh, bbox);

    // Setup geometry control point mesh
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints = patches.patch(0).coefs().transpose();
    participant.addMesh(GeometryControlPointMesh, geometryControlPoints, geometryControlPointIDs);


    // Setup the geometry knot mesh
    gsMatrix<> geometryKnots = knotsToMatrix(bases.basis(0)); 


    participant.addMesh(GeometryKnotMesh,geometryKnots);

    real_t precice_dt = participant.initialize();

    //G[et the force mesh from the coupling interface
    gsVector<index_t> forceKnotIDs;
    gsMatrix<> forceKnots;
    participant.getMeshVertexIDsAndCoordinates(ForceKnotMesh,forceKnotIDs,forceKnots);

    gsBasis<> * basis = knotMatrixToBasis<real_t>(forceKnots).get();

    gsVector<index_t> forceControlPointIDs;
    gsMatrix<> forceControlPoints;
    participant.getMeshVertexIDsAndCoordinates(ForceControlPointMesh, forceControlPointIDs,forceControlPoints);


    // // Step 2: Regenerate the geometry
    gsMultiPatch<> forceMesh; //Geometry object belongs to gsFunctionSet
    forceMesh.addPatch(give(basis->makeGeometry(forceControlPoints.transpose())));


    // Define boundary condition for solid mesh
    gsBoundaryConditions<> bcInfo;
    //Bottom Side
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, -1);
    bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr, 2);

    // Assign geometry map
    bcInfo.setGeoMap(patches);

    // Surface force equals to the force mesh

//------------------------------------------ gsKLShell Module ------------------------------------------------
    //Setup the material properties
    gsFunctionExpr<> E_modulus(std::to_string(E),3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu),3);
    gsFunctionExpr<> Density(std::to_string(rho),3);

    gsOptionList options;


    //Define thickness
    real_t thickness = 0.1;

    gsFunctionExpr<> t(std::to_string(thickness),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsMaterialMatrixBase<real_t>* materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3, real_t>(patches, t, parameters, Density, options);

    // REMOVE IT LATER
    
    // forceMesh.embed(3);

    gsThinShellAssembler<3, real_t, true> assembler(patches, bases, bcInfo, forceMesh, materialMatrix);
    gsOptionList assemblerOptions = options.wrapIntoGroup("Assembler");

    assembler.assemble();
    // forceMesh.patch(0).coefs().setZero();
    // assembler.assemble();
    // gsDebugVar(assembler.rhs());

    assembler.setOptions(assemblerOptions);

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;


    // Compute the mass matrix (since it is constant over time)
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.massMatrix();
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();


    // Define the solution collection for Paraview
    std::string dirname = "./output";

    gsFileManager::mkdir(dirname); 
    gsParaviewCollection collection(dirname + "/solution");

    // Time step
    real_t dt = precice_dt;
    real_t t_read = 0;
    real_t t_write = 0;

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
        solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
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

    // gsDebugVar(dt);

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
        participant.readData(ForceControlPointMesh,ForceControlPointData,forceControlPointIDs,forceControlPoints);

        if (get_readTime)
            t_read += participant.readTime();
        forceMesh.patch(0).coefs() = forceControlPoints.transpose();
        // forceMesh.embed(3);
        assembler.assemble();
        F = assembler.rhs();
        // gsDebugVar(F);

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
        // write heat fluxes to interface
        geometryControlPoints = solution.patch(0).coefs().transpose();
        participant.writeData(GeometryControlPointMesh,GeometryControlPointData,geometryControlPointIDs,geometryControlPoints);
        
        if (get_writeTime)
            t_write +=participant.writeTime();


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
            // gsTimeIntegrator advances the time step                   
            // advance variables
            time += dt;
            timestep++;

            gsField<> solField(patches,solution);
            if (timestep % plotmod==0 && plot)
            {
                // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                std::string fileName = dirname + "/solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName,time,".vts");
            }
            solution.patch(0).eval_into(points,pointDataMatrix);
            otherDataMatrix<<time;
            writer.add(pointDataMatrix,otherDataMatrix);
        }
    }
    if (get_readTime)
    {
        gsInfo << "Solid Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Solid Write time: " << t_write << "\n";
    }

    if (plot)
    {
        collection.save();
    }

    delete timeIntegrator;
    return  EXIT_SUCCESS;
}
