/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

    Give knots/control points for the boundary curve
*/

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

using namespace gismo;

template <short_t DIM, class T, template<short_t _DIM, class _T> class basis_type>
inline void getKnots(const gsBasis<T> & source, std::vector<gsKnotVector<T>> & tensorKnots)
{
    if ( const basis_type<DIM,T> * basis = dynamic_cast<const basis_type<DIM,T>*>(&source) )
        for (index_t d=0; d!=DIM; d++)
            tensorKnots[d] = basis->knots(d);
}


// Helper function check if points are in counterclockwise order

bool isCounterclockwise(const gsMatrix<>& points, const index_t patch_index) 
{
    if (points.rows() < 2) {
        return true;
    }
    // Check if all X coordinates are the same (vertical line case)
    bool allXSame = true;
    for (int i = 1; i < points.rows(); ++i) 
    {
        if (points(i, 0) != points(0, 0)) 
        {
            allXSame = false;
            break;
        }
    }

    if (allXSame && patch_index > 1) 
    {
        gsDebugVar(points(0, 1)>points(1,1));
        return points(0, 1)>points(1,1);
    }

    // Normal case for non-vertical lines
    double sum = 0.0;
    int n = points.rows();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double x1 = points(j, 0) - points(i, 0);
        double y1 = points(j, 1) - points(i, 1);
        double x2 = points((j + 1) % n, 0) - points(j, 0);
        double y2 = points((j + 1) % n, 1) - points(j, 1);
        sum += x1 * y2 - y1 * x2;
    }
    return sum > 0; // Positive sum indicates counterclockwise order
}


// Function to write coupling interface to Paraview (check the boundaries)
void writeCouplingInterfaceParaview(const gsMultiPatch<> &patches, const std::vector<patchSide> &couplingInterfaces, const std::string &filename)
{
    gsMultiPatch<> couplingPatches;

    for (size_t i = 0; i < couplingInterfaces.size(); ++i)
    {
        std::unique_ptr<gsGeometry<> > boundaryGeom(patches.patch(couplingInterfaces[i].patch).boundary(couplingInterfaces[i].side())->clone());
        couplingPatches.addPatch(std::move(boundaryGeom));
    }

    gsWriteParaview<>(couplingPatches, filename);
}


// Helper function, if the points are not in counterclockwise order, reverse them

void ensureCounterclockwise(gsMatrix<>& controlPoints, index_t patch_index) {
    if (!isCounterclockwise(controlPoints, patch_index)) 
    {
        gsDebugVar(controlPoints);
        std::cout << "Reversing points to make counterclockwise for patch index " << patch_index << std::endl;
        for (int i = 0, j = controlPoints.rows() - 1; i < j; ++i, --j) {
            controlPoints.row(i).swap(controlPoints.row(j));
        }
    } else {
        std::cout << "Points are already counterclockwise for patch index " << patch_index << std::endl;
    }
}



int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    std::string precice_config;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("M", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6: RK4",method);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate domain
    // gsMultiPatch<> patches;
    // patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,0.1));
    // 
    // 
    std::string filenameSolid = "pde/TUDFlame_5p.xml";
    gsMultiPatch<real_t> patches;
    gsReadFile<real_t>( filenameSolid, patches);
    

    //! [GetGeometryData]
    gsInfo << "The domain is a "<< patches <<"\n";
    // patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.5,1.0));
    // gsWriteParaview(patches,"output_flame",1000);
    // return 0;


    //Embed the 2D geometry to 3D
    gsMultiPatch<> solutions;
    patches.addAutoBoundaries();
    // p-refine
    if (numElevate!=0)
        patches.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();

    // patches.embed(3);
    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    real_t rho = 12.5;
    real_t E = 400.0;
    real_t nu = 0.3;

    // Set the interface for the precice coupling


    std::vector<patchSide> couplingInterfaces(12);

    couplingInterfaces[0] = patchSide(0,boundary::south);
    couplingInterfaces[1] = patchSide(1,boundary::south);

    couplingInterfaces[2] = patchSide(2,boundary::south);
    couplingInterfaces[3] = patchSide(2,boundary::east);
    couplingInterfaces[4] = patchSide(2,boundary::north);

    couplingInterfaces[5] = patchSide(1,boundary::north);

    couplingInterfaces[6] = patchSide(3,boundary::south);
    couplingInterfaces[7] = patchSide(4,boundary::south);
    couplingInterfaces[8] = patchSide(4,boundary::east);
    couplingInterfaces[9] = patchSide(4,boundary::north);

    couplingInterfaces[10] = patchSide(3,boundary::north);

    couplingInterfaces[11] = patchSide(0,boundary::west);


    // Generate the multi-patch boundary pvd file to check if the boundaries are correct
    // writeCouplingInterfaceParaview(patches, couplingInterfaces, "coupling_interface");

    // return 0;


    /*
     * Initialize the preCICE participant
     *
     *
     */
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);
    gsInfo << "I am participant "<< participantName <<"\n";
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



    // Step 3: initialize the participant
    // real_t precice_dt = participant.initialize();

//----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Coupling side
    // gsFunctionExpr<> g_C("1","0",patches.geoDim());
    gsPreCICEFunction<real_t> g_C(&participant,ForceMesh,ForceData,patches,patches.geoDim(),false);
    // Add all BCs
    // Coupling interface
    // bcInfo.addCondition(0, boundary::north,  condition_type::neumann , &g_C);
    bcInfo.addCondition(0, couplingInterfaces[1],  condition_type::neumann  , &g_C, -1, true);
    bcInfo.addCondition(0, couplingInterfaces[2],  condition_type::neumann  , &g_C, -1, true);
    bcInfo.addCondition(0, couplingInterfaces[3],  condition_type::neumann  , &g_C, -1, true);

    bcInfo.addCondition(0, couplingInterfaces[0],  condition_type::dirichlet  , &g_D, -1, true);
    // bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_C);
    // East side (prescribed temp)
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &g_D, 0);
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &g_D, 1);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

    std::vector<gsGeometry<>::uPtr> boundaries(couplingInterfaces.size());
    // gsMultiPatch<> boundaryGeometry;

    // Derive the basis for boundaries

    // Finish this for both knot vectors and control points
    gsMultiBasis<> boundaryBases;

    size_t vector_size;
    size_t N = 0;

    // Pre-calculate total size N to resize derived_knots just once.
    for(index_t i= 0; i < couplingInterfaces.size();++i)
    {
        boundaries[i] = patches.patch(0).boundary(couplingInterfaces[i].side()); // Add boundary coefficients to a vector
        auto& basis_temp = boundaries[i]->basis(); // Assuming this returns a pointer to gsBasis
        boundaryBases.addBasis(&basis_temp);
        N += knotsToVector(boundaryBases.basis(i)).size();
    }    

    gsVector<> derived_knots;
    derived_knots.resize(N);
    // gsWriteParaview<>(boundaryBases,"boundaryBases",1000);
    // return 0;


    // boundaryBases.basis(i).knots().asMatrix()

    // static_cast<gsBSplineBasis<>>(boundaryBases.basis(i)).knots().asMatrix()

    size_t offset = 0; // Tracks the current fill position in derived_knots.
    for (index_t i = 0; i < couplingInterfaces.size(); ++i) 
    {
        // No need to call knotsToVector twice for the same basis.
        const gsVector<> temp_knots = knotsToVector(boundaryBases.basis(i));
        derived_knots.segment(offset, temp_knots.size()) = temp_knots;
        offset += temp_knots.size(); // Increment offset by the size of the current knot vector.
    }
    gsDebugVar(derived_knots);

    // Control Points
    // In one matrix 
    // gsMatrix<> derived_control_points;

    size_t totalControlPoints = 0;
    size_t controlpoint_offset = 0;
    for (index_t i = 0; i < couplingInterfaces.size(); ++i) 
    {
        totalControlPoints += boundaries[i]->coefs().rows();
    }
    gsMatrix<> derived_control_points(totalControlPoints, boundaries[0]->coefs().cols());
    derived_control_points.setZero();

    for (index_t i = 0; i < couplingInterfaces.size(); ++i) 
    {
        gsMatrix<> temp_control_points = boundaries[i]->coefs();

        if (i!=0 && i<couplingInterfaces.size()-1)
            ensureCounterclockwise(temp_control_points,i);

        for (index_t j = 0; j < temp_control_points.rows(); ++j) 
        {
            derived_control_points.row(controlpoint_offset + j) << temp_control_points.row(j);
        }
        controlpoint_offset += temp_control_points.rows();
    }   

    gsDebugVar(derived_control_points); 

    // Step 1: write the meshes to PreCICE
    // Step 1a: KnotMesh
    // get the knots in a matrix, using the utility function knotsToMatrix
    gsMatrix<> knots = knotsToMatrix(bases.basis(0)); //
    participant.addMesh(KnotMesh,derived_knots.transpose());

    // Step 1b: ControlPointMesh
    // get the control points, in the format where every column is a control point
    gsVector<index_t> controlPointIDs; // needed for writing
    gsMatrix<> controlPoints = derived_control_points;
    // gsDebugVar(controlPoints);

    participant.addMesh(ControlPointMesh,derived_control_points.transpose(),controlPointIDs);

    // Step 1c: ForceMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsAssembler<>::defaultOptions();

    // Get the quadrature points
    gsMatrix<> quadPoints = gsQuadrature::getAllNodes(bases.basis(0),quadOptions,couplingInterfaces);
    gsDebugVar(quadPoints.dim());
    participant.addMesh(ForceMesh,quadPoints);
    gsDebugVar(quadPoints);

    // Step 2 (not needed???)
    // Needed for direct mesh coupling
    gsMatrix<> bbox(2,2);
    bbox.col(0).setConstant(-1e300);
    bbox.col(1).setConstant(1e300);
    bbox.transposeInPlace();
    participant.setMeshAccessRegion(ForceMesh,bbox);
    gsInfo << "Mesh access region set\n";

    real_t precice_dt = participant.initialize();
    gsInfo << "PreCICE initialized\n";



    // for (index_t i = 0; i < couplingInterfaces.size(); ++i) {
    //     derived_control_points.push_back(boundaries[i]->coefs().transpose());
    // }
    // std::cout << "Control Points: " << derived_control_points[0]<< std::endl;


    // gsDebugVar(boundaryBases.basis(0));
    // derived_knots.push_back(knotsToVector(boundaryBases.basis(0)));

    // gsDebugVar(derived_knots);


    // Maybe also make the fluid part as well (for testing).

//----------------------------------------------------------------------------------------------
   // // source function, rhs
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
    assembler.options().setInt("MaterialLaw",material_law::hooke);
    assembler.assemble();
    gsInfo << "Hello, there!\n";


    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> K_T;

    // Time step
    real_t dt = precice_dt;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    // Assemble the RHS
    gsVector<> F = assembler.rhs();
    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(),1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(),1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(),1);

    // Define the solution collection for Paraview
    gsParaviewCollection collection("./solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) {
        // to do: add time dependency of forcing
        assembler.assemble(x, fixedDofs);
        m = assembler.matrix();
        return true;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        assembler.assemble(x,fixedDofs);
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

    real_t time = 0;
    // Plot initial solution
    if (plot)
    {
        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,fixedDofs,solution);

        // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);
        std::string fileName = ".//solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileName,time,".vts");
    }

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

    gsStructuralAnalysisOutput<real_t> writer("pointData.csv",points);
    writer.init({"x","y"},{"time"}); // point1 - x, point1 - y, time

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);

    gsDebugVar(dt);
    gsInfo<<"Starting Simulation\n";

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

        // solve gismo timestep
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time,dt,U,V,A);
        solVector = U;
        gsInfo<<"Finished\n";

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,fixedDofs,solution);
        // write heat fluxes to interface
        // solution.embed(3);
        controlPoints = solution.patch(0).coefs().transpose();

        gsMatrix<> controlPointsOutput(totalControlPoints, boundaries[0]->coefs().cols());
        controlPointsOutput.setZero();
        
        size_t controlpoint_offset = 0;
        for (index_t i = 0; i < couplingInterfaces.size(); ++i) 
        {
            const gsMatrix<> temp_control_points = solution.patch(0).boundary(couplingInterfaces[i].side())->coefs();
            // const gsMatrix<> temp_control_points = patches.patch(0).boundary(couplingInterfaces[i].side())->coefs();
            for (index_t j = 0; j < temp_control_points.rows(); ++j) 
            {
                controlPointsOutput.row(controlpoint_offset + j) << temp_control_points.row(j);
            }
            controlpoint_offset += temp_control_points.rows();
        }  

        gsDebugVar(controlPointsOutput); 

        // gsDebugVar(derived_control_points);
        participant.writeData(ControlPointMesh,ControlPointData,controlPointIDs,controlPointsOutput.transpose());
        gsInfo<<"Wrote data\n";
        
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
                std::string fileName = "./solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = "solution" + util::to_string(timestep) + "0";
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
    