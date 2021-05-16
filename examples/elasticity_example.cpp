/** @file elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
// #include <gsElasticity/gsLinearMaterial.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    index_t testCase = 1;
    bool last = false;
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("????.");
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine-1; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setOptions(Aopt);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    variable ff = A.getCoeff(f, G);

    // Recover manufactured solution
    variable u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);


    /// --------------------------------------------------

    //gsVector<> pt(2);
    //pt.setConstant(0.5);

    // gsLinearMaterial<> matf(1,0.3);
    // variable mat = A.getCoeff(matf, G);

    /* plane strain
                  		  1-nu           , nu	     , 0
       E/(1+nu)(1-2nu)		* nu  		 , 1-nu      , 0
                  		  0          	 , 0         , (1-2nu)/2
       
       plane stress
                  		  1 	         , nu	     , 0
       E/(1-nu^2)		* nu  		 , 1         , 0
                  		  0          	 , 0         , (1-nu)/2
     */

    bool bl_plane_stress = True;

    real_t E = 1;
    real_t nu = 0.3;

    if (bl_plane_stress)
    {
       real_t mm_factor = E/(1-pow(nu,2));
       gsMatrix<> C(3,3);
       C.row(0)<<mm_factor   , mm_factor*nu  ,0;
       C.row(1)<<mm_factor*nu, mm_factor     ,0;
       C.row(2)<<0           , 0             ,mm_factor*((1-nu)/2);
       C.resize(9,1);
    }


    gsConstantFunction<> Cfun(C,2);
    variable mat = A.getCoeff(Cfun, G);

    //gsInfo<<ev.eval(reshape(mat,3,3),pt)<<"\n";

    //space v = A.getSpace(dbasis,2);


    //gsInfo<<ev.eval(v,pt)<<"\n";


    //gsInfo<<ev.eval(jac(v).tr()*jac(v),pt)<<"\n";

    /* assemble proposal
       assemble weak form: int( w K u )dOmega = int( w f )dOmega
     			   K = B^T D B

       
 				N_1,x , 0 
       2D Case:		B_i = 	0, N_1,y	, with i being shape function index
    				N_1,y , N_1,x


       assembly of B-matrix for multiple trial functions
	
		N_1,x , 0 	N_2,x , 0
	B = 	0 , N_1,y	0 , N_2,y	...
    		N_1,y , N_1,x 	N_2,y , N_2,x

       this should probably go into an _expr as in kirchhoff-love example.

       possible to get derivatives from jac(u) within that expression and then just rearrange them accordingly? 
	
    */ 

    // A.assemble(setB(u,G).transpose() * mat * setB(u,G) * u * meas(G), u * f * meas(G));

    /// --------------------------------------------------

    return 0;

    //! [Problem setup]

    gsSparseSolver<>::CGDiagonal solver;

    u.setup(bc, dirichlet::interpolation, 0);

    // Initialize the system
    A.initSystem(false);

    gsInfo<< A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    variable g_N = A.getBdrFunction();
    A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
