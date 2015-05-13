
#include <iostream>

#include <gsAssembler/gsPoissonAssembler2.h>

#include <gismo.h>

using namespace std;
using namespace gismo;


int main(int argc, char *argv[])
{   
    // Input options
    int numElevate  = 0;
    int numHref     = 1;
    int basisDegree = 0;
    bool plot       = false;

    // Multipatch object
    gsMultiPatch<> mp;

    // Pde
    gsPoissonPde<> * pde = NULL;

    int result = 0;
    std::string fn(GISMO_DDATA_DIR "/planar/two_squares.xml");
    std::string fn_pde("");
    
    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", 
                   "Number of dyadic h-refinement (bisection) steps to perform before solving",
                   numHref);
    cmd.addInt("p","degree",
                   "Degree of the basis functions to use for solving (will elevate or reduce the input)",
                   basisDegree);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",
                      fn);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    bool ok = cmd.getValues(argc,argv);
    if ( !ok ) 
    {
        std::cout << "Something went wrong when reading the command line. Exiting.\n";
        return 1;
    }
        
    try // Read input
    {
        gsMultiPatch<> * mpptr = gsReadFile<>(fn);
        mp = * safe(mpptr);
        // Read PDE (and known exact solution)
        fn = fn_pde;

    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; return -1; }
    
    if ( fn.empty() )
        switch ( mp.geoDim() )
        {
        case 1:
            fn = GISMO_DDATA_DIR "/pde/poisson1d_sin.xml" ;
            break;
        case 2:
            fn = GISMO_DDATA_DIR "/pde/poisson2d_sin.xml";
            break;
        case 3:
            fn = GISMO_DDATA_DIR "/pde/poisson3d_sin.xml";
            break;
        default:
            cout <<"Got "<< mp;
        }
    pde = gsReadFile<>(fn) ;

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator 
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, pde->solution() );
    }
    
    // Set up and elevate discretization bases
    gsMultiBasis<> bases(mp);

    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);
//               bases.degreeReduce(1);	

    // Run tests
    gsMatrix<> testsL2(numHref+1,7);
    testsL2.setZero();

    gsMatrix<> testsH1(numHref+1,7);
    testsH1.setZero();
    
    //Eigen::ConjugateGradient< gsSparseMatrix<> > solver; // identity preconditioner
    Eigen::ConjugateGradient< gsSparseMatrix<>, Eigen::Lower, Eigen::DiagonalPreconditioner<real_t> > solver;
    //Eigen::ConjugateGradient< gsSparseMatrix<>, Eigen::Lower, Eigen::IncompleteLUT<real_t> > solver;
    //Eigen::BiCGSTAB< gsSparseMatrix<>, Eigen::IncompleteLUT<real_t> > solver( galerkin.matrix() );
    //Eigen::BiCGSTAB< gsSparseMatrix<> > solver( galerkin.fullMatrix() );

    int i = 0;
    do
    {
        // Setup the assemblers
        gsPoissonAssembler2<real_t> galerkin(mp,bases,BCs,*pde->rhs(),
                                            dirichlet::elimination,iFace::glue);
        gsPoissonAssembler2<real_t> galerkin_wBC(mp,bases,BCs,*pde->rhs(),
                                                dirichlet::nitsche,iFace::glue);
        gsPoissonAssembler2<real_t> galerkin_dg(mp,bases,BCs,*pde->rhs(),
                                               dirichlet::nitsche,iFace::dg);
        //gsPoissonAssemler<> gp( bvp, bases, elimination, dg  );

        cout<<"Discretization Space for patch 0: \n"<< bases[0] << endl;
        cout<<"Initial DoFs             : "<< galerkin_dg.numDofs() << endl;  
        cout<<"Penalty constant         : "<< galerkin_dg.penalty(0) << endl;  
        cout<<"Gauss nodes per direction: "<< endl;  
        cout << "---------------------------------------\n";
        cout<<"System  size (elim. BCs)  : "<< galerkin.numDofs() << endl;
        cout<<"Coupled size (elim. BCs)  : "<< galerkin.dofMapper().coupledSize() << endl;  
        cout<<"System size (Nitsche BCs) : "<< galerkin_wBC.numDofs() << endl;  
        cout<<"Coupled size(Nitsche BCs) : "<< galerkin_wBC.dofMapper().coupledSize() << endl;  
        cout<<"System size (dg)          : "<< galerkin_dg.numDofs() << endl;
        cout << "---------------------------------------\n";
        
        cout<< "Computing conforming C^0 solution..\n";
        gsStopwatch time;
        galerkin.assemble();
        const double assTime = time.stop();
        cout << "Assembly time (elim. BCs): " << assTime << " s" << "\n";
        time.restart();
        
        gsMatrix<> solVector;

        double solvTime(0.0);
        if ( galerkin.numDofs() )
        {
            solver.compute( galerkin.matrix() );
            solVector = solver.solve( galerkin.rhs() );
            solvTime = time.stop();
            cout << "Solving time (elim. BCs): " << solvTime << " s" << "\n";
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        }
            gsField<> * sol = galerkin.constructSolution(solVector);
            
        cout << "Computing solution with weakly imposed BCs..\n";
        time.restart();
        galerkin_wBC.assemble();
        const double assTime_wBC = time.stop();
        cout << "Assembly time (weak  BCs): " << assTime_wBC << " s" << "\n";
        time.restart();
        solver.compute( galerkin_wBC.matrix() );
        solVector = solver.solve( galerkin_wBC.rhs() );
        const double solvTime_wBC = time.stop();
        cout << "Solving time (weak  BCs): " << solvTime_wBC << " s" << "\n";
        gsInfo << "residual error: " << solver.error() << "\n";
        gsInfo << "    iterations: " << solver.iterations() << "\n";
        gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        gsField<> * sol_wBC = galerkin_wBC.constructSolution(solVector);
        
        gsField<> * sol_dg = 0;
        double solvTime_dg(0), assTime_dg(0);
        if ( mp.nPatches() > 1 )
        {
            cout << "Computing solution with patch-wise Disc. Galerkin method..\n";
            time.restart();
            galerkin_dg.assemble();
            assTime_dg = time.stop();
            cout << "Assembly time (full DG): " << assTime_dg << " s" << "\n";
            time.restart();
            solver.compute( galerkin_dg.matrix() );
            solVector = solver.solve( galerkin_dg.rhs() );
            solvTime_dg = time.stop();
            cout << "Solving time (full DG): " << solvTime_dg << " s" << "\n";
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
            sol_dg = galerkin_dg.constructSolution(solVector);
        }

        // Collect data
        testsL2(i,0)  =
            testsH1(i,0)= bases.totalSize();
        if ( pde->solutionGiven() )
        {
            testsL2(i,1)= sol->distanceL2    ( *pde->solution() ) ;
            testsL2(i,3)= sol_wBC->distanceL2( *pde->solution() ) ;
            
            testsH1(i,1) = sol->distanceH1    ( *pde->solution() ) ;
            testsH1(i,3) = sol_wBC->distanceH1( *pde->solution() ) ;
        }
                                          
        if ( mp.nPatches() > 1 )
        {
            testsL2(i,5)= sol_dg->distanceL2( *pde->solution() ) ;
            testsH1(i,5)= igaFieldDGDistance(*sol_dg, *pde->solution(),false);
        }

        if (i > 0)
        {
            testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);
            testsL2(i,4)= testsL2(i-1,3) / testsL2(i,3);
            testsL2(i,6)= testsL2(i-1,5) / testsL2(i,5);

            testsH1(i,2)= testsH1(i-1,1) / testsH1(i,1);
            testsH1(i,4)= testsH1(i-1,3) / testsH1(i,3);
            testsH1(i,6)= testsH1(i-1,5) / testsH1(i,5);
        }
        
        cout << "Total time (elim. BCs): " << assTime+solvTime   << " s" << "\n";
        cout << "Total time (weak  BCs): " << assTime_wBC+solvTime_wBC << " s" << "\n";
        if ( mp.nPatches() > 1 )
            cout << "Total time (full DG)  : " << assTime_dg+solvTime_dg  << " s" << "\n";    
        cout << "---------------------------------------\n";    
        cout << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
        cout << "    Dofs   |  L2 error  | err. ratio|  L2 error  | err. ratio|  L2 error  | err. ratio    \n" << testsL2.row(i)  << endl;  
        cout << "           |  H1 error  | err. ratio|  H1 error  | err. ratio|  H1 error  | err. ratio    \n" << testsH1.row(i)  << endl;  

        if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            std::cout<<"Plotting in Paraview...\n";
            gsWriteParaview<>( *sol, "poisson_problem", 1000);
            
            // Run paraview
                result = system("paraview poisson_problem.pvd &");
        }
        
        delete sol;
        delete sol_wBC;
        if ( sol_dg )
            delete sol_dg;
        
        bases.uniformRefine();
    } 
    while ( i++ < numHref );
    
    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= std::log(testsL2(i,2))/std::log(2.0);
        testsL2(i,4)= std::log(testsL2(i,4))/std::log(2.0);
        testsL2(i,6)= std::log(testsL2(i,6))/std::log(2.0);

        testsH1(i,2)= std::log(testsH1(i,2))/std::log(2.0);
        testsH1(i,4)= std::log(testsH1(i,4))/std::log(2.0);
        testsH1(i,6)= std::log(testsH1(i,6))/std::log(2.0);
    }
    
    if ( mp.nPatches() == 1 )
    {
        testsL2.col(6).setZero();
        testsH1.col(6).setZero();
    }
    
    cout << "Summary:\n\n";
    cout << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    cout << "    Dofs   |  L2 error  | conv. rate|  L2 error  | conv. rate|  L2 error  | conv. rate    \n" << testsL2  << endl;


    cout << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    cout << "    Dofs   |  H1 error  | conv. rate|  H1 error  | conv. rate|  H1 error  | conv. rate    \n" << testsH1  << endl;


    delete pde;

    return result;
}
