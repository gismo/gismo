/** @file poisson.cpp

    @brief Reads the data for a Poisson boundary value problem
     (coefficients, domain, boundary conditions) from an XML file and
     solves the BVP

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

bool read_input( int argc, char *argv[], int & numRefine, int & numElevate,  
                 int & Dirichlet, bool & plot, int & plot_pts, 
                 gsMultiPatch<> & geo, gsPoissonPde<> * & ppde,
                 gsBoundaryConditions<> & BCs);

int main(int argc, char *argv[])
{

    /////////////////// Input ///////////////////
    int numRefine;  // defaults to 2
    int numElevate; // defaults to -1
    int Dirichlet;  // defaults to 0
    bool plot;      // defaults to false
    int plot_pts;   // defaults to 1000
    gsMultiPatch<> geo ; // defaults to BSplineCube
    gsPoissonPde<> * ppde ;
    gsBoundaryConditions<> BCs;

    bool success = read_input(argc, argv, numRefine, numElevate, Dirichlet, 
                              plot, plot_pts, geo, ppde, BCs);
    if ( ! success )
    {
        gsWarn <<"Input reading failed, quitting.\n";
        return 1;
    }
    /////////////////// Print info ///////////////////
    gsInfo<<"Type "<< argv[0]<< " -h, to get the list of command line options.\n\n";
    gsInfo<<"Domain: "<< geo <<"\n";
    gsInfo<<"p-refinent steps before solving: "<< numElevate <<"\n";  
    gsInfo<<"h-refinent steps before solving: "<< numRefine <<"\n";  

    gsInfo<< * ppde <<"\n"; 

    gsMultiBasis<> bases(geo);
  
    // Elevate and refine the solution space
    if ( numElevate > -1 )
    {
        // get maximun degree
        int tmp = bases.maxDegree(0);
      
        // Elevate all degrees uniformly
        tmp += numElevate;
        for (size_t j = 0; j < bases.nBases(); ++j )
            bases[j].setDegree(tmp);
    }

    for (size_t j = 0; j < bases.nBases(); ++j )
        for (int i = 0; i < numRefine; ++i)
            bases[j].uniformRefine();
  
    gsInfo << "Discret. Space 0: "<< bases[0] << "\n";
    gsInfo << BCs<<"\n";

    double time, totalTime = 0.0;
    gsStopwatch watch;
  
    gsInfo << "Setup problem.. \n";
    gsPoissonAssembler<real_t> poisson(geo, bases, BCs, *ppde->rhs(),
                                       ( Dirichlet==1 ? dirichlet::nitsche 
                                         : dirichlet::elimination) );
    time = watch.stop();
    totalTime += time;
    gsInfo << "time: " <<  time << " s" << "\n";
    gsInfo<<"System size: "<< poisson.numDofs() << "\n";
  
    // Assemble and solve

    gsInfo << "Assembling.. \n";
    watch.restart();
    poisson.assemble();
    time = watch.stop();
    totalTime += time;
    gsInfo << "time: " <<  time << " s" << "\n";

    gsInfo << "Solving.. \n";
    //gsDebugVar(poisson.matrix().isCompressed());
    watch.restart();
    gsSparseSolver<>::CGDiagonal solver( poisson.matrix() );
    gsMatrix<> solVector = solver.solve( poisson.rhs() );
    time = watch.stop();
    totalTime += time;
    gsInfo << "time: " <<  time << " s" << "\n";    

    gsInfo << "Constructing solution.. \n";
    watch.restart();
    gsField<>* x = poisson.constructSolution(solVector);
    time = watch.stop();
    totalTime += time;
    gsInfo << "time: " <<  time << " s" << "\n";    

    //gsInfo << "Solution: " << x->function()  << "\n";

    if ( ppde->solution() )
    {
        gsInfo << "Computing L2 error.. \n";
        watch.restart();
        const real_t L2error = x->distanceL2(*ppde->solution() );
        gsInfo << "time: " <<  time << " s" << "\n";
        totalTime += time;
        gsInfo << "L2 error: " << L2error<< "\n";
    }

    gsInfo << "Total time: " <<  totalTime << " s" << "\n";

    // Optionally plot solution in paraview
    if (plot)
    {
        gsInfo<<"Plotting in Paraview..." << "\n";
        gsWriteParaview<>( *x, "poisson_sol" , plot_pts) ;
      
        //Plot exact solution in paraview
        gsField<> exact( geo , *ppde->solution() , false ) ;
        gsWriteParaview<>( exact, "poisson_sol_exact", plot_pts) ;
    }

    // bool save = false;
    // if (save)

    delete x;
    delete ppde;

    return ( plot ? system("paraview poisson_sol.pvd&") : 0 ); 
}

bool read_input( int argc, char *argv[], 
                 int & numRefine, // opts
                 int & numElevate, // opts  
                 int & Dirichlet, // opts
                 bool & plot, int & plot_pts, //args
                 gsMultiPatch<> & geo, 
                 gsPoissonPde<> *& ppde,
                 gsBoundaryConditions<> & BCs)
{
    bool arg_dirich = false;
    plot = false;
    plot_pts = 1000;
    numElevate = -1;
    numRefine = 2;
    std::string fn( GISMO_DATA_DIR "/planar/two_squares.xml");
  
    gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");
    cmd.addSwitch("nitsche", "Use the Nitsche's method for Dirichlet sides", arg_dirich);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("s", "plotSamples", "Number of sample points to use for plotting", plot_pts);
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)",
               numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", 
               numRefine);
    cmd.addString("g", "geometry", "File containing Geometry (.xml, .axl, .txt)", fn);
      
    bool ok = cmd.getValues(argc,argv);
    if (!ok) 
    {
        gsWarn << "Error during parsing!";
        return 0;
    }
    
    if ( arg_dirich )
        Dirichlet  = 1;
    else
        Dirichlet  = 0;

    if (numRefine<0)
    { 
        gsWarn << "Number of refinements must be non-negative, setting to zero.\n"; 
        numRefine = 0;
    }

    if (numElevate<-1)
    { 
        gsWarn << "Number of elevations must be non-negative, ignoring parameter.\n"; 
        numElevate = -1;
    }

    // Load all data from the XML file
    gsFileData<> xmlData(fn);

    // Check if all needed data exist in the file
    if ( !xmlData.has< gsMultiPatch<> >() )
    {
        gsWarn <<"Did not find a multipatch domain in this file.\n";  
        return false;
    }

    if ( !xmlData.has<gsPoissonPde<> >() )
    {
        gsWarn <<"Did not find a poisson PDE domain in this file.\n"; 
        return false;
    }

    if ( !xmlData.has<gsBoundaryConditions<> >() )
    {
        gsWarn <<"Did not find boundary conditions in this this file.\n"; 
        // default to homogeneous dirichlet ?
        return false;
    }

    ppde = xmlData.getFirst<gsPoissonPde<> >();
    xmlData.getFirst(geo);
    xmlData.getFirst(BCs);

    return true;
}
