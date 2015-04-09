// example solving the poisson equation with gismo

#include <iostream>
#include <ctime>

#include <gismo.h>


using namespace gismo;
using std::cout;
using  std::endl;

// Getting the input (defined after the main function)
bool parse_input( int argc, char *argv[], int & numRefine, int & numElevate,  
                  int & Dirichlet, bool & plot, int & plot_pts, 
                  gsMultiPatch<> * & geo, gsPoissonPde<> * & ppde);


int main(int argc, char *argv[])
{

  /////////////////// Input ///////////////////
  int numRefine;  // defaults to 2
  int numElevate; // defaults to -1
  int Dirichlet;  // defaults to 0
  bool plot;      // defaults to false
  int plot_pts;   // defaults to 1000
  gsMultiPatch<> * geo ; // defaults to BSplineCube
  gsPoissonPde<> * ppde ;
  
  bool success = parse_input(argc, argv, numRefine, numElevate, Dirichlet, 
                             plot, plot_pts, geo, ppde);
  if ( ! success )
  {
      cout<<"Input failed, quitting.\n";
      return 1;
  }
  /////////////////// Print info ///////////////////
  cout<<"Type "<< argv[0]<< " -h, to get the list of command line options.\n\n";
  cout<<"Domain: "<< *geo <<"\n";
  cout<<"p-refinent steps before solving: "<< numElevate <<"\n";  
  cout<<"h-refinent steps before solving: "<< numRefine <<"\n";  

  cout<< * ppde <<"\n"; 

  /////////////////// Setup boundary value problem ///////////////////  
  // Boundary conditions
  gsBoundaryConditions<> BCs;
  // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator 
             bit = geo->bBegin(); bit != geo->bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, ppde->solution() );
    }
  
    gsMultiBasis<> bases(*geo);
  
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
  
  cout << "Discret. Space 0: "<< bases[0] << endl;

  
  /////////////////// Setup solver ///////////////////
  gsPoissonAssembler<real_t> poisson(*geo,bases,BCs,*ppde->rhs(),
                                     ( Dirichlet==1 ? dirichlet::nitsche 
                                       : dirichlet::elimination) );

  cout<<"System size: "<< poisson.numDofs() << endl;
  
  // Assemble and solve
  gsStopwatch time;
  poisson.assemble();
  const double assmTime = time.stop();
  Eigen::ConjugateGradient< gsSparseMatrix<> > solver( poisson.matrix() );
  gsMatrix<> solVector = solver.solve( poisson.rhs() );
  const double solveTime = time.stop();
  gsField<>* x = poisson.constructSolution(solVector);

  //cout <<  poisson.linearSystem() <<"\n";

  cout << "Assembling time: " <<  assmTime << " s" << "\n";    
  cout << "Solver time:     " <<  solveTime - assmTime << " s" << "\n";    
  cout << "Total time:      " <<  solveTime << " s" << "\n";    
  cout << "Solution: " << x->function()  << endl;
  if ( ppde->solution() )
      cout << "L2 error: " << x->distanceL2(*ppde->solution() ) 
           << endl;

  // Optionally plot solution in paraview
  if (plot)
  {
      cout<<"Plotting in Paraview..." << endl;
      gsWriteParaview<>( *x, "poisson_sol" , plot_pts) ;
      
      //Plot exact solution in paraview
      gsField<> exact( *geo , *ppde->solution() , false ) ;
      gsWriteParaview<>( exact, "poisson_sol_exact", plot_pts) ;
  }

  delete x;
  delete geo;
  delete ppde;

  return ( plot ? system("paraview poisson_sol.pvd&") : 0 ); 
}


bool parse_input( int argc, char *argv[], int & numRefine, int & numElevate,  
                  int & Dirichlet, bool & plot, int & plot_pts, 
                  gsMultiPatch<> *& geo, gsPoissonPde<> *& ppde)
{
  std::string fn_pde("");
  bool arg_dirich = false;
  plot = false;
  plot_pts = 1000;
  numElevate = -1;
  numRefine = 2;
  std::string fn_basis("");
  std::string fn("");
  
  gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");
  cmd.addString("p", "pde", "File containing a poisson PDE (.xml)", fn_pde);
  cmd.addSwitch("nitsche", "Use the Nitsche's method for Dirichlet sides", arg_dirich);
  cmd.addSwitch("plot", "Plot result in ParaView format", plot);
  cmd.addInt("s", "plotSamples", "Number of sample points to use for plotting", plot_pts);
  cmd.addInt("e", "degreeElevation", "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)",
             numElevate);
  cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving", 
             numRefine);
  cmd.addString("b", "basis", "File containing basis for discretization (.xml)", 
                fn_basis); 
  cmd.addString("g", "geometry", "File containing Geometry (.xml, .axl, .txt)",fn);
  
    
    bool ok = cmd.getValues(argc,argv);
    if (!ok) {
      cout << "Error during parsing!";
      return false;
    }
    
    if ( arg_dirich )
        Dirichlet  = 1;
    else
        Dirichlet  = 0;

    if ( ! fn_basis.empty() )
	{
        cout << "basis reading not implemented.\n"; 
	}

    if (numRefine<0)
    { 
      cout << "Number of refinements must be non-negative, setting to zero.\n"; 
      numRefine = 0;
    }
    if (numElevate<-1)
    { 
      cout << "Number of elevations must be non-negative, ignoring parameter.\n"; 
      numElevate = -1;
    }

    if ( fn_pde.empty() )
	{
	    fn_pde = GISMO_DATA_DIR;
	    if ( !fn.empty() )
		{
		    geo = gsReadFile<>( fn );
		    if ( ! geo )
		      {
			gsWarn<< "Did not find any geometry in "<<fn<<", quitting.\n";
			return false;
		      }
		    switch ( geo->geoDim() )
			{
			case 1:
			    fn_pde+="/pde/poisson1d_sin.xml";
			    break;
			case 2:
			    fn_pde+="/pde/poisson2d_sin.xml";
			    break;
			case 3:
			    fn_pde+="/pde/poisson3d_sin.xml";
			    break;
			default:
			    return false;
			}
		}
	    else
            fn_pde+="/pde/poisson2d_sin.xml";
	}
    ppde = gsReadFile<>(fn_pde);
    if ( !ppde )
      {
	gsWarn<< "Did not find any PDE in "<< fn<<", quitting.\n";
	return false;
      }
    
    if ( fn.empty() )
	{
	  fn = GISMO_DATA_DIR;
      switch ( ppde->m_compat_dim )
	    {
	    case 1:
	      fn+= "domain1d/bspline1d_01.xml";
	      break;
	    case 2:
	      fn+= "domain2d/square.xml";
	      break;
	    case 3:
	      fn+= "volumes/cube.xml";
	      break;
	    default:
	      return false;
	    }
	}
    geo = gsReadFile<>( fn );
    if ( !geo )
      {
	cout << "Did not find any geometries in "<< fn<<", quitting.\n";
	return false;
      }
  return true;
}
