// example solving the poisson equation with gismo

#include <iostream>
#include <ctime>

#include <gismo.h>



using namespace gismo;
using std::cout;
using  std::endl;

// Getting the input (defined after the main function)
bool parse_input( int argc, char *argv[], int & numRefine, int & numElevate,  
                  int & Dirichlet, int & DG, bool & plot, int & plot_pts, 
                  gsMultiPatch<> * & geo, gsPoissonPde<> * & ppde, 
                  std::vector< gsBasis<>* > & bases );


int main(int argc, char *argv[])
{

  /////////////////// Input ///////////////////
  int numRefine;  // defaults to 2
  int numElevate; // defaults to -1
  int Dirichlet;  // defaults to 0
  int DG;         // defaults to 0
  bool plot;      // defaults to false
  int plot_pts;   // defaults to 1000
  gsMultiPatch<> * geo ; // defaults to BSplineCube
  gsPoissonPde<> * ppde ;
  std::vector< gsBasis<>* > bases;// not yet given by input
  
  bool success = parse_input(argc, argv, numRefine, numElevate, Dirichlet, 
                             DG, plot, plot_pts, geo, ppde, bases);
  if ( ! success )
  {
      cout<<"Input failed, quitting.\n";
      return 0;
  }
  /////////////////// Print info ///////////////////
  cout<<"Type "<< argv[0]<< " -h, to get the list of command line options.\n\n";
  cout<<"Domain: "<< *geo <<"\n";
  cout<<"p-refinent steps before solving: "<< numElevate <<"\n";  
  cout<<"h-refinent steps before solving: "<< numRefine <<"\n";  

  cout<< * ppde <<"\n"; 

  /////////////////// Setup boundary value problem ///////////////////  
  gsBVProblem<> bvp( *geo, ppde ) ;
  // Create Dirichlet boundary conditions for all boundaries
  for (gsMultiPatch<>::const_biterator 
       bit = geo->bBegin(); bit != geo->bEnd(); ++bit)
  {
      bvp.addCondition( *bit, boundary::dirichlet, ppde->solution() );
  }

  cout<< bvp <<"\n"; 
  
  if ( bases.size() == 0 )
      bases = geo->basesCopy();
  
  // Elevate and refine the solution space
  if ( numElevate > -1 )
  {
      // get maximun degree
      int tmp = bases[0]->degree(0);
      for (int j = 0; j < geo->parDim(); ++j )
          if ( tmp < bases[0]->degree(j) )
              tmp = bases[0]->degree(j);
      
      // Elevate all degrees uniformly
      tmp += numElevate;
      for (size_t j = 0; j < bases.size(); ++j )
              bases[j]->setDegree(tmp);
  }

  for (size_t j = 0; j < bases.size(); ++j )
      for (int i = 0; i < numRefine; ++i)
          bases[j]->uniformRefine();
  
  cout << "Discret. Space 0: "<< *bases[0] << endl;

  
  /////////////////// Setup solver ///////////////////
  gsGalerkinMethod<> poisson( bvp, bases);
  if ( Dirichlet == 1)
  { 
      cout<<"Using Nitsche's method for Dirichlet boundaries.\n";  
      poisson.setDirichletStrategy(dirichlet::nitsche);
  }

  if ( DG == 1)
  { 
      cout<<"Using DG method for patch interfaces.\n";  
      poisson.setInterfaceStrategy(iFace::dg);
  }

  cout<<"Dofs: "<< poisson.dofs() << "\n";  
  cout<<"System size: "<< poisson.freeDofs() << endl;
  
  // Assemble and solve
  gsStopwatch time;
  poisson.assemble();
  const double assmTime = time.stop();
  gsField<>* x = poisson.solve();
  const double solveTime = time.stop();

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
      
      char cmd[100];
      if (geo->nPatches() == 1)
	  strcpy(cmd,"paraview poisson_sol.vts\0");
      else
	  strcpy(cmd,"paraview poisson_sol.pvd\0");
      strcat(cmd," &");

      delete x;
      delete geo;
      freeAll(bases);
      return system(cmd); 
  }

  delete x;
  delete geo;
  freeAll(bases);

  return 0; 
}


bool parse_input( int argc, char *argv[], int & numRefine, int & numElevate,  
                  int & Dirichlet, int & DG, bool & plot, int & plot_pts, 
                  gsMultiPatch<> *& geo, gsPoissonPde<> *& ppde, 
                  std::vector< gsBasis<>* > &  bases )
{
  try 
  {
    gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");    

    gsArgVal<std::string> a6("p","pde","File containing a poisson PDE (.xml)", 
			     false,"", "PDE file", cmd );
    gsArgSwitch arg_dirich("n", "nitsche", "Use the Nitsche's method for Dirichlet sides", cmd);
    gsArgSwitch arg_dg("d", "discGalerkin", "Use Discontiouous Galerkin method for patch interfaces", cmd);
    gsArgSwitch a4("", "plot", "Plot result in ParaView format", cmd);
    gsArgVal<int> arg_plot_pts("s","plotSamples", 
		     "Number of sample points to use for plotting", 
		     false,1000, "sampling points for plotting", cmd );

    gsArgVal<int> a3("e","degreeElevation", 
		     "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", 
		     false,-1, "degree elevation steps", cmd );
    gsArgVal<int> a2("r","uniformRefine", 
		     "Number of Uniform h-refinement steps to perform before solving", 
		     false,2, "refinement steps", cmd );
    gsArgVal<std::string> a5("b","basis","File containing basis for discretization (.xml)", 
			     false,"", "basis file", cmd );
    gsArgVal<std::string> a1("g","geometry","File containing Geometry (.xml, .axl, .txt)", 
			     false, "", "geometry file", cmd );
    
    cmd.parse(argc,argv);
    std::string fn = a1.getValue();
    std::string fn_pde = a6.getValue();
    std::string fn_basis = a5.getValue();

    numRefine  = a2.getValue();
    numElevate = a3.getValue();
    if ( arg_dirich.getValue() )
        Dirichlet  = 1;
    else
        Dirichlet  = 0;

    if ( arg_dg.getValue() )
    {
        DG  = 1;
        Dirichlet  = 1;// use Nitsche for boundary
    }
    else
        DG  = 0;

    plot       = a4.getValue();
    plot_pts   = arg_plot_pts.getValue();

    if ( ! fn_basis.empty() )
	{
    gsBasis<> * bb = gsReadFile<>( fn_basis );
    cout << "Got basis: "<< * bb<<"\n"; 
    bases.push_back(bb);
    //cout << "Warning: basis ignored.\n"; 
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
	    fn_pde = GISMO_SOURCE_DIR;
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
			    fn_pde+="/filedata/pde/poisson1d_sin.xml";
			    break;
			case 2:
			    fn_pde+="/filedata/pde/poisson2d_sin.xml";
			    break;
			case 3:
			    fn_pde+="/filedata/pde/poisson3d_sin.xml";
			    break;
			default:
			    return false;
			}
		}
	    else
		fn_pde+="/filedata/pde/poisson2d_sin.xml";
	}
    ppde = gsReadFile<>(fn_pde);
    if ( !ppde )
      {
	gsWarn<< "Did not find any PDE in "<< fn<<", quitting.\n";
	return false;
      }
    
    if ( fn.empty() )
	{
	  fn = GISMO_SOURCE_DIR;
      switch ( ppde->m_compat_dim )
	    {
	    case 1:
	      fn+= "/filedata/bspline1d.xml";
	      break;
	    case 2:
	      fn+= "/filedata/square.xml";
	      break;
	    case 3:
	      fn+= "/filedata/cube.xml";
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
    
  } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << "\n"; return false; }
  return true;
}
