// example solving the poisson equation with gismo

#include <iostream>
#include <ctime>

#include <gismo.h>


using namespace std;
using namespace gismo;

int main(int argc, char *argv[])
{
    int numRefine; // defaults to 2
    int numElevate; // defaults to 0
    gsMultiPatch<> * geo; // defaults to BSplineCube
    bool plot; // If set to true, paraview file is generated and launched on exit
    //     try 
    //     {
    // 	gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    // 	gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
    // 	cmd.parse(argc,argv);
    // 	plot       = ap.getValue();
    //     } catch ( gsArgException& e )
    //     { cout << "Error: " << e.error() << " " << e.argId() << endl; return -1; }

    //---------------------------------------------------------------------------//
    try 
    {
        gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");    
        gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
        gsArgVal<int> a1("r","uniformRefine", 
                "Number of Uniform h-refinement steps to perform before solving", 
                false,2, "refinement steps", cmd );
        std::string fn;
        gsArgVal<std::string> a2("g","geometry","File containing Geometry (.axl, .txt)", 
                false,"", "geometry file", cmd );
        gsArgVal<int> a3("e","degreeElevation", 
                "Number of degree eleveation steps to perform before solving", 
                false,0, "degree elevation steps", cmd );
        gsArgVal<std::string> a4("p","pde","File containing a poisson PDE (.xml)", 
			     false,"", "PDE file", cmd );	

        cmd.parse(argc,argv);
        numRefine = a1.getValue();
        fn        = a2.getValue();
        numElevate = a3.getValue();
	std::string fn_pde = a4.getValue();
        plot = ap.getValue();
        if (fn.empty() )
            //             geo = gsMultiPatch<>( gsNurbsCreator<>::BSplineHalfCube(2,0.5,0.5,0.5) );
            geo = new gsMultiPatch<>( *safe(gsNurbsCreator<>::NurbsSphere()) );
        else
            geo = gsReadFile<>( fn ) ;

        if (numRefine<0)
        { 
            std::cout<<"Number of refinements must be non-negative, quitting.\n"; 
            return -1;
        }
        if (numElevate<0)
        { 
            std::cout<<"Number of elevations must be non-negative, quitting.\n"; 
            return -1;
        }

    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; return -1; }


    // Source function
    //  The source function is taken from the Spherical harmonic table
    //   on Wikipedia. We consider the order(I) = 3
    gsFunctionExpr<> f("3*sqrt(7/pi)*z*(2*z^2 -3*x^2-3*y^2)") ;
    cout<<"Source function "<< f <<".\n" << endl;

    // Exact solution
    // The exact must be computed with the quadrature node on the control points
    // This is a NOT exact RHS
    // currently needs the exact RHS to 
    gsFunctionExpr<> g("0.25*sqrt(7/pi)*z*(2*z^2 -3*x^2-3*y^2)") ;
    //gsFunctionExpr<> g("x+y") ;
    cout<<"Exact solution "<< g <<".\n" << endl;

    //---------------------------------------------------------//

    gsBVProblem<> bvp( *geo, new gsSurfacePoissonPde<real_t>( f, geo->geoDim() ) );
    //     for (gsMultiPatch<>::const_biterator bit = geo->bBegin(); bit != geo->bEnd(); ++bit)
    //     {
    //         cout<<"add_boundary.\n" << endl;
    //         bvp.addCondition( *bit, boundary::dirichlet, &g);
    //     }

    gsBasis<> * tbasis = geo->patch(0).basis().clone();

    // The numRefine is for the number of refinements
    for (int i = 0; i < numRefine; ++i)
        tbasis->uniformRefine();

    std::vector<gsBasis<> * > all_bases;
    for ( unsigned i = 0; i!= geo->nPatches(); ++i )
        all_bases.push_back( tbasis );

    gsGalerkinMethod<>  surfacepoisson( bvp , all_bases );

    // Refine the solution space
    //solver.basis().uniformRefine();

    cout<<"-Discretization Space: \n"<< *tbasis << endl;  

    // Compute solution field
    gsField<> * x = surfacepoisson.solve();

    if (plot)
    {
        std::cout<<"Plotting in Paraview...\n";
        gsWriteParaview<>( *x, "surfacepoisson_sol",1000) ;

        //Plot exact solution in paraview
        gsField<> exact( *geo, g, false ) ;
        gsWriteParaview<>( exact, "surfacepoisson_exact_sol", 1000) ;

        //run: paraview paraview_poisson.vts
        char cmd[100];
        if ( geo->nPatches() > 1 )
            strcpy(cmd,"paraview surfacepoisson_sol.pvd \0");
        else
            strcpy(cmd,"paraview surfacepoisson_sol.vts \0");
        strcat(cmd," &");
        return system(cmd);
    }
    delete tbasis;
    delete x;
    delete geo;

    return 0;
}
