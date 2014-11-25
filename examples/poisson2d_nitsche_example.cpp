// example solving the poisson equation with gismo

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>


using namespace std;
using namespace gismo;


int main(int argc, char *argv[])
{
    int numRefine;// defaults to 3
    int degree;
    bool plot;
    bool noNitsche;

    try 
    {
        gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization and Nitsche's method for Dirichlet boundaries.");    
        gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
        gsArgSwitch ano("", "no-nitsche", "Disable Nitsche method and use elimination of Dirichlet dofs", cmd);
        gsArgVal<int> a1("r","uniformRefine", 
                         "Number of Uniform h-refinement steps to perform before solving", 
                         false,3, "refinement steps", cmd );
        gsArgVal<int> arg_deg("p", "degree",
                "Degree of the B-spline discretization space",
                false, 2, "int", cmd);
        cmd.parse(argc,argv);
        numRefine  = a1.getValue();
        plot       = ap.getValue();
        noNitsche  = ano.getValue();
        degree     = arg_deg.getValue();
        if (numRefine<0)
        { 
            std::cout<<"Number of refinements must be non-negative, quitting.\n"; 
            return -1;
        }
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; return -1; }

/*
    // Source function
    gsFunctionExpr<> f("(8-9*sqrt(x*x + y*y)*sin(2*atan(y/x))) / (x*x+y*y) ");
    // Exact solution
    gsFunctionExpr<> g    ("(x*x+y*y-3*sqrt(x*x+y*y)+2)*sin(2*atan(y/x))");
    // Boundary condition
    gsFunctionExpr<> g_dir("(x*x+y*y-3*sqrt(x*x+y*y)+2)*sin(2*atan(y/x))");// zero for nurbs ring
*/

    // Source function
    gsFunctionExpr<> f("2*pi^2*sin(pi*x)*sin(pi*y)") ;
    // Exact solution
    gsFunctionExpr<> g    ("sin(pi*x) * sin(pi*y)");
    // Boundary condition
    gsFunctionExpr<> g_dir("sin(pi*x) * sin(pi*y)");

    cout<<"Source function "<< f <<".\n" << endl;
    cout<<"Exact solution "<< g <<".\n" << endl;

    // Define Geometry
    gsKnotVector<> KV (0.0,1.0, 0,3) ;
    gsMatrix<> C(9,2);
    C  << 1 , 0 , 1.5,  0  , 2 , 0  
        , 1 , 1 , 1.5, 1.5 , 2 , 2 
        , 0 , 1 ,  0 , 1.5 , 0 , 2 ;

    // gsKnotVector<> KV (0,1,0,2) ;
    // gsMatrix<> C(4,2) ;
    // C << 0,0,  1,0,
    //      0,1,  1,1;

    gsGeometry<> * geo = new gsTensorBSpline<2,real_t>(KV,KV, give(C));
    
    // test for boundary moment
    //geo->uniformRefine(3);
    //gsGaussAssembler<real_t >Angl(geo);
//     gsVector<>  *BM= Angl.boundaryMoments(&geo->basis().component(),gsFunctionExpr<>("1"), boundary::west);
//     cout<<"Boundary Moments" << *BM << endl;
//     return 0;
//     
    //
       
    // Setup the Poisson Boundary Value problem
    gsBVProblem<> bvp( geo, new gsPoissonPde<>( f, geo->geoDim() ) );

    // Dirichlet
//     bvp.addCondition( boundary::west,  boundary::dirichlet, 0 );
//     addCondition(int p, boundary::side s, boundary::type t, gsFunction<> * f)
    bvp.addCondition( boundary::west,  boundary::dirichlet, &g_dir);
    bvp.addCondition( boundary::east,  boundary::dirichlet, &g_dir);
    bvp.addCondition( boundary::north, boundary::dirichlet, &g_dir);
    bvp.addCondition( boundary::south, boundary::dirichlet, &g_dir);

    // Neumann
    //gsFunctionExpr<> gx("0") ;//-pi*sin(pi*y) * cos(pi*y)
    //gsFunctionExpr<> gy("0") ;//-pi*sin(pi*x) * cos(pi*y)
    //bvp.addCondition( boundary::north, boundary::neumann, &gx );
    //bvp.addCondition( boundary::south, boundary::neumann, &gy );

    // Define descritization space by refining the basis of the geometry  
    gsBasis<> * tbasis =  geo->basis().clone() ;
    tbasis->setDegree( degree );
    for (int i = 0; i < numRefine; ++i)
      tbasis->uniformRefine();

    cout<<"Discretization Space: dim=" << tbasis->dim() << " deg=" << tbasis->minDegree() << " dofs=" << tbasis->size() << endl;

    // Use Nitsche method
    gsGalerkinMethod<>  poisson( bvp, *tbasis, noNitsche ? dirichlet::elimination 
                                                         : dirichlet::nitsche );

    // Compute solution field
    gsStopwatch time;
    poisson.assemble();
    const double assmTime = time.stop();
    gsField<>* x = poisson.solve();
    const double solveTime = time.stop();
    cout << "Assembling time: " <<  assmTime << " s" << endl;    
    cout << "Solver time:     " <<  solveTime - assmTime << " s" << endl;    
    cout << "Total time:      " <<  solveTime << " s" << endl;    

    cout << "L2 error: " << 
          x->distanceL2( g ) << endl;

    // Optionally plot solution in paraview
    if (plot)
    {
        std::cout<<"Plotting in Paraview...\n";
        gsWriteParaview<>( *x, "poisson2d",1000) ;

        //Plot exact solution in paraview
        gsField<> exact( bvp.patches(), g, false ) ;
        gsWriteParaview<>( exact, "poisson2d_exact", 1000) ;

        delete tbasis;
        delete x;

        //run: paraview paraview_poisson.vts
        char cmd[100];
        strcpy(cmd,"paraview poisson2d.vts \0");
        strcat(cmd," &");
        return system(cmd);
    }

    delete tbasis;
    delete x;

    return 0;
}
