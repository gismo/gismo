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
    bool plot;

    try 
    {
        gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");
        gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
        gsArgVal<int> a1("r","uniformRefine", 
                "Number of Uniform h-refinement steps to perform before solving",
                false, 3, "refinement steps", cmd );
        cmd.parse(argc,argv);

        numRefine = a1.getValue();
        plot       = ap.getValue();
        if (numRefine<0)
        { 
            std::cout<<"Number of refinements must be non-negative, quitting.\n"; 
            return -1;
        }
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; return -1; }

    // Source function
    gsFunctionExpr<> f("2*pi^2*sin(pi*x)*sin(pi*y)") ; 
    // Exact solution
    gsFunctionExpr<> g    ("sin(pi*x) * sin(pi*y)"); 

    cout<<"Source function "<< f <<".\n" << endl;
    cout<<"Exact solution "<< g <<".\n" << endl;

    // Define Geometry
    gsGeometry<> * geo = gsNurbsCreator<>::BSplineSquare(2);
    //gsGeometry<> * geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus();
    //gsGeometry<> * geo = gsNurbsCreator<>::BSplineQuarterAnnulus();
    //gsGeometry<> * geo = gsNurbsCreator<>::NurbsQuarterAnnulus();

    // Setup the Poisson Boundary Value problem
    gsBVProblem<> bvp( geo, new gsPoissonPde<>(f, geo->geoDim()) );
    // Dirichlet BCs
    gsFunctionExpr<> g_dir("sin(pi*x) * sin(pi*y)");
    bvp.addCondition( boundary::west,  boundary::dirichlet, &g_dir);
    //bvp.addCondition( boundary::east,  boundary::dirichlet, &g_dir);
    //bvp.addCondition( boundary::north, boundary::dirichlet, &g_dir);
    //bvp.addCondition( boundary::south, boundary::dirichlet, &g_dir);

    // Neumann BCs (Note: pure Neumann problem is singular.)
    gsFunctionExpr<> gx (" pi*sin(pi*y) * cos(pi*x)");
    gsFunctionExpr<> mgx("-pi*sin(pi*y) * cos(pi*x)");
    gsFunctionExpr<> gy (" pi*sin(pi*x) * cos(pi*y)");
    gsFunctionExpr<> mgy("-pi*sin(pi*x) * cos(pi*y)");

    //bvp.addCondition( boundary::west,  boundary::neumann, &mgx );
    bvp.addCondition( boundary::east,  boundary::neumann, &gx  );
    bvp.addCondition( boundary::north, boundary::neumann, &gy  );
    bvp.addCondition( boundary::south, boundary::neumann, &mgy );

    // Define discretization space by refining the basis of the geometry
    gsBasis<> * tbasis =  geo->basis().clone() ;
    for (int i = 0; i < numRefine; ++i)
      tbasis->uniformRefine();

    cout<<"Discretization Space: dim=" << tbasis->dim() << " deg=" << geo->basis().minDegree() << " dofs=" << tbasis->size() << endl;

    gsGalerkinMethod<>  poisson( bvp, *tbasis );

    cout << "Dofs: " << poisson.dofs() << endl;
    cout << "System size: " << poisson.freeDofs() << endl;

    // Compute solution field
    gsStopwatch time;
    poisson.assemble();
    const double assmTime = time.stop();
    gsField<>* sol = poisson.solve();
    const double solveTime = time.stop();
    cout << "Assembling time: " <<  assmTime << " s" << endl;
    cout << "Solver time:     " <<  solveTime - assmTime << " s" << endl;
    cout << "Total time:      " <<  solveTime << " s" << endl;

    cout << "L2 error: " << computeL2Distance(*sol, g, false, 3*tbasis->size()) << endl;

    // Optionally plot solution in paraview
    int result = 0;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        std::cout<<"Plotting in Paraview...\n";
        gsWriteParaview<>( *sol, "poisson2d", 1000);
        gsField<> exact( bvp.patches(), g, false );
        gsWriteParaview<>( exact, "poisson2d_exact", 1000);

        // Run paraview
        result = system("paraview poisson2d.vts &");
    }

    delete tbasis;
    delete sol;

    return result;
}
