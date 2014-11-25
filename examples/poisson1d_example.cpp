// example solving the poisson equation with gismo

#include <iostream>
#include <ctime>

#include <gismo.h>


using namespace std;
using namespace gismo;

int main()
{

    // We solve the equation
    // -u''= f on [0,2.5]
    // with boundary condition
    // u(0)=g(0), u(2.5)=g(2.5)=8.09..
    // where f = -g''
    // and g(x)=(7*x+2*x^2-3*x^3)*cos(x)*sin(x)


    // Exact solution
    gsFunctionExpr<> g("(7*x+2*x^2-3*x^3)*cos(x)*sin(x)") ;
    //gsFunctionExpr<> g("x") ;
    cout<<"Exact function "<< g <<".\n" << endl;

    // Source function
    gsFunctionExpr<> f("(6*x^3-4*x^2-23*x+2)*cos(2*x)-(18*x^2-8*x-14)*sin(2*x)") ;
    //gsFunctionExpr<> f("0") ;
    cout<<"Source function "<< f <<".\n" << endl;

    // Geometry
    // We have a parametrization F:[0,1]->[0, 2.5]
    // that describe the geometry using BSpline of degree
    // p=6
    // on the open grid of step 0.05
    int p = 6;
    gsMatrix<> C(11,1) ;
    C<< 0, .5, 1.1, 1.01, 1.5, 1.6, 1.7, 1.9, 2.2, 2.3,  2.5 ;
    gsKnotVector<> KV (0,1,4,p+1) ;
    gsBSpline<> * bsp = new gsBSpline<>(KV, give(C));

    // Setup the Poisson Boundary Value problem
    gsBVProblem<> bvp( bsp, new gsPoissonPde<>( f, bsp->geoDim() ) );
    bvp.addCondition( boundary::left, boundary::dirichlet, &g );
    bvp.addCondition( boundary::right,  boundary::dirichlet, &g );

    std::auto_ptr< gsBSplineBasis<> > b = safe( bsp->basis().clone() );
    // Refine the solution space
    b->uniformRefine();
    b->uniformRefine();
    
    gsGalerkinMethod<>  poisson( bvp, *b );

    cout<<"Discretization space on the parametric domain: \n"<< *b << endl;

    // Compute solution field
    double time;
    time = (double)clock()/CLOCKS_PER_SEC;
    std::auto_ptr< gsField<> > x = safe( poisson.solve() );
    time = (double)clock()/CLOCKS_PER_SEC - time;
    cout<<"Solver time: " <<  time << " ms" << endl;    

    gsGeometry<> * bs = dynamic_cast<gsGeometry<> *>(&x->function()) ;
    cout<<"Solution: \n" <<  bs->coefs()  << endl;    
    //delete bs; // bs is just a type-cast, no memory allocated for it
    delete bsp;
    return 0;
}
