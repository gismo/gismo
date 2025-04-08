
#include <iostream>
#include <gismo.h>

#include "gsSemiRegularBasis.h"

using namespace gismo;


// Generates p+1 evenly spaced points on a straight line from lineStart to lineEnd
gsMatrix<real_t> StraightLinePoints(index_t p, const gsVector<real_t> &lineStart, const gsVector<real_t> &lineEnd)
{
    
    gsMatrix<real_t> points(2, p + 1);

    for (index_t k = 0; k <= p; ++k)
    {
        real_t t = (1.0 * k) / p;
        gsVector<real_t> pt = (1.0 - t) * lineStart + t * lineEnd; //Straight line equation

        points.col(k) = pt;
    }

    return points;
}


int main(int argc, char* argv[])
{
    index_t p = 3, l = 2, s = 3;
    bool plot = false;
    gsCmdLine cmd("Tutorial on gsBasis class.");
    cmd.addInt   ("p", "degree", "Degree", p);
    cmd.addInt   ("l", "layers", "layers", l);
    cmd.addInt   ("s", "samples", "number of sample points", s);
    cmd.addSwitch("plot"   , "Plot the result", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // making the basis
    // ======================================================================

    gsSemiRegularBasis<2> BB(p, l);

    // ======================================================================
    // printing some information about the basis
    // ======================================================================

    // printing the basis
    gsInfo << "The basis: \n" << BB << "\n";


    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << BB.dim() << "\n"
           << "Number of basis functions: "        << BB.size() << "\n"
           << "\n";


    // support of the basis
    // (dim x 2 matrix, the parametric domain)
    //gsMatrix<> support = BB.support();


    // ======================================================================
    // basis evaluation
    // ======================================================================
    gsVector<> a(2); a.setZero();
    gsVector<> b(2); b.setOnes();
    gsVector<unsigned> np = uniformSampleCount(a,b,s);
    gsMatrix<> u = gsPointGrid(a, b, np);
    if (s<100)
        gsInfo << "Points (" << u.cols() << "): \n" << u << "\n\n";

    // indices of active (nonzero) functions at parameter u
    gsMatrix<index_t> active = BB.active(u);
    if (s<100)
        gsInfo << "#Active basis functions at u: \n"
               << active << "\n\n";

    // values of all active functions at u
    gsMatrix<> values = BB.eval(u);
    if (s<100)
        gsInfo << "Values at u : \n"
               << values << "\n\n";

    // values of single basis functions
    if (s<100)
    {
        for (index_t i = 0; i != active.rows(); i++)
        {
            gsMatrix<> val = BB.evalSingle(active(i), u);
            
            gsInfo << "basis fun. index:  " << active(i)
                   << "   value: " << val(0, 0) << "\n";
        }
        gsInfo << "\n";
    }

    // ----------------------------------------------------------------------
    // derivatives
    // ----------------------------------------------------------------------


    gsMatrix<> derivs = BB.deriv(u);
    if (s<100)
        gsInfo << "Derivatives at u : \n"
               << derivs << "\n\n";

    // ----------------------------------------------------------------------
    // second derivatives
    // ----------------------------------------------------------------------

    gsMatrix<> derivs2 = BB.deriv2(u);
    if (s<100)
        gsInfo << "Second derivatives at u : \n"
               << derivs2 << "\n\n";

    gsInfo << "\nFor more information about evaluation "
           << "(and order of derivatives) look at doxygen documentation."
           << "\n\n";

    // -----------------------------------------------------------------------------------------------------
    // Geometry parametrization: Parametrization of the collapsed triangle as a Bezier surface
    // -----------------------------------------------------------------------------------------------------

    gsVector<real_t> A(2), B(2);
    A << 1.0, 0.0;
    B << 0.0, 1.0;
    gsMatrix<real_t> BoundaryCP = StraightLinePoints(3, A, B); //Boundary Control Points
    gsInfo << "Line points:\n" << BoundaryCP << "\n";

    //Perform the linear scaling
    const index_t rows = BoundaryCP.rows();
    const index_t cols = BoundaryCP.cols();

    gsMatrix<real_t> ScaledCP(rows, cols * (p + 1));

    for (index_t k = 0; k <= p; ++k)
    {
        // Perform linear scaling of BoundaryCP
        real_t scale = (1.0 * k) / p;
        gsMatrix<real_t> scaled = scale * BoundaryCP;

        ScaledCP.block(0, k * cols, rows, cols) = scaled;
    }


    //Create the Bezier basis
    gsKnotVector<> knotU_bezier(0,1,0,p+1);
    gsKnotVector<> knotV_bezier(0,1,0,p+1);

    gsTensorBSplineBasis<2, real_t> bezier_basis(knotU_bezier, knotV_bezier);

    gsMatrix<real_t> coefs(ScaledCP.cols(), 3);  // coefs will store the ScaledCP transposed (as accepted by the bSplineSurface constructor)

    for (index_t col = 0; col != ScaledCP.cols(); col++) {
        real_t x = ScaledCP(0, col);  // x-coordinate from ScaledCP
        real_t y = ScaledCP(1, col);  // y-coordinate from ScaledCP
    
        coefs(col, 0) = x; 
        coefs(col, 1) = y;  
        coefs(col, 2) = 0; // Set z-coordinate to 0 in the third column to have a planar collapsed triangle
    }

    gsTensorBSpline<2, real_t>  collapsed_triangle(bezier_basis, coefs);

    gsMatrix<> mapped_grid = collapsed_triangle.eval(u);

    //mapped_grid_2D is needed for function evaluation
    gsMatrix<real_t> mapped_grid_2D(2, mapped_grid.cols());

    for (index_t i = 0; i < mapped_grid.cols(); ++i) {
        mapped_grid_2D(0, i) = mapped_grid(0, i);  // x-coordinate
        mapped_grid_2D(1, i) = mapped_grid(1, i);  // y-coordinate
    }


    // ----------------------------------------------------------------------
    // fitting
    // ----------------------------------------------------------------------
    gsFunctionExpr<real_t> ff("x", "y", "x^2 + y^3", 2);
    gsMatrix<> mapped_ff_values = ff.eval(mapped_grid_2D);


    /*
    //print mapped_grid_2D in visible format
    gsInfo << "Mapped grid 2D: \n" << mapped_grid_2D << "\n";
    gsInfo << "Mapped ff values: \n" << mapped_ff_values << "\n";

    gsFitting<> fit(mapped_grid_2D, mapped_ff_values, BB);
    fit.compute();
    gsInfo << *fit.result() <<"\n";

    auto res = fit.result()->eval(mapped_grid_2D);

    auto error_L2 = ( mapped_ff_values - res ).norm();

    gsInfo<<"Error is: "<<error_L2<< "\n";
    */
  


   if (plot)
   {
       gsWriteParaview(bezier_basis, "bezier_basis", 10000);
       gsWriteParaview(collapsed_triangle, "collapsed_triangle", 10000);

       gsMesh<> mesh;
       collapsed_triangle.controlNet(mesh);
       gsWriteParaview(mesh, "bezier_control_net", 10000);

       gsWriteParaviewPoints(mapped_grid, "mapped_grid");

       gsWriteParaviewPoints(mapped_ff_values, "mapped_ff_values");

       //gsWriteParaview(*fit.result(), "srfitting", 10000);
       
   }


    return EXIT_SUCCESS;
}
