
#include <iostream>
#include <gismo.h>

#include "gsSemiRegularBasis.h"

#include "gsSemiRegularBezier.h"

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

gsMatrix<real_t> DistortInnerPoints(const gsMatrix<real_t>& cp, real_t amount)
{
    gsMatrix<real_t> distortedCP = cp;

    for (index_t i = 1; i < distortedCP.cols() - 1; ++i) // Skip first and last point
    {
        if ((i - 1) % 2 == 0) // First inner point: add, second: subtract, etc.
        {
            distortedCP(0, i) += amount;
            distortedCP(1, i) += amount;
        }
        else
        {
            distortedCP(0, i) -= amount;
            distortedCP(1, i) -= amount;
        }
    }

    return distortedCP;
}

std::vector<double> generateKnotValues(int ll, int p)
{
    int numIntervals = std::pow(2, ll);
    std::vector<double> knots;

    // Append p+1 zeros.
    knots.insert(knots.end(), p + 1, 0.0);

    // Append the interior knots: i/(2^ll) for i=1 to 2^ll - 1.
    for (int i = 1; i < numIntervals; ++i)
    {
        knots.push_back(i / static_cast<double>(numIntervals));
    }

    // Append p+1 ones.
    knots.insert(knots.end(), p + 1, 1.0);
    return knots;
}




int main(int argc, char* argv[])
{
    index_t p = 3, l = 3, s = 5041;
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

    /*
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


    */

    // -----------------------------------------------------------------------------------------------------
    // Geometry parametrization: Parametrization of the collapsed triangle as a Bezier surface
    // -----------------------------------------------------------------------------------------------------

    gsVector<real_t> A(2), B(2);
    A << 1.0, 0.0;
    B << 0.0, 1.0;
    gsMatrix<real_t> BoundaryCP = StraightLinePoints(3, A, B); //Boundary Control Points
    gsInfo << "Line points:\n" << BoundaryCP << "\n";

    //========================= Experiment 1 =========================
    //Scale the control points sampled from a straight line

    //Perform the linear scaling with the correct ordering of the COntrol Points
    const index_t rows = BoundaryCP.rows();
    const index_t cols = BoundaryCP.cols();

    gsMatrix<real_t> ScaledCP(rows, cols * (p + 1));


    for (index_t j = 0; j < cols; ++j)
    {
        for (index_t k = 0; k < (p + 1); ++k)
        {
            real_t scale = (1.0 * k) / p;
            gsVector<real_t> scaledCol = scale * BoundaryCP.col(j);


            ScaledCP.col(j * (p + 1) + k) = scaledCol;
        }
    }

    gsInfo << "Scaled Control Points:\n" << ScaledCP << "\n";

    //========================= Experiment 2 =========================
    //Distort only the outer boundary control points from 2nd to pre-last
    gsMatrix<real_t> DistortedCP = ScaledCP;

    real_t amount = 0.15;
    gsMatrix<real_t> OuterCP(2, 4);
    for (index_t j=0; j<((p +1)* (p + 1)); j=j+(p+1))
    {
        OuterCP.col(j/(p+1)) = ScaledCP.col(j + p);
    }

 
    gsMatrix<real_t> DistortedOuterCP = DistortInnerPoints(OuterCP, amount);
    //Replace the outer control points in DistortedCP with the distorted ones
    for (index_t j=0; j<((p +1)* (p + 1)); j=j+(p+1))
    {
        DistortedCP.col(j + p) = DistortedOuterCP.col(j/(p+1));
    }


    //+========================== Experiment 3 =========================
    //Distort and afterwards scale
    //Scale the DistortedOuterCP
    gsMatrix<real_t> ScaledDistortedOuterCP(rows, cols * (p + 1));
    for (index_t j = 0; j < DistortedOuterCP.cols(); ++j)
    {
        for (index_t k = 0; k < (p + 1); ++k)
        {
            real_t scale = (1.0 * k) / p;
            gsVector<real_t> scaledCol = scale * DistortedOuterCP.col(j);

            ScaledDistortedOuterCP.col(j * (p + 1) + k) = scaledCol;
        }
    }


    //=================== Experiment 4 =========================
    //Reduced Bezier control net with equidistant points
    //allocate gsmatrix ReducedCP of size 2x10
    gsMatrix<real_t> ReducedCP(2, 10);
    ReducedCP << 0.0, 0.333333, 0.666667, 1.0, 0.333333, 0.666667, 0.333333, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.333333, 0.333333, 0.666667, 0.333333, 0.666667, 1.0;



    //create the semi regular Bezier basis
    gsSemiRegularBezier<2> semibezier(3);


    gsInfo << "The basis semi Bezier: \n" << semibezier << "\n";


    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << semibezier.dim() << "\n"
           << "Number of basis functions: "        << semibezier.size() << "\n"
           << "\n";

    //========================= Choose the control points from the desired Experiment as the FinalCP =========================

    gsMatrix<real_t> FinalCP = ScaledCP;


    //Create the Bezier basis
    gsKnotVector<> knotU_bezier(0,1,0,p+1);
    gsKnotVector<> knotV_bezier(0,1,0,p+1);

    gsTensorBSplineBasis<2, real_t> bezier_basis(knotU_bezier, knotV_bezier);

    gsMatrix<real_t> coefs(FinalCP.cols(), 3);  // coefs will store the ScaledCP transposed (as accepted by the bSplineSurface constructor)

    for (index_t col = 0; col != FinalCP.cols(); col++) {
        real_t x = FinalCP(0, col);  // x-coordinate from ScaledCP
        real_t y = FinalCP(1, col);  // y-coordinate from ScaledCP
    
        coefs(col, 0) = x; 
        coefs(col, 1) = y;  
        coefs(col, 2) = 0; // Set z-coordinate to 0 in the third column to have a planar collapsed triangle
    }

    gsInfo << "Coefs:\n" << coefs << "\n";

    gsTensorBSpline<2, real_t>  collapsed_triangle(bezier_basis, coefs);


    gsMatrix<real_t> Reducedcoefs(ReducedCP.cols(), 3);  // coefs will store the ScaledCP transposed (as accepted by the bSplineSurface constructor)

    for (index_t col = 0; col != ReducedCP.cols(); col++) {
        real_t x = ReducedCP(0, col);  // x-coordinate from ScaledCP
        real_t y = ReducedCP(1, col);  // y-coordinate from ScaledCP
    
        Reducedcoefs(col, 0) = x; 
        Reducedcoefs(col, 1) = y;  
        Reducedcoefs(col, 2) = 0; // Set z-coordinate to 0 in the third column to have a planar collapsed triangle
    }


    auto reduced_triangle = semibezier.makeGeometry( Reducedcoefs );
    // collapsed_triangle is now a std::unique_ptr<gsGeometry<real_t>>
    // pointing to the correct concrete type (e.g. a custom gsBezierGeometry)
    // and you can use it like any other gsGeometry:
    gsInfo 
    << " REDUCED param dim = " << reduced_triangle->parDim()
    << ", REDUCED geo dim = "   << reduced_triangle->geoDim() << "\n"
    << " REDUCED coefs  =\n"      << reduced_triangle->coefs()  << "\n";



    

    // A general 2d Geometry object (gsBSpline< T >) not a TensorBSpline object

    gsMatrix<> mapped_grid = collapsed_triangle.eval(u);

    //mapped_grid_2D is needed for function evaluation
    gsMatrix<real_t> mapped_grid_2D(2, mapped_grid.cols());

    for (index_t i = 0; i < mapped_grid.cols(); ++i) {
        mapped_grid_2D.row(0).col(i) = mapped_grid.row(0).col(i);
        mapped_grid_2D.row(1).col(i) = mapped_grid.row(1).col(i);
    }

    //gsInfo<<mapped_grid_2D << "\n";


    // ----------------------------------------------------------------------
    // fitting
    // ----------------------------------------------------------------------
    //gsFunctionExpr<real_t> ff("x", "y", "sin(10*x)*sin(10*y)", 2);
    gsFunctionExpr<real_t> ff("x", "y", "sin(10*x) + sin(10*y)", 2);
    gsMatrix<> mapped_ff_values = ff.eval(mapped_grid_2D);


    

    gsFitting<> fit(u, mapped_ff_values, BB);
    fit.compute();

    gsInfo << *fit.result() << "\n";
    gsInfo << fit.result()->coefs() <<"\n";

    auto res = fit.result()->eval(u);

    //auto error_L2 = ( mapped_ff_values - res ).norm();
    auto max_error = ( mapped_ff_values - res ).cwiseAbs().maxCoeff();

    gsInfo<<"Error in max norm is: "<<max_error<< "\n";

    //get the basis of fit object

    std::vector<double> dyadic_knots = generateKnotValues(l, p);
    gismo::gsKnotVector<double> kv(dyadic_knots, p);
    gsTensorBSplineBasis<2, real_t> TB(kv, kv);


  


   if (plot)
   {
       gsWriteParaview(bezier_basis, "bezier_basis", 10000);
       gsWriteParaview(collapsed_triangle, "collapsed_triangle", 10000);

       gsMesh<> mesh;
       collapsed_triangle.controlNet(mesh);
       gsWriteParaview(mesh, "bezier_control_net", 10000);

       gsWriteParaviewPoints(mapped_grid, "mapped_grid");

       gsWriteParaviewPoints(mapped_ff_values, "mapped_ff_values");

       gsWriteParaview(*fit.result(), "srfitting", 10000);

       gsWriteParaview(semibezier, "semibezier_basis", 10000);

       gsWriteParaview(TB,"tensor_basis", 10000);

       gsWriteParaview(BB, "semiregular_basis", 10000);

   }


    return EXIT_SUCCESS;
}
