
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
            distortedCP(0, i) += amount;
            distortedCP(1, i) += amount;
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
    index_t p = 3, l = 2, s = 300000;
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



    gsVector<real_t> A(2), B(2);
    A << 1.0, 0.0;
    B << 0.0, 1.0;
    gsMatrix<real_t> BoundaryCP = StraightLinePoints(3, A, B); //Boundary Control Points
    gsInfo << "Line points:\n" << BoundaryCP << "\n";


    //Experiment 1
    //Distort BoundaryCP from 2nd to pre-last then scale
    real_t amount = 0.15;
    gsMatrix<real_t> DistortedCP = DistortInnerPoints(BoundaryCP, amount);
    gsInfo << "Distorted Control Points:\n" << DistortedCP << "\n";

    //Perform the linear scaling with the correct ordering of the COntrol Points
    const index_t rows = BoundaryCP.rows();
    const index_t cols = BoundaryCP.cols();
    gsMatrix<real_t> ScaledCP(rows, cols * (p + 1));
    for (index_t j = 0; j < cols; ++j)
    {
        for (index_t k = 0; k < (p + 1); ++k)
        {
            real_t scale = (1.0 * k) / p;
            gsVector<real_t> scaledCol = scale * DistortedCP.col(j);

            ScaledCP.col(j * (p + 1) + k) = scaledCol;
        }
    }
    gsInfo << "Scaled Control Points:\n" << ScaledCP << "\n";

    /*
        d12── d13── d14── d15   ← Collapsed d0, d4, d8, d12
        │      │      │      │
        d8  ── d9  ── d10 ─ d11
        │      │      │      │
        d4  ── d5  ── d6  ─ d7
        │      │      │      │
        d0  ── d1  ── d2  ── d3

    */

    //Experiment 2
    //Leave the the distorted BoundaryCP unchanged and smoothen the inner rows of the control points
    gsMatrix<real_t> SmoothedCP(2, 16);
    SmoothedCP << 0.0, 0.333333, 0.666667, 1.0, 0.0, 0.222222, 0.444444, 0.816667, 0.0, 0.111111, 0.222222, 0.483333, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.111111, 0.222222, 0.483333 , 0.0, 0.222222, 0.444444, 0.816667, 0.0, 0.333333, 0.666667, 1.0;

    //Create the Bezier basis
    gsKnotVector<> knotU_bezier(0,1,0,p+1);
    gsKnotVector<> knotV_bezier(0,1,0,p+1);

    gsTensorBSplineBasis<2, real_t> bezier_basis(knotU_bezier, knotV_bezier);

    gsMatrix<real_t> scaled_coefs(ScaledCP.cols(), 3);  // coefs will store the ScaledCP transposed (as accepted by the bSplineSurface constructor)

    for (index_t col = 0; col != ScaledCP.cols(); col++) {
        real_t x = ScaledCP(0, col);  // x-coordinate from ScaledCP
        real_t y = ScaledCP(1, col);  // y-coordinate from ScaledCP
    
        scaled_coefs(col, 0) = x; 
        scaled_coefs(col, 1) = y;  
        scaled_coefs(col, 2) = 0; // Set z-coordinate to 0 in the third column to have a planar collapsed triangle
    }

    gsMatrix<real_t> smoothed_coefs(SmoothedCP.cols(), 3);  // coefs will store the ScaledCP transposed (as accepted by the bSplineSurface constructor)

    for (index_t col = 0; col != SmoothedCP.cols(); col++) {
        real_t x = SmoothedCP(0, col);  // x-coordinate from ScaledCP
        real_t y = SmoothedCP(1, col);  // y-coordinate from ScaledCP
    
        smoothed_coefs(col, 0) = x; 
        smoothed_coefs(col, 1) = y;  
        smoothed_coefs(col, 2) = 0; // Set z-coordinate to 0 in the third column to have a planar collapsed triangle
    }

    gsTensorBSpline<2, real_t>  scaled_triangle(bezier_basis, scaled_coefs);
    gsMatrix<> scaled_mapped_grid = scaled_triangle.eval(u);


    gsTensorBSpline<2, real_t>  smoothed_triangle(bezier_basis, smoothed_coefs);
    gsMatrix<> smoothed_mapped_grid = smoothed_triangle.eval(u);

    //mapped_grid_2D is needed for function evaluation
    gsMatrix<real_t> scaled_grid_2D(2, scaled_mapped_grid.cols());

    for (index_t i = 0; i < scaled_mapped_grid.cols(); ++i) {
        scaled_grid_2D.row(0).col(i) = scaled_mapped_grid.row(0).col(i);
        scaled_grid_2D.row(1).col(i) = scaled_mapped_grid.row(1).col(i);
    }

    gsMatrix<real_t> smoothed_grid_2D(2, smoothed_mapped_grid.cols());

    for (index_t i = 0; i < smoothed_mapped_grid.cols(); ++i) {
        smoothed_grid_2D.row(0).col(i) = smoothed_mapped_grid.row(0).col(i);
        smoothed_grid_2D.row(1).col(i) = smoothed_mapped_grid.row(1).col(i);
    }

     // ----------------------------------------------------------------------
    // fitting
    // ----------------------------------------------------------------------
    gsFunctionExpr<real_t> ff("x", "y", "sin(10*x + 10*y)", 2);
    gsMatrix<> scaled_ff_values = ff.eval(scaled_grid_2D);
    gsMatrix<> smoothed_ff_values = ff.eval(smoothed_grid_2D);

    gsInfo << "====================Fitting with Tensor-Product B-Splines====================\n";

    std::ofstream errorFile1("scaled_cp_tensor_prod_errors.txt");
    if (!errorFile1.is_open())
    {
        gsWarn << "Could not open file for writing!\n";
        return 1;
    }

    std::vector<int> dyadic_levels = {1, 2, 3, 4, 5, 6, 7, 8};
    for (int l:dyadic_levels) {

        std::vector<double> dyadic_knots = generateKnotValues(l, p);
        gismo::gsKnotVector<double> kv(dyadic_knots, p);
        gsTensorBSplineBasis<2, real_t> TB(kv, kv);

        //gsInfo << "The basis for level l=" << l << "\n" << TB << "\n";
        gsFitting<> fit(u, scaled_ff_values, TB);
        fit.compute();

        auto res = fit.result()->eval(u);
        auto error_max = ( scaled_ff_values - res ).cwiseAbs().maxCoeff();
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile1 << l << " " << error_max << "\n";


    }
    errorFile1.close();
    gsInfo << "All errors written to scaled_cp_tensor_prod_errors.txt\n";

    std::ofstream errorFile2("smoothed_cp_tensor_prod_errors.txt");
    if (!errorFile2.is_open())
    {
        gsWarn << "Could not open file for writing!\n";
        return 1;
    }

    for (int l:dyadic_levels) {

        std::vector<double> dyadic_knots = generateKnotValues(l, p);
        gismo::gsKnotVector<double> kv(dyadic_knots, p);
        gsTensorBSplineBasis<2, real_t> TB(kv, kv);

        //gsInfo << "The basis for level l=" << l << "\n" << TB << "\n";
        gsFitting<> fit(u, smoothed_ff_values, TB);
        fit.compute();

        auto res = fit.result()->eval(u);
        auto error_max = ( smoothed_ff_values - res ).cwiseAbs().maxCoeff();
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile2 << l << " " << error_max << "\n";


    }
    errorFile2.close();
    gsInfo << "All errors written to smoothed_cp_tensor_prod_errors.txt\n";



  
    gsInfo << "\n \n ====================Fitting with Semi-Regular B-Splines====================\n";
    std::ofstream errorFile3("scaled_cp_semi_regular_errors.txt");
    if (!errorFile3.is_open())
    {
        gsWarn << "Could not open file for writing!\n";
        return 1;
    }
    std::vector<int> semi_levels = {1, 2, 3, 4, 5, 6, 7, 8};
    for (int l:semi_levels)
    {

        gsSemiRegularBasis<2> BB(p, l); //create the semi-regular basis for each level l
        //printing the basis for each level l
        //gsInfo << "The basis for level l=" << l << "\n" << BB << "\n";
        gsFitting<> fit(u, scaled_ff_values, BB);
        fit.compute();
 
    
        auto res = fit.result()->eval(u);
    
        auto error_max = ( scaled_ff_values - res ).cwiseAbs().maxCoeff();
        fit.computeMaxNormErrors();
        //auto error_new = fit.maxPointError(); 
    
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile3 << l << " " << error_max << "\n";

    }

    errorFile3.close();
    gsInfo << "All errors written to scaled_cp_semi_regular_errors.txt\n";


    std::ofstream errorFile4("smoothed_cp_semi_regular_errors.txt");
    if (!errorFile4.is_open())
    {
        gsWarn << "Could not open file for writing!\n";
        return 1;
    }

    for (int l:semi_levels)
    {

        gsSemiRegularBasis<2> BB(p, l); //create the semi-regular basis for each level l
        //printing the basis for each level l
        //gsInfo << "The basis for level l=" << l << "\n" << BB << "\n";
        gsFitting<> fit(u, smoothed_ff_values, BB);
        fit.compute();
 
    
        auto res = fit.result()->eval(u);
    
        auto error_max = ( smoothed_ff_values - res ).cwiseAbs().maxCoeff();
        fit.computeMaxNormErrors();
        //auto error_new = fit.maxPointError(); 
    
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile4 << l << " " << error_max << "\n";

    }

    errorFile4.close();
    gsInfo << "All errors written to smoothed_cp_semi_regular_errors.txt\n";


    if (plot)
   {
       gsWriteParaview(bezier_basis, "bezier_basis", 10000);
       gsWriteParaview(scaled_triangle, "scaled_triangle", 10000);
       gsWriteParaview(smoothed_triangle, "smoothed_triangle", 10000);

       gsMesh<> mesh;
       scaled_triangle.controlNet(mesh);
       gsWriteParaview(mesh, "scaled_bezier_control_net", 10000);

       gsMesh<> mesh_smoothed;
       smoothed_triangle.controlNet(mesh_smoothed);
       gsWriteParaview(mesh_smoothed, "smoothed_bezier_control_net", 10000);


   }

   if (plot)
   {

       gsWriteParaviewPoints(scaled_mapped_grid, "scaled_mapped_grid");

       gsWriteParaviewPoints(smoothed_mapped_grid, "smoothed_mapped_grid");

       gsWriteParaviewPoints(scaled_ff_values, "scaled_ff_values");

       gsWriteParaviewPoints(smoothed_ff_values, "smoothed_ff_values");




   }




 

  




    return EXIT_SUCCESS;

}
