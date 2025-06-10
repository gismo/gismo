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


int main(int argc, char* argv[])
{
    index_t p = 3, s = 300000;

    //create the point grid
    gsVector<> a(2); a.setZero();
    gsVector<> b(2); b.setOnes();
    gsVector<unsigned> np = uniformSampleCount(a,b,s);
    gsMatrix<> u = gsPointGrid(a, b, np);


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



    //========================= Choose the control points from the desired Experiment as the FinalCP =========================

    gsMatrix<real_t> FinalCP = ScaledDistortedOuterCP;



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

    gsInfo << "Control points:\n" << coefs << "\n";

    gsTensorBSpline<2, real_t>  collapsed_triangle(bezier_basis, coefs);

    gsMatrix<> mapped_grid = collapsed_triangle.eval(u);

    //mapped_grid_2D is needed for function evaluation
    gsMatrix<real_t> mapped_grid_2D(2, mapped_grid.cols());

    for (index_t i = 0; i < mapped_grid.cols(); ++i) {
        mapped_grid_2D.row(0).col(i) = mapped_grid.row(0).col(i);
        mapped_grid_2D.row(1).col(i) = mapped_grid.row(1).col(i);
    }


    // ----------------------------------------------------------------------
    // fitting
    // ----------------------------------------------------------------------
    gsFunctionExpr<real_t> ff("x", "y", "sin(10*x) + sin(10*y)", 2);
    gsMatrix<> mapped_ff_values = ff.eval(mapped_grid_2D);

    /*
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
        gsFitting<> fit(u, mapped_ff_values, TB);
        fit.compute();

        auto res = fit.result()->eval(u);
        auto error_max = ( mapped_ff_values - res ).cwiseAbs().maxCoeff();
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile1 << l << " " << error_max << "\n";


    }
    errorFile1.close();
    gsInfo << "All errors written to scaled_cp_tensor_prod_errors.txt\n";

    */


  
    gsInfo << "\n \n ====================Fitting with Semi-Regular B-Splines====================\n";
    std::ofstream errorFile2("scaled_cp_semi_regular_errors.txt");
    if (!errorFile2.is_open())
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
        gsFitting<> fit(u, mapped_ff_values, BB);
        fit.compute();
 
    
        auto res = fit.result()->eval(u);
    
        auto error_max = ( mapped_ff_values - res ).cwiseAbs().maxCoeff();
        fit.computeMaxNormErrors();
        //auto error_new = fit.maxPointError(); 
    
        gsInfo<<"Error in max norm is: "<<error_max<< " for level l="<<l<< "\n\n";
        errorFile2 << l << " " << error_max << "\n";

    }

    errorFile2.close();
    gsInfo << "All errors written to scaled_cp_semi_regular_errors.txt\n";


 

  


    return EXIT_SUCCESS;
}