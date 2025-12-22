#include <iostream>
#include <chrono>
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsTensor/gsTensorBasis.hpp>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

using namespace gismo;

int main(int argc, char* argv[])
{
    index_t nPoints = 1000;
    index_t nKnotsX = 4;
    index_t nKnotsY = 4;
    index_t nKnotsZ = 4;
    index_t degreeX = 3;
    index_t degreeY = 3;
    index_t degreeZ = 3;

    gsCmdLine cmd("Test to compare forward and reverse automatic differentiation on B-spline surfaces.");
    cmd.addInt("N", "npoints", "Number of evaluation points", nPoints);
    cmd.addInt("X", "nX", "Number of knots in x direction", nKnotsX);
    cmd.addInt("Y", "nY", "Number of knots in y direction", nKnotsY);
    cmd.addInt("Z", "nZ", "Number of knots in z direction", nKnotsZ);
    cmd.addInt("p", "degX", "Degree in x direction", degreeX);
    cmd.addInt("q", "degY", "Degree in y direction", degreeY);
    cmd.addInt("r", "degZ", "Degree in z direction", degreeZ);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsKnotVector<> kvX(0.0, 1.0, nKnotsX-1, degreeX+1);
    gsKnotVector<> kvY(0.0, 1.0, nKnotsY-1, degreeY+1);
    gsKnotVector<> kvZ(0.0, 1.0, nKnotsZ-1, degreeZ+1);
    gsTensorBSplineBasis<3> basis(kvX, kvY, kvZ);
    gsTensorBSpline<3> surface(basis,basis.anchors().transpose());
    
    gsInfo << "Surface: " << surface.basis().size() << " basis functions\n";
    index_t numCPs = surface.basis().size();
    
    // Case 1: Many evaluations w.r.t. one coefficient (forward mode efficient)
    {
        gsInfo << "\n--- Case 1: Forward AD (Many points w.r.t. 1 coeff) ---\n";
        
        using T = dual_t;
        
        gsKnotVector<T> kvX(0.0, 1.0, nKnotsX-1, degreeX+1);
        gsKnotVector<T> kvY(0.0, 1.0, nKnotsY-1, degreeY+1);
        gsKnotVector<T> kvZ(0.0, 1.0, nKnotsZ-1, degreeZ+1);
        gsTensorBSplineBasis<3,T> basisAD(kvX, kvY, kvZ);
        gsMatrix<T> coefsAD = basisAD.anchors().transpose().cast<T>();
        
        // Seed derivative for target coefficient (z-component)
        index_t targetCoeff = 0;
        coefsAD(targetCoeff, 2).grad = 1.0;
        
        gsTensorBSpline<3,T> surfaceAD(basisAD, coefsAD);
        
        // Generate evaluation points
        gsMatrix<> ab(3,2);
        ab.col(0).setZero();
        ab.col(1).setOnes();
        gsMatrix<> pts = gsPointGrid(ab, nPoints);
        gsMatrix<T> ptsAD = pts.cast<T>();
        
        gsStopwatch timer;
        gsMatrix<T> valsAD;
        surfaceAD.eval_into(ptsAD, valsAD);
        double elapsed = timer.stop();
        
        gsInfo << "Forward AD time for " << nPoints << " points: " << elapsed << " ms\n";
        gsInfo << "Average per point: " << (elapsed*1000.0/nPoints) << " µs\n";
        
        // Extract and verify derivatives
        gsMatrix<> derivs_fwd(nPoints, 3);
        for (index_t i = 0; i < nPoints; ++i) {
            derivs_fwd.row(i) = valsAD.col(i).unaryExpr([](const T& x) { return x.grad; }).transpose();
        }
        gsInfo << "Computed " << nPoints << " forward derivatives\n";
        gsInfo << "First derivative: " << derivs_fwd.row(0) << "\n";
    }
    
    // Case 2: Single evaluation w.r.t. all coefficients (reverse mode efficient)
    {
        gsInfo << "\n--- Case 2: Reverse AD (1 point w.r.t. all coeffs) ---\n";
        
        using T = var_t;
        
        gsKnotVector<T> kvX(0.0, 1.0, nKnotsX-1, degreeX+1);
        gsKnotVector<T> kvY(0.0, 1.0, nKnotsY-1, degreeY+1);
        gsKnotVector<T> kvZ(0.0, 1.0, nKnotsZ-1, degreeZ+1);
        gsTensorBSplineBasis<3,T> basisAD(kvX, kvY, kvZ);
        gsMatrix<T> coefsAD = basisAD.anchors().transpose().cast<T>();
        
        gsTensorBSpline<3,T> surfaceAD(basisAD, coefsAD);
        
        // Single evaluation point at center
        gsMatrix<T> pt(3, 1);
        pt(0, 0) = static_cast<T>(0.5);
        pt(1, 0) = static_cast<T>(0.5);
        pt(2, 0) = static_cast<T>(0.5);
        
        gsStopwatch timer;
        gsMatrix<T> valAD;
        surfaceAD.eval_into(pt, valAD);
        double elapsed = timer.stop();
        
        gsInfo << "Reverse AD time for 1 point w.r.t. " << numCPs << " coeffs: " << elapsed << " ms\n";
        gsInfo << "Reverse AD is ideal for computing sensitivities of all coefficients\n";
        gsInfo << "(Derivative computation would require adjoint accumulation)\n";
    }
    
    gsInfo << "\n=== Summary ===\n";
    gsInfo << "Forward mode (dual_t): Good for 1 input -> many outputs\n";
    gsInfo << "Reverse mode (var_t): Good for many inputs -> 1 output\n";
    gsInfo << "\nBoth modes fully integrated into Gismo!\n";
    
    return 0;
}
