/** @file bSplineBasis_knot_derivative.cpp

    @brief Example demonstrating automatic differentiation of B-spline basis evaluations w.r.t. knot positions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#include <iostream>
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    index_t nPoints = 10;
    index_t degree = 2;
    index_t nKnots = 3;

    gsCmdLine cmd("Compute derivatives of B-spline basis evaluations w.r.t. knot positions.");
    cmd.addInt("N", "npoints", "Number of evaluation points", nPoints);
    cmd.addInt("p", "degree", "Degree of B-spline", degree);
    cmd.addInt("k", "nknots", "Number of interior knots", nKnots);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

#ifdef GISMO_AUTODIFF_FORWARD
    gsInfo << "=== B-spline Basis: Knot Position Derivatives ===\n";
    gsInfo << "Forward AD with a finite-difference cross-check\n\n";

    // Create base knot vector (quadratic example)
    typedef dual_t T;
    typedef GISMO_COEFF_TYPE double_t;
    const index_t basisDegree = 2;
    std::vector<T> knots({0.0,0.0,0.0,0.25,0.5,0.75,1.0,1.0,1.0});

    // Generate evaluation points
    gsMatrix<T> pts(1, nPoints);
    pts.row(0).setLinSpaced(nPoints,0.0,1.0);
    gsInfo<<"Points: "<<pts<<"\n";
    gsMatrix<T> vals;
    std::vector<gsMatrix<double_t>> ders(knots.size());
    gsMatrix<index_t> actives;
    for (size_t k = basisDegree+1; k != knots.size()-basisDegree-1; k++)
    {
        std::vector<T> knots_ad(knots.size());
        for (size_t i = 0; i != knots.size(); ++i)
            knots_ad[i] = autodiff::val(knots[i]);

        gsInfo<<"Derivative w.r.t. knot "<<k<<" with value "<<knots[k]<<"\n";
        autodiff::detail::seed<1>(knots_ad[k], 1.0);

        gsKnotVector<T> kv(knots_ad,basisDegree );
        gsBSplineBasis<T> basis(kv);

        basis.active_into(pts,actives);
        basis.eval_into(pts,vals);

        ders[k].resize(basis.size(),nPoints);
        ders[k].setZero();
        for (index_t i=0; i!=actives.rows(); i++)
            for (index_t j=0; j!= actives.cols(); j++)
                ders[k](actives(i,j),j) = autodiff::detail::derivative<1>(vals(i,j));

        gsInfo<<"Derivative:\n"<<ders[k]<<"\n";
    }

    // ===== Finite difference verification =====
    {
        real_t h = 1e-6;
        gsInfo << "\n=== Finite Difference Verification ===\n";
        gsInfo << "Central differences with step h = " << h << "\n";

        // Base knot vector and basis in real arithmetic
        std::vector<real_t> knots_real(knots.size());
        for (size_t k = 0; k != knots.size(); ++k)
            knots_real[k] = autodiff::val(knots[k]);

        gsKnotVector<real_t> kv_base(knots_real, basisDegree);
        gsBSplineBasis<real_t> basis_base(kv_base);

        gsMatrix<real_t> pts_real(1, nPoints);
        pts_real.row(0).setLinSpaced(nPoints, 0.0, 1.0);

        real_t max_error = 0.0;
        for (size_t k = basisDegree+1; k != knots.size()-basisDegree-1; ++k)
        {
            std::vector<real_t> knots_plus  = knots_real;
            std::vector<real_t> knots_minus = knots_real;
            knots_plus[k]  += h;
            knots_minus[k] -= h;

            gsMatrix<real_t> ders_fd(basis_base.size(), nPoints);
            ders_fd.setZero();
            gsBSplineBasis<real_t> basis_plus (gsKnotVector<real_t>(knots_plus,  basisDegree));
            gsBSplineBasis<real_t> basis_minus(gsKnotVector<real_t>(knots_minus, basisDegree));
            for (index_t i = 0; i != basis_base.size(); ++i)
            {
                const gsMatrix<real_t> vals_plus_i  = basis_plus.evalSingle(i, pts_real);
                const gsMatrix<real_t> vals_minus_i = basis_minus.evalSingle(i, pts_real);
                ders_fd.row(i) = (vals_plus_i.row(0) - vals_minus_i.row(0)) / (2.0 * h);
            }

            real_t err = (ders_fd - ders[k].template cast<real_t>()).cwiseAbs().maxCoeff();
            max_error = std::max(max_error, err);
            gsInfo << "Max |AD - FD| for knot " << k << ": " << err << "\n";
        }

        gsInfo << "Global max error: " << max_error << "\n";
        if (max_error < 1e-4)
            gsInfo << "AD derivatives verified against finite differences!\n";
        else
            gsInfo << "WARNING: AD and FD derivatives deviate significantly.\n";
    }
#else
    gsInfo << "=== B-spline Basis: Knot Position Derivatives ===\n";
    gsInfo << "This example requires GISMO_AUTODIFF_FORWARD.\n";
#endif
    return 0;
}
