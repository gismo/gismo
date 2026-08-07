#include <iostream>
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>

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

    gsInfo << "=== B-spline Basis: Knot Position Derivatives ===\n";
    gsInfo << "Computing derivatives of basis values w.r.t. knot positions using AD\n\n";

    // Create base knot vector
    typedef dual_t T;
    typedef GISMO_COEFF_TYPE double_t;
    std::vector<T> knots({0.0,0.0,0.0,0.25,0.5,0.75,1.0,1.0,1.0});

    // Generate evaluation points
    gsMatrix<T> pts(1, nPoints);
    pts.row(0).setLinSpaced(nPoints,0.0,1.0);
    gsInfo<<"Points: "<<pts<<"\n";
    gsMatrix<T> vals;
    std::vector<gsMatrix<double_t>> ders(knots.size());
    gsMatrix<index_t> actives;
    for (index_t k = degree+1; k != knots.size()-degree-1; k++)
    {
        gsInfo<<"Derivative w.r.t. knot "<<k<<" with value "<<knots[k]<<"\n";
        autodiff::detail::seed<1>(knots[k], 1.0);

        gsKnotVector<T> kv(knots,2 );
        gsBSplineBasis<T> basis(kv);

        basis.active_into(pts,actives);
        basis.eval_into(pts,vals);

        ders[k].resize(basis.size(),nPoints);
        for (index_t i=0; i!=actives.rows(); i++)
            for (index_t j=0; j!= actives.cols(); j++)
                ders[k](actives(i,j),j) = autodiff::detail::derivative<1>(vals(i,j));

        gsInfo<<"Derivative:\n"<<ders[k]<<"\n";
    }
    return 0;
}
