/** @file quadrature_circle_area_example.cpp

    @brief Pure quadrature test: integrate the area of a disk immersed in
           [0,1]^2 using the three FCM quadrature rules (CutCell, Algoim,
           Octree) and compare against the exact value. No PDE is solved, so
           the reported error is purely the quadrature/geometry error of each
           rule.

    Disk:   Omega = { (x,y) : (x-0.5)^2 + (y-0.5)^2 < R^2 },  R = 0.4
    Level set:  phi(x,y) = (x-0.5)^2 + (y-0.5)^2 - R^2   ( < 0 inside )
    Exact area  = pi * R^2

    The embedding geometry is the unit square with identity map, so the
    physical area equals the sum of quadrature weights over the integration
    domain (interior + cut cells).

    Example:
    ./build/bin/quadrature_circle_area_example -r 5 -L 3

    This file is part of the G+Smo library.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <iomanip>
#include <string>
#include <vector>

using namespace gismo;

// Integrate 1 over the trimmed domain (= area of Omega) using the quadrature
// rule selected through `options`. The embedding map is the identity unit
// square, hence the physical area equals the sum of quadrature weights over
// the integration domain (interior + cut cells).
template <class T>
T integrateArea(const gsImplicitTrimmedDomain<2,T> & tr_domain,
                const gsOptionList & options)
{
    typename gsQuadRule<T>::uPtr QuRule = gsQuadrature::getPtr(tr_domain, options);

    gsMatrix<T> pts;
    gsVector<T> wts;
    T area = 0;

    for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
    {
        QuRule->mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
        for (index_t k = 0; k < wts.size(); ++k)
            area += wts[k];                 // identity map => |det J| = 1
    }
    return area;
}

int main(int argc, char *argv[])
{
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t octLevels  = 3;
    real_t  radius     = 0.4;

    gsCmdLine cmd("Quadrature test: area of an immersed disk in [0,1]^2.");
    cmd.addInt ("e", "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt ("r", "uniformRefine",   "Number of uniform refinement steps", numRefine);
    cmd.addInt ("L", "octLevels",       "Number of octree subdivision levels", octLevels);
    cmd.addReal("R", "radius",          "Disk radius", radius);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Embedding box [0,1]^2 (identity geometry map)
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());

    // Level set: negative inside the disk
    const std::string phiStr =
        "(x-0.5)^2 + (y-0.5)^2 - " + std::to_string(radius*radius);
    gsFunctionExpr<> impl_fun(phiStr, 2);

    const real_t exactArea = EIGEN_PI * radius * radius;

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);

    // The three rules under test (name + enum value)
    struct RuleInfo { const char * name; index_t id; };
    const std::vector<RuleInfo> rules = {
        { "CutCell", gsQuadrature::CutCellRule },
        { "Algoim ", gsQuadrature::AlgoimRule  },
        { "Octree ", gsQuadrature::OctreeRule  }
    };

    gsInfo << "Disk radius R   = " << radius << "\n";
    gsInfo << "Exact area      = " << std::setprecision(12) << exactArea << "\n";
    gsInfo << "Octree levels   = " << octLevels << "\n\n";

    // Store the |error| of each rule at each refinement to print convergence
    gsMatrix<real_t> err(rules.size(), numRefine + 1);

    for (index_t r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<2,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");

        gsInfo << "Refinement " << r << "\n";
        gsInfo << "  rule      area              abs. error\n";

        for (size_t i = 0; i < rules.size(); ++i)
        {
            // Options carrying the quadrature-rule selection
            gsOptionList options;
            options.addInt("quRule", "Quadrature rule id", rules[i].id);
            options.addReal("quA", "quA: nodes = quA*deg + quB", 1.0);
            options.addInt ("quB", "quB: nodes = quA*deg + quB", 1);
            options.addInt ("octLevels",
                "Number of octree subdivision levels for OctreeRule", octLevels);

            gsImplicitTrimmedDomain<2,real_t> tr_domain(impl_fun, *tbsPtr);

            const real_t area = integrateArea(tr_domain, options);
            err(i, r) = std::abs(area - exactArea);

            gsInfo << "  " << rules[i].name << "  "
                   << std::scientific << std::setprecision(8)
                   << area << "    " << err(i, r) << "\n";
        }
        gsInfo << "\n";
    }

    // Experimental orders of convergence (per rule)
    if (numRefine > 0)
    {
        gsInfo << "Experimental order of convergence (area error):\n";
        for (size_t i = 0; i < rules.size(); ++i)
        {
            gsInfo << "  " << rules[i].name << "  " << std::fixed << std::setprecision(2);
            for (index_t r = 1; r <= numRefine; ++r)
            {
                const real_t e0 = err(i, r-1), e1 = err(i, r);
                const real_t rate = (e0 > 0 && e1 > 0)
                    ? std::log(e0 / e1) / std::log(2.0) : 0.0;
                gsInfo << rate << "  ";
            }
            gsInfo << "\n";
        }
    }

    return EXIT_SUCCESS;
}
