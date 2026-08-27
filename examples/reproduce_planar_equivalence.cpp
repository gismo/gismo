/** @file reproduce_planar_equivalence.cpp

    @brief Reproduces Section 4.1 (planar parameterization) of the
           composite-spline-relocation paper, i.e. the four L2 distances
           reported in the sub-captions of Figure 5.

    The same trapezoidal geometry is re-parameterised twice with the
    gradient-based monitor (paper Eq. 12, smoothing theta = 10) and the star
    indicator (paper Eq. 33):

      * PLANAR  formulation: the geometry is a 2 -> 2 map, so the area element
        is the Jacobian determinant det(J)          (paper Eq. planar_area);
      * SURFACE formulation: the SAME trapezoid embedded in 3D (z = 0), so the
        area element is || d1 G x d2 G ||           (paper Eq. area_element).

    Both formulations should yield (nearly) the same deformation sigma. The
    paper quantifies this with the average L2 distance between the two control
    nets of sigma, reported per configuration in Figure 5:

        d=1,  8x8   ->  1.51e-10
        d=1, 16x16  ->  6.88e-03
        d=2,  8x8   ->  8.15e-03
        d=2, 16x16  ->  1.08e-02

    This driver recomputes exactly those four numbers.

    Usage:
        reproduce_planar_equivalence [-p <planar.xml>] [-s <surface.xml>]

    Both inputs ship with the paper's data set as
        parametrization/monitor_example_planar.xml          (2D geometry, 2D indicator)
        parametrization/monitor_example_planar_surface.xml  (3D geometry, 3D indicator)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji (y.ji-1@tudelft.nl)
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <iomanip>

using namespace gismo;

// Read the 'geometry' patch and the 'function' indicator from one of the paper's xml files.
static void readCase(const std::string & file,
                     gsMultiPatch<real_t> & mp,
                     gsFunctionExpr<real_t> & fun)
{
    GISMO_ENSURE(gsFileManager::fileExists(file),
                 "Input file not found: '"<<file<<"'.\n"
                 "  The default path is relative to the current directory; pass an "
                 "absolute path with -p / -s.");
    gsFileData<real_t> fd(file);
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",mp);
    else                         mp.addPatch(*fd.getFirst<gsGeometry<real_t>>());
    GISMO_ENSURE(fd.template has<gsFunctionExpr<real_t>>(),
                 "No indicator function in "<<file);
    fd.getFirst<gsFunctionExpr<real_t>>(fun);
}

// Re-parameterise `geo` with the gradient-based monitor and return sigma's control vector.
static gsVector<real_t> relocate(const gsGeometry<real_t> & geo,
                                 const gsFunctionExpr<real_t> & fun,
                                 index_t deg, index_t N,
                                 real_t theta, real_t penalty, index_t maxIt)
{
    gsGeometry<real_t>::uPtr comp = gsNurbsCreator<real_t>::BSplineSquareDeg(deg);
    const index_t nref = (index_t)std::round(std::log2((double)N));
    GISMO_ENSURE((1<<nref)==N, "N must be a power of two, got "<<N);
    for (index_t r=0;r<nref;++r) comp->uniformRefine();

    gsSquareDomain<real_t> domain(*comp);
    domain.options().addSwitch("Slide","",true);   // boundary slides (paper Fig. 2)
    domain.applyOptions();

    gsOptim<real_t>::LBFGS opt;
    opt.options().setInt("MaxIterations",maxIt);
    gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>
        relo(domain, geo, fun, opt, /*parametric*/false);
    relo.options().setReal("Smoothing",theta);
    relo.options().setReal("Penalty",penalty);

    const gsVector<real_t> before = domain.getControls();
    relo.solve();
    const gsVector<real_t> after = domain.getControls();

    const real_t minJs = relo.computeMinJacobian();
    GISMO_ENSURE(math::isfinite(minJs) && minJs > 0,
                 "Re-parameterisation folded/degenerate (min det Jsigma = "<<minJs<<").");
    GISMO_ENSURE((after-before).norm() > 1e-12,
                 "L-BFGS did not move (line search failed at iteration 0). Try a larger --penalty.");
    return after;
}

int main(int argc, char** argv)
{
    std::string planarFile  = "parametrization/monitor_example_planar.xml";
    std::string surfaceFile = "parametrization/monitor_example_planar_surface.xml";
    real_t theta   = 1e1;    // paper Section 4.1: theta = 10
    real_t penalty = 1e-2;
    index_t maxIt  = 250;

    gsCmdLine cmd("Reproduce Section 4.1 (planar vs surface formulation, Figure 5).");
    cmd.addString("p","planar" ,"Planar (2->2) trapezoid xml" ,planarFile);
    cmd.addString("s","surface","Surface (2->3) trapezoid xml",surfaceFile);
    cmd.addReal ("S","smoothing","Smoothing parameter theta",theta);
    cmd.addReal ("P","penalty","Fold-barrier penalty",penalty);
    cmd.addInt  ("","maxIt","Maximum L-BFGS iterations",maxIt);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> mpP, mpS;
    gsFunctionExpr<real_t> funP, funS;
    readCase(planarFile , mpP, funP);
    readCase(surfaceFile, mpS, funS);

    const gsGeometry<real_t> & gP = mpP.patch(0);
    const gsGeometry<real_t> & gS = mpS.patch(0);
    GISMO_ENSURE(gP.domainDim()==2 && gP.targetDim()==2, "Planar input must be a 2->2 map.");
    GISMO_ENSURE(gS.domainDim()==2 && gS.targetDim()==3, "Surface input must be a 2->3 map.");

    gsInfo<<"Section 4.1: planar vs surface formulation, theta = "<<theta<<"\n";
    gsInfo<<"Paper Figure 5 sub-captions give the average L2 distance between the two\n"
            "control nets of sigma. Reproducing them:\n\n";

    gsInfo<<std::left<<std::setw(14)<<"config"
          <<std::right<<std::setw(16)<<"avg L2 dist"
          <<std::setw(16)<<"paper"
          <<std::setw(14)<<"rel. diff"<<"\n";
    gsInfo<<std::string(60,'-')<<"\n";

    struct Cfg { index_t d, N; real_t paper; };
    const std::vector<Cfg> cfgs = { {1, 8, 1.51e-10}, {1,16, 6.88e-3},
                                    {2, 8, 8.15e-3 }, {2,16, 1.08e-2} };

    for (const Cfg & c : cfgs)
    {
        const gsVector<real_t> cP = relocate(gP, funP, c.d, c.N, theta, penalty, maxIt);
        const gsVector<real_t> cS = relocate(gS, funS, c.d, c.N, theta, penalty, maxIt);
        GISMO_ENSURE(cP.size()==cS.size(), "Control nets have different sizes.");

        // sigma's controls are stored as [x-coords ; y-coords] of the free control points;
        // the "average L2 distance between the control nets" is the mean over control
        // points of the Euclidean distance between the planar and surface nets.
        const index_t nc = cP.size()/2;
        real_t acc = 0;
        for (index_t i=0;i<nc;++i)
        {
            const real_t dx = cP[i]      - cS[i];
            const real_t dy = cP[i+nc]   - cS[i+nc];
            acc += math::sqrt(dx*dx + dy*dy);
        }
        const real_t avgL2 = acc / (real_t)nc;

        std::ostringstream lab; lab<<"d="<<c.d<<" "<<c.N<<"x"<<c.N;
        const real_t rel = math::abs(avgL2 - c.paper) / math::max(c.paper,(real_t)1e-30);
        gsInfo<<std::left<<std::setw(14)<<lab.str()
              <<std::right<<std::scientific<<std::setprecision(3)
              <<std::setw(16)<<avgL2
              <<std::setw(16)<<c.paper
              <<std::setw(14)<<rel<<"\n";
    }

    gsInfo<<"\nBoth formulations should give nearly identical sigma: a small L2 distance\n"
            "confirms the planar and surface formulations are numerically equivalent.\n";
    return 0;
}
