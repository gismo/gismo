/** @file reproduce_surface_metrics.cpp

    @brief Reproduces the surface-parameterization quality tables (scaled
           Jacobian and uniformity) of the composite-spline-relocation paper,
           i.e. Tables 1 (human face) and 2 (ogre face).

    For a given input surface it re-parameterises with the indicator-free
    monitor (paper Eq. 14) using composite splines of bi-degree (d,d) with
    N x N elements, and evaluates on a dense grid the metrics

        m_SJ   = || G_xi x G_eta || / ( ||G_xi|| ||G_eta|| )   (scaled Jacobian)
        m_unif = ( area / A - 1 )^2 ,   A = mean(area)         (uniformity)

    where G = S o sigma is the composed geometry.  It prints, for the initial
    geometry and for each (d, N) configuration, min/avg m_SJ and max/avg m_unif,
    matching the columns of Tables 1 and 2.

    Usage:
        reproduce_surface_metrics -i <surface.xml> [-g <grid>] [--csv out.csv]

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

using namespace gismo;

// Evaluate the paper metrics of a (composed) surface map on an (ngrid x ngrid)
// uniform grid of the unit parametric square.
template<class Geo>
void metrics(const Geo & G, index_t ngrid,
             real_t & minSJ, real_t & avgSJ, real_t & maxU, real_t & avgU)
{
    gsVector<index_t> np(2); np << ngrid, ngrid;
    gsMatrix<> pts = gsPointGrid<real_t>(gsVector<real_t>::Zero(2),
                                         gsVector<real_t>::Ones(2), np);
    gsMatrix<> D; G.deriv_into(pts, D); // rows: [G0_xi,G0_eta, G1_xi,G1_eta, G2_xi,G2_eta]
    const index_t N = pts.cols();
    gsVector<> area(N), sj(N);
    for (index_t k = 0; k != N; ++k)
    {
        gsVector<real_t,3> gx, gy;
        gx << D(0,k), D(2,k), D(4,k);
        gy << D(1,k), D(3,k), D(5,k);
        const real_t a = gx.cross(gy).norm();
        area(k) = a;
        const real_t nx = gx.norm(), ny = gy.norm();
        sj(k) = (nx*ny > 1e-14) ? a/(nx*ny) : 0.0;
    }
    const real_t A = area.mean();          // parametric area == 1, so A == total area
    minSJ = sj.minCoeff();  avgSJ = sj.mean();
    gsVector<> U(N);
    for (index_t k = 0; k != N; ++k) { const real_t r = area(k)/A - 1.0; U(k) = r*r; }
    maxU = U.maxCoeff();  avgU = U.mean();
}

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_face.xml";
    std::string csv;
    index_t ngrid = 1001;         // paper uses 1001 x 1001
    // Penalty is BOTH the fold-barrier epsilon and the Tikhonov shift. Too small a value
    // makes the barrier so stiff that L-BFGS' first trial step lands in a region where the
    // objective is astronomically large, the line search cannot bracket, and it silently
    // returns the identity map (observed for the ogre at d=2, 4x4/8x8 with 1e-4).
    // 1e-2 is the library default and converges for every configuration.
    real_t penalty = 1e-2, smoothing = 1e1;
    index_t maxIt = 250;          // cap L-BFGS iterations (paper: "a few dozen")

    gsCmdLine cmd("Reproduce the surface-metric Tables 1/2 of the paper.");
    cmd.addString("i","input","Input file with a 'geometry' surface (or first geometry)",input);
    cmd.addString("","csv","Optional CSV output file",csv);
    cmd.addInt("g","grid","Sampling points per direction",ngrid);
    cmd.addInt("","maxIt","Maximum L-BFGS iterations",maxIt);
    cmd.addReal("P","penalty","Fold-barrier penalty",penalty);
    cmd.addReal("S","smoothing","Smoothing (unused without indicator)",smoothing);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Fail loudly on a missing input: the default path is relative to the CWD, and a
    // silent empty run is far worse than an error message.
    GISMO_ENSURE(gsFileManager::fileExists(input),
                 "Input file not found: '"<<input<<"'.\n"
                 "  The default path is relative to the current directory; pass an "
                 "absolute path with -i <surface.xml>.");

    gsFileData<> fd(input);
    gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    const gsGeometry<>& S = geomMP.patch(0);
    GISMO_ENSURE(S.domainDim()==2 && S.targetDim()==3, "Expected a surface (2->3).");
    gsInfo << "Geometry: "<<S.basis().numElements()<<" elements, grid "<<ngrid<<"x"<<ngrid<<"\n\n";

    gsInfo << std::left << std::setw(12) << "config"
           << std::right << std::setw(12) << "minSJ" << std::setw(12) << "avgSJ"
           << std::setw(12) << "maxUnif" << std::setw(12) << "avgUnif"
           << std::setw(11) << "minDetJs" << "\n";

    // Initial (identity composition = the geometry itself)
    {
        real_t a,b,c,d; metrics(S, ngrid, a,b,c,d);
        gsInfo << std::left << std::setw(12) << "Initial" << std::right
               << std::setw(12) << a << std::setw(12) << b
               << std::setw(12) << c << std::setw(12) << d << std::setw(11) << "-" << "\n" << std::flush;
    }

    struct Cfg { index_t d, N; };
    const std::vector<Cfg> cfgs = { {1,2},{1,4},{1,8},
                                    {2,1},{2,2},{2,4},{2,8},
                                    {3,1},{3,2},{3,4},{3,8} };
    std::ofstream csvfile;
    if (!csv.empty()) { csvfile.open(csv); csvfile<<"d,N,minSJ,avgSJ,maxUnif,avgUnif,minDetJs\n"; }

    for (const Cfg & cfg : cfgs)
    {
        // Identity composition: degree-d unit square, refined to N x N elements.
        gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(cfg.d);
        const index_t nref = (index_t)std::round(std::log2((double)cfg.N));
        for (index_t r = 0; r < nref; ++r) comp->uniformRefine();
        gsSquareDomain<real_t> domain(*comp);
        domain.options().addSwitch("Slide","",true);   // boundary slides (paper Fig. 2)
        domain.applyOptions();
        std::ostringstream label; label<<"d="<<cfg.d<<" "<<cfg.N<<"x"<<cfg.N;
        if (domain.nControls()==0)
        {
            gsWarn<<"["<<label.str()<<"] skipped: no free control points.\n";
            continue;                                   // no interior DoF (e.g. d=1, 1x1)
        }

        gsOptim<real_t>::LBFGS opt;
        opt.options().setInt("MaxIterations",maxIt);
        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relo(domain, S, opt, /*parametric*/false);
        relo.options().setReal("Penalty",penalty);
        relo.options().setReal("Smoothing",smoothing);

        // L-BFGS can fail its very first line search and return the INITIAL design
        // (the identity sigma) without any error. That looks like a converged run whose
        // metrics happen to equal the Initial row -- a silently wrong result. Detect it.
        const gsVector<real_t> ctrlBefore = domain.getControls();
        relo.solve();
        const gsVector<real_t> ctrlAfter  = domain.getControls();
        const bool moved = (ctrlAfter-ctrlBefore).norm() > 1e-12;

        const real_t minJs = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(minJs) && minJs > 0,
                     "["<<label.str()<<"] re-parameterisation folded/degenerate: "
                     "min det Jsigma = "<<minJs);
        if (!moved)
            gsWarn<<"["<<label.str()<<"] L-BFGS did not move (line search failed at iteration 0); "
                    "the metrics below are just the Initial ones. Try a larger --penalty (-P).\n";

        gsComposedGeometry<real_t> G(domain.domain(), S);
        real_t a,b,c,d; metrics(G, ngrid, a,b,c,d);

        gsInfo << std::left << std::setw(12) << label.str() << std::right
               << std::setw(12) << a << std::setw(12) << b
               << std::setw(12) << c << std::setw(12) << d << std::setw(11) << minJs << "\n" << std::flush;
        if (csvfile.is_open())
            csvfile<<cfg.d<<","<<cfg.N<<","<<a<<","<<b<<","<<c<<","<<d<<","<<minJs<<"\n";
    }

    gsInfo << "\nLegend: m_SJ = scaled Jacobian (orthogonality, max 1), "
           << "m_unif = uniformity (min 0). Compare against paper Tables 1/2.\n";
    return 0;
}
