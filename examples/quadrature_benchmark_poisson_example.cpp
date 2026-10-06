/** @file quadrature_benchmark_poisson_example.cpp

    @brief Poisson equation with Nitsche BC, sweeping all cut-cell quadrature
           methods over one geometry per run.

    Per (geometry, method, refinement): assemble, solve, report L2/H1/quadPts/timing.
    The surface (Nitsche) rule stays the same for every method; only the volume
    rule changes. Moment_qA1 is expected to produce NaN (mass matrix indefinite).

    Usage:
      ./bin/quadrature_benchmark_poisson_example --geo cow -f obj/spot.obj -r 2 --quadDepth 3
      ./bin/quadrature_benchmark_poisson_example --geo sphere -r 3 --quadDepth 3

    This file is part of the G+Smo library.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  Method table (same as measure driver, used for volume rule only)
// =============================================================================

enum MethodType { MT_CUTCELL, MT_OCTREE, MT_ALGOIM, MT_UNIFORM, MT_DIVI, MT_MOMENT };

struct MethodConfig
{
    const char*  name;
    MethodType   type;
    int          octLevels     = 0;
    int          maxDepth      = 3;
    real_t       indicatorTol  = 1e-2;
    real_t       quA           = 1;
};

static const std::vector<MethodConfig> g_methods = {
    {"CutCell",         MT_CUTCELL},
    {"Octree_L1",        MT_OCTREE,  1},
    {"Octree_L2",        MT_OCTREE,  2},
    {"Octree_L3",        MT_OCTREE,  3},
    {"Algoim",          MT_ALGOIM},
    {"Uniform_d2",      MT_UNIFORM, 0, 2},
    {"Uniform_d3",      MT_UNIFORM, 0, 3},
    {"Divi_d2_t1e-2",   MT_DIVI,    0, 2, 1e-2},
    {"Divi_d2_t1e-3",   MT_DIVI,    0, 2, 1e-3},
    {"Divi_d3_t1e-2",   MT_DIVI,    0, 3, 1e-2},
    {"Divi_d3_t1e-3",   MT_DIVI,    0, 3, 1e-3},
    {"Moment_qA1",      MT_MOMENT,  0, 3, 1e-2, 1},
    {"Moment_qA2",      MT_MOMENT,  0, 3, 1e-2, 2},
};

// =============================================================================
//  Build a volume quadrature rule for one method.
// =============================================================================

template<class T>
static typename gsQuadRule<T>::uPtr buildVolumeRule(
    MethodType type, const gsFunction<T>& phi,
    gsTensorBSplineBasis<3,T>& bkgBasis,
    const MethodConfig& mc)
{
    const index_t deg = bkgBasis.maxDegree();

    switch (type)
    {
    case MT_CUTCELL:
    {
        gsOptionList opts;
        opts.addInt  ("quRule", "Quadrature rule id", gsQuadrature::CutCellRule);
        opts.addInt  ("quB",    "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal ("quA",    "quA: nodes = quA*deg + quB", 1.0);
        gsImplicitTrimmedDomain<3,T> trd(phi, bkgBasis);
        return gsQuadrature::getPtr<T>(trd, opts);
    }
    case MT_OCTREE:
    {
        gsOptionList opts;
        opts.addInt  ("quRule",    "Quadrature rule id", gsQuadrature::OctreeRule);
        opts.addInt  ("octLevels", "Number of octree subdivision levels", mc.octLevels);
        opts.addInt  ("quB",       "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal ("quA",       "quA: nodes = quA*deg + quB", 1.0);
        gsImplicitTrimmedDomain<3,T> trd(phi, bkgBasis);
        return gsQuadrature::getPtr<T>(trd, opts);
    }
    case MT_ALGOIM:
    {
        gsOptionList opts;
        opts.addInt ("quB",  "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal("quA",  "quA: nodes = quA*deg + quB", 1.0);
        return memory::make_unique(new gsAlgoimGenericRule<T>(phi, bkgBasis, opts));
    }
    case MT_UNIFORM:
    case MT_DIVI:
    {
        gsOptionList opts = gsAlgoimAdaptiveRule<T>::defaultOptions();
        opts.setInt  ("maxDepth",     mc.maxDepth);
        opts.setInt  ("nFallback",    deg + 1);
        opts.setString("indicator",   type == MT_UNIFORM ? "uniform" : "integralChange");
        if (type == MT_DIVI)
            opts.setReal("indicatorTol", mc.indicatorTol);
        return memory::make_unique(new gsAlgoimAdaptiveRule<T>(phi, bkgBasis, opts));
    }
    case MT_MOMENT:
    {
        gsOptionList opts = gsAlgoimAdaptiveRule<T>::defaultOptions();
        opts.setInt  ("maxDepth",     mc.maxDepth);
        opts.setInt  ("nFallback",    deg + 1);
        opts.setString("indicator",   "integralChange");
        opts.setReal ("indicatorTol", mc.indicatorTol);
        opts.setReal ("quA",          mc.quA);

        // Same output order as the deleted gsAlgoimMomentFittingRule::_numNodes1d():
        // n = round(quA*maxDegree) + quB, clamped to >= 1. NOT exactnessOrder(),
        // which returns 2*deg+1 and would redefine the Moment_qA1/qA2 columns.
        const long nRaw = std::lround(opts.askReal("quA", 1.0)
                                      * static_cast<double>(deg))
                        + static_cast<long>(opts.askInt("quB", 1));
        const index_t mfOrder1d = nRaw > 0 ? static_cast<index_t>(nRaw) : 1;

        // The compressor OWNS the adaptive rule (no clone() in the hierarchy).
        return gsMomentRule<T>::make(
            typename gsQuadRule<T>::uPtr(
                new gsAlgoimAdaptiveRule<T>(phi, bkgBasis, opts)),
            mfOrder1d);
    }
    }
    GISMO_ERROR("Unknown method type");
    return nullptr;
}

// =============================================================================
//  Filter methods by csv list
// =============================================================================

static std::vector<size_t> filterMethods(const std::string& csv)
{
    std::vector<size_t> idx;
    if (csv.empty())
    {
        for (size_t i = 0; i < g_methods.size(); ++i) idx.push_back(i);
        return idx;
    }
    std::istringstream iss(csv);
    std::string tok;
    while (std::getline(iss, tok, ','))
    {
        bool found = false;
        for (size_t i = 0; i < g_methods.size(); ++i)
            if (g_methods[i].name == tok) { idx.push_back(i); found = true; break; }
        if (!found) gsWarn << "Unknown method '" << tok << "' in --methods, ignored.\n";
    }
    return idx;
}

// =============================================================================
//  main
// =============================================================================

int main(int argc, char* argv[])
{
    std::string geoKey    = "cow";
    std::string filename  = "";
    std::string outDir    = "output_poisson_benchmark";
    std::string methodsCsv= "";
    index_t     numRefine = 3;
    index_t     numElevate= 0;
    index_t     quadDepth = 3;
    std::string indicator = "integralChange";
    real_t      indicatorTol = 1e-2;
    real_t      fill      = 0.9;
    real_t      gamma     = 1e3;
    bool        plot      = false;
    bool        quadStats = false;

    gsCmdLine cmd("Poisson benchmark: sweep all cut-cell quadrature methods "
                  "for the volume rule, keeping the surface rule fixed.");
    cmd.addString("",  "geo",     "Geometry: cow, sphere",       geoKey);
    cmd.addString("f", "file",    "Input .obj mesh file (for cow)", filename);
    cmd.addInt   ("r", "refine",  "Number of uniform refinement steps", numRefine);
    cmd.addInt   ("e", "degreeElevation", "Degree elevation steps", numElevate);
    cmd.addInt   ("",  "quadDepth", "Adaptive octree depth for surface rule and "
                                    "for the underlying adaptive rule of "
                                    "moment/divi/uniform methods",  quadDepth);
    cmd.addString("",  "indicator", "Adaptive indicator: uniform, fallback, "
                                    "integralChange",               indicator);
    cmd.addReal  ("",  "indicatorTol", "integralChange tolerance", indicatorTol);
    cmd.addReal  ("",  "fill",    "Fill fraction of [0,1]^3",      fill);
    cmd.addReal  ("g", "gamma",   "Nitsche penalty parameter",      gamma);
    cmd.addString("o", "output",  "Output folder",                  outDir);
    cmd.addString("",  "methods", "CSV of method names (empty=all)", methodsCsv);
    cmd.addSwitch("quadStats", "Print quadrature-point split and CG diagnostics", quadStats);
    cmd.addSwitch("plot", "Write ParaView output", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    outDir = gsFileManager::getCanonicRepresentation(outDir);

    const std::vector<size_t> methodIdx = filterMethods(methodsCsv);
    GISMO_ENSURE(!methodIdx.empty(), "No methods selected.");

    const int D = 3;
    const index_t nFallback = 4;

    // -------------------------------------------------------------------------
    //  Geometry setup
    // -------------------------------------------------------------------------
    std::unique_ptr<gsFunction<real_t>> phiPtr;

    if (geoKey == "sphere")
    {
        real_t R = 0.4;
        std::string phiStr = "sqrt((x-0.5)^2+(y-0.5)^2+(z-0.5)^2)-" + std::to_string(R);
        phiPtr = memory::make_unique(new gsFunctionExpr<real_t>(phiStr, 3));
        gsInfo << "Geometry: sphere (signed distance, R=" << R << ", fill=" << fill << ")\n";
    }
    else
    {
        geoKey = "cow"; // default
        GISMO_ENSURE(!filename.empty(), "Geometry 'cow' requires -f/--file.");
    }

    // Mesh loading + rescaling for cow
    gsSurfMesh mesh;
    std::size_t nVert = 0;
    if (geoKey == "cow")
    {
        GISMO_ENSURE(gsReadSurfMesh(filename, mesh),
                     "Failed to read triangle mesh: " << filename);
        nVert = mesh.n_vertices();
        gsInfo << "Loaded mesh '" << filename << "': "
               << nVert << " vertices, " << mesh.n_faces() << " triangles.\n";

        const real_t extent =
            (gsSurfMeshBoundingBox(mesh).col(1) - gsSurfMeshBoundingBox(mesh).col(0)).maxCoeff();
        const real_t scale = gsNormalizeToUnitBox(mesh, fill);

        phiPtr = memory::make_unique(new gsMeshSignedDist<real_t>(mesh, gsUnitBox3()));

        gsInfo << "  fill=" << fill << ", extent=" << extent << ", scale=" << scale << "\n";
    }

    // Background geometry [0,1]^3
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    const index_t deg = dbasis.maxCwiseDegree();

    // Manufactured solution
    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
    gsFunctionExpr<> f_rhs("3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);

    gsFileManager::mkdir(outDir);

    gsInfo << "Degree: " << deg << ", quadDepth: " << quadDepth
           << ", indicator: " << indicator << ", indicatorTol: " << indicatorTol
           << ", gamma: " << gamma << "\n\n";

    // -------------------------------------------------------------------------
    //  Header
    // -------------------------------------------------------------------------
    gsInfo << std::setw(14) << "method"
           << std::setw(13) << "L2"
           << std::setw(13) << "H1"
           << std::setw(10) << "EoC(L2)"
           << std::setw(10) << "EoC(H1)"
           << std::setw(12) << "volQPs"
           << std::setw(12) << "surfQPs"
           << std::setw(12) << "cgIters"
           << std::setw(12) << "time_ms\n";

    // ParaView collections (only from last method if plot enabled)
    std::unique_ptr<gsParaviewCollection<real_t>> colCut, colAll;
    if (plot)
    {
        gsFileManager::mkdir(outDir + "/points_cut");
        gsFileManager::mkdir(outDir + "/points_all");
        colCut.reset(new gsParaviewCollection<real_t>(outDir + "/points_cut/cut"));
        colAll.reset(new gsParaviewCollection<real_t>(outDir + "/points_all/all"));
    }

    // Per-method error vectors
    gsMatrix<real_t> l2err(methodIdx.size(), numRefine + 1);
    gsMatrix<real_t> h1err(methodIdx.size(), numRefine + 1);

    // -------------------------------------------------------------------------
    //  Refinement loop
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<D,real_t>* tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<D,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        real_t hmax = 0;
        for (std::size_t p = 0; p != dbasis.nBases(); ++p)
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

        gsImplicitTrimmedDomain<D,real_t> tr_domain(*phiPtr, *tbsPtr);
        gsGaussRule<real_t> gauss(gsVector<index_t,D>::Constant(deg + 1));

        // Surface rule: same for ALL methods. Adaptive with given indicator.
        gsOptionList surfOpts = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        surfOpts.setInt   ("maxDepth",  quadDepth);
        surfOpts.setInt   ("nFallback", nFallback);
        surfOpts.setString("indicator", indicator);
        surfOpts.setReal  ("indicatorTol", indicatorTol);
        surfOpts.setInt   ("dim", D);
        gsAlgoimAdaptiveRule<real_t> surfaceRule(*phiPtr, *tbsPtr, surfOpts);

        // -----------------------------------------------------------------
        //  Sweep volume methods
        // -----------------------------------------------------------------
        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            const auto& mc = g_methods[methodIdx[mi]];
            double tMs = 0;
            index_t nVolQP = 0, nSurfQP = 0, nCG = -1;
            real_t l2 = 0, h1 = 0;
            bool failed = false;

            typename gsQuadRule<real_t>::uPtr volRule;
            gsStopwatch clock;

            try
            {
                volRule = buildVolumeRule<real_t>(mc.type, *phiPtr, *tbsPtr, mc);

                const index_t nDofs = tbsPtr->size();
                gsVector<real_t> F_vec(nDofs);
                F_vec.setZero();
                gsSparseEntries<real_t> triplets;

                // Assembly: volume
                auto addVolume = [&](const gsMatrix<real_t>& pts, const gsVector<real_t>& wts,
                                      real_t coeff = 1, bool withRhs = true)
                {
                    if (pts.cols() == 0) return;
                    gsMatrix<index_t> act;
                    gsMatrix<real_t>  bv, bd, fv;
                    tbsPtr->active_into(pts, act);
                    tbsPtr->eval_into  (pts, bv);
                    tbsPtr->deriv_into (pts, bd);
                    if (withRhs) f_rhs.eval_into(pts, fv);
                    const index_t na = act.rows();
                    for (index_t q = 0; q < pts.cols(); ++q)
                    {
                        const real_t w = wts(q) * coeff;
                        for (index_t i = 0; i < na; ++i)
                        {
                            const index_t di = act(i,q);
                            if (withRhs)
                                F_vec(di) += w * bv(i,q) * fv(0,q);
                            for (index_t j = 0; j < na; ++j)
                            {
                                real_t kij = 0;
                                for (int d = 0; d < D; ++d)
                                    kij += bd(i*D+d,q) * bd(j*D+d,q);
                                triplets.add(di, act(j,q), w * kij);
                            }
                        }
                    }
                };

                // Regularization: tiny full-background pass
                {
                    const real_t alpha = 1e-10;
                    for (auto it = tbsPtr->domain()->beginAll();
                         it != tbsPtr->domain()->endAll(); ++it)
                    {
                        gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                        gauss.mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                        addVolume(ptmp, wtmp, alpha, false);
                    }
                }

                // Volume integration: interior cells + cut cells
                for (auto it = tr_domain.beginInterior();
                     it != tr_domain.end<InteriorSign>(); ++it)
                {
                    gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                    gauss.mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                    addVolume(ptmp, wtmp);
                    nVolQP += ptmp.cols();
                }

                for (auto it = tr_domain.beginBdr(boundary::none);
                     it != tr_domain.endBdr(boundary::none); ++it)
                {
                    gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                    volRule->mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                    addVolume(ptmp, wtmp);
                    nVolQP += ptmp.cols();
                }

                // Nitsche boundary terms (same surface rule for all methods)
                surfaceRule.resetStats();
                for (auto it = tr_domain.beginBdr(boundary::none);
                     it != tr_domain.endBdr(boundary::none); ++it)
                {
                    gsMatrix<real_t> sInterior(3,0), spts(3,0);
                    gsVector<real_t> sInteriorW, swts;
                    surfaceRule.mapToSeparated(it.lowerCorner(), it.upperCorner(),
                                               sInterior, sInteriorW, spts, swts, nullptr);

                    gsMatrix<index_t> act;
                    gsMatrix<real_t>  bv, bd, phi_d, gDv;
                    tbsPtr->active_into(spts, act);
                    tbsPtr->eval_into  (spts, bv);
                    tbsPtr->deriv_into (spts, bd);
                    phiPtr->deriv_into (spts, phi_d);
                    u_exact.eval_into  (spts, gDv);
                    const index_t na = act.rows();
                    nSurfQP += spts.cols();

                    for (index_t q = 0; q < spts.cols(); ++q)
                    {
                        const real_t w = swts(q);
                        real_t nNorm2 = 0;
                        for (int d = 0; d < D; ++d)
                            nNorm2 += phi_d(d,q) * phi_d(d,q);
                        const real_t nNorm = std::sqrt(nNorm2);
                        if (nNorm < 1e-14) continue;

                        for (index_t i = 0; i < na; ++i)
                        {
                            const index_t di = act(i,q);
                            const real_t  vi = bv(i,q);
                            real_t gi_n = 0;
                            for (int d = 0; d < D; ++d)
                                gi_n += bd(i*D+d,q) * phi_d(d,q) / nNorm;

                            F_vec(di) += w * (-gi_n * gDv(0,q) + (gamma/hmax) * vi * gDv(0,q));

                            for (index_t j = 0; j < na; ++j)
                            {
                                const index_t dj = act(j,q);
                                const real_t  vj = bv(j,q);
                                real_t gj_n = 0;
                                for (int d = 0; d < D; ++d)
                                    gj_n += bd(j*D+d,q) * phi_d(d,q) / nNorm;

                                const real_t kij = -gi_n*vj - vi*gj_n + (gamma/hmax)*vi*vj;
                                triplets.add(di, dj, w * kij);
                            }
                        }
                    }
                }

                // Solve
                gsSparseMatrix<real_t> K(nDofs, nDofs);
                K.setFrom(triplets);

                gsVector<real_t> solVector(nDofs);
                bool solveOk = true;
                real_t cgError = -1;
                {
#                   ifdef GISMO_WITH_PARDISO
                    gsSparseSolver<real_t>::PardisoLDLT slvr;
                    slvr.compute(K);
                    solVector = slvr.solve(F_vec);
#                   else
                    gsSparseSolver<real_t>::CGDiagonal slvr;
                    slvr.compute(K);
                    solVector = slvr.solve(F_vec);
                    nCG = static_cast<index_t>(slvr.iterations());
                    cgError = slvr.error();
                    solveOk = slvr.succeed();
#                   endif
                }
                tMs = clock.stop() * 1000.0;

                // Error norms using the same quadrature
                real_t l2sq = 0, h1sq = 0;
                auto addError = [&](const gsMatrix<real_t>& pts, const gsVector<real_t>& wts)
                {
                    if (pts.cols() == 0) return;
                    gsMatrix<index_t> act;
                    gsMatrix<real_t>  bv, bd, uExV, uExD;
                    tbsPtr->active_into(pts, act);
                    tbsPtr->eval_into  (pts, bv);
                    tbsPtr->deriv_into (pts, bd);
                    u_exact.eval_into  (pts, uExV);
                    u_exact.deriv_into (pts, uExD);
                    const index_t na = act.rows();
                    for (index_t q = 0; q < pts.cols(); ++q)
                    {
                        const real_t w = wts(q);
                        real_t uh = 0;
                        for (index_t i = 0; i < na; ++i)
                            uh += bv(i,q) * solVector(act(i,q));
                        const real_t dl2 = uh - uExV(0,q);
                        l2sq += w * dl2 * dl2;
                        for (int d = 0; d < D; ++d)
                        {
                            real_t guh_d = 0;
                            for (index_t i = 0; i < na; ++i)
                                guh_d += bd(i*D+d,q) * solVector(act(i,q));
                            const real_t dh1 = guh_d - uExD(d,q);
                            h1sq += w * dh1 * dh1;
                        }
                    }
                };

                for (auto it = tr_domain.beginInterior();
                     it != tr_domain.end<InteriorSign>(); ++it)
                {
                    gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                    gauss.mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                    addError(ptmp, wtmp);
                }
                for (auto it = tr_domain.beginBdr(boundary::none);
                     it != tr_domain.endBdr(boundary::none); ++it)
                {
                    gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                    volRule->mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                    addError(ptmp, wtmp);
                }

                l2 = math::sqrt(l2sq);
                h1 = l2 + math::sqrt(h1sq);

                // Check for NaN in error (indicative of junk solve)
                if (!std::isfinite(l2) || !std::isfinite(h1))
                    failed = true;
            }
            catch (std::exception& e)
            {
                failed = true;
                gsInfo << "  " << std::setw(12) << mc.name << " FAILED: " << e.what();
            }

            l2err(mi, r) = l2;
            h1err(mi, r) = h1;

            gsInfo << "  " << std::setw(12) << mc.name;
            if (failed)
            {
                gsInfo << std::setw(13) << "nan"
                       << std::setw(13) << "nan"
                       << std::setw(10) << "-"
                       << std::setw(10) << "-"
                       << std::setw(12) << nVolQP
                       << std::setw(12) << nSurfQP
                       << std::setw(12) << "-"
                       << std::setw(12) << std::fixed << std::setprecision(2) << tMs
                       << "  (FAILED)\n";
                continue;
            }

            // EoC requires previous refinement; print N/A at r==0
            gsInfo << std::setw(13) << std::scientific << std::setprecision(3) << l2
                   << std::setw(13) << std::scientific << std::setprecision(3) << h1;
            if (r > 0)
            {
                const real_t l0 = l2err(mi, r-1), l1 = l2;
                const real_t h0 = h1err(mi, r-1), h1v = h1;
                const real_t eocL2 = (l0 > 0 && l1 > 0) ? std::log(l0/l1)/std::log(2.0) : 0;
                const real_t eocH1 = (h0 > 0 && h1v > 0) ? std::log(h0/h1v)/std::log(2.0) : 0;
                gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << eocL2
                       << std::setw(10) << std::fixed << std::setprecision(2) << eocH1;
            }
            else
            {
                gsInfo << std::setw(10) << "N/A"
                       << std::setw(10) << "N/A";
            }
            gsInfo << std::setw(12) << nVolQP
                   << std::setw(12) << nSurfQP
                   << std::setw(12) << nCG
                   << std::setw(12) << std::fixed << std::setprecision(2) << tMs << "\n";

            gsInfo << std::flush;

            // Plot from last method only
            if (plot && mi == methodIdx.size() - 1 && volRule)
            {
                const std::string rs = std::to_string(r);
                gsMatrix<real_t> quadCut(4,0), quadAll(4,0);
                auto appendCloud = [&](gsMatrix<real_t>& acc,
                                       const gsMatrix<real_t>& phys,
                                       const gsVector<real_t>& scalar)
                {
                    const index_t c = acc.cols();
                    acc.conservativeResize(4, c + phys.cols());
                    acc.block(0, c, 3, phys.cols()) = phys;
                    acc.row(3).segment(c, phys.cols()) = scalar.transpose();
                };

                for (auto it = tr_domain.beginBdr(boundary::none);
                     it != tr_domain.endBdr(boundary::none); ++it)
                {
                    gsMatrix<real_t> ptmp; gsVector<real_t> wtmp;
                    volRule->mapTo(it.lowerCorner(), it.upperCorner(), ptmp, wtmp);
                    if (ptmp.cols() > 0)
                    {
                        gsMatrix<real_t> phys;
                        mp.patch(0).eval_into(ptmp, phys);
                        appendCloud(quadCut, phys, wtmp);
                        appendCloud(quadAll, phys, wtmp);
                    }
                }
                if (quadCut.cols() > 0)
                {
                    gsWriteParaviewPoints(quadCut, outDir + "/points_cut/cut_r" + rs);
                    colCut->addPart("cut_r" + rs + ".vtp", r);
                }
                if (quadAll.cols() > 0)
                {
                    gsWriteParaviewPoints(quadAll, outDir + "/points_all/all_r" + rs);
                    colAll->addPart("all_r" + rs + ".vtp", r);
                }
            }
        }

        gsInfo << "(h=" << hmax << ", gamma/h=" << gamma/hmax << ")\n";
    }

    // -------------------------------------------------------------------------
    //  Full EoC table
    // -------------------------------------------------------------------------
    if (numRefine > 0)
    {
        gsInfo << "\n\n--- Convergence summary ---\n";
        gsInfo << "L2 error history:\n";
        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            gsInfo << "  " << std::setw(12) << g_methods[methodIdx[mi]].name
                   << "  " << std::scientific << std::setprecision(3)
                   << l2err.row(mi).transpose() << "\n";
        }
        gsInfo << "\nH1 error history:\n";
        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            gsInfo << "  " << std::setw(12) << g_methods[methodIdx[mi]].name
                   << "  " << std::scientific << std::setprecision(3)
                   << h1err.row(mi).transpose() << "\n";
        }

        gsInfo << "\nEoC (L2):\n" << std::setw(14) << "method";
        for (int r = 1; r <= numRefine; ++r)
        {
            std::ostringstream oss;
            oss << "r" << (r-1) << "->" << r;
            gsInfo << std::setw(10) << oss.str();
        }
        gsInfo << "\n";
        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            gsInfo << "  " << std::setw(12) << g_methods[methodIdx[mi]].name;
            for (int r = 1; r <= numRefine; ++r)
            {
                const real_t e0 = l2err(mi, r-1), e1 = l2err(mi, r);
                if (e0 > 0 && e1 > 0 && std::isfinite(e0) && std::isfinite(e1))
                {
                    const real_t rate = std::log(e0/e1)/std::log(2.0);
                    gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << rate;
                }
                else
                    gsInfo << std::setw(10) << "  -";
            }
            gsInfo << "\n";
        }

        gsInfo << "\nEoC (H1):\n" << std::setw(14) << "method";
        for (int r = 1; r <= numRefine; ++r)
        {
            std::ostringstream oss;
            oss << "r" << (r-1) << "->" << r;
            gsInfo << std::setw(10) << oss.str();
        }
        gsInfo << "\n";
        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            gsInfo << "  " << std::setw(12) << g_methods[methodIdx[mi]].name;
            for (int r = 1; r <= numRefine; ++r)
            {
                const real_t e0 = h1err(mi, r-1), e1 = h1err(mi, r);
                if (e0 > 0 && e1 > 0 && std::isfinite(e0) && std::isfinite(e1))
                {
                    const real_t rate = std::log(e0/e1)/std::log(2.0);
                    gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << rate;
                }
                else
                    gsInfo << std::setw(10) << "  -";
            }
            gsInfo << "\n";
        }
    }

    if (plot)
    {
        if (geoKey == "cow")
        {
            gsWriteParaview(mesh, outDir + "/geometry");
        }
        if (colCut) colCut->save();
        if (colAll) colAll->save();
        gsInfo << "\nParaView files written to: " << outDir << "\n";
    }

    return EXIT_SUCCESS;
}
