/** @file poisson2_nitsche_immersed_geometry_octree_adaptive_example.cpp

    @brief Poisson equation on an immersed 3D domain (a triangle mesh loaded
           from a Wavefront .obj file) inside [0,1]^3, solved with the Finite
           Cell Method, using Algoim combined with an ADAPTIVE octree
           subdivision of the cut cells for the immersed quadrature.

    This is the cost-optimal variant of
    poisson2_nitsche_immersed_geometry_octree_example.cpp. Instead of splitting
    every cut cell uniformly to --quadDepth, the octree only refines where the
    implicit interface actually is:

      For each background cut cell (from gsImplicitTrimmedDomain<3>), the box is
      recursively subdivided. Each child box is classified with a Lipschitz
      signed-distance test (phi is ~1-Lipschitz):
        phi(center) < -radius  => child fully INSIDE  {phi<0}  -> cheap Gauss
        phi(center) >  radius  => child fully OUTSIDE          -> skipped
        otherwise              => child still CUT              -> subdivide
      where radius = half the box diagonal. Cut children are subdivided until
      --quadDepth; a child that is still cut at that depth gets Algoim (+ the
      uniform-grid fallback).

    Consequently the expensive Algoim work is spent only in the narrow band
    around {phi=0}, while the fully-interior sub-cells use a plain Gauss rule
    and the fully-exterior ones are dropped. quadDepth==0 reproduces the
    original behaviour (Algoim once on the whole cut cell).

    In the ParaView output the Gauss points from inside sub-cells are tagged as
    "interior" and the Algoim points as "cut" accordingly. The total number of
    volume quadrature points is printed per refinement level.

    This is adaptive QUADRATURE only. The B-spline basis, the trimmed-domain
    classification, the solver, and the Poisson problem are unchanged.

    Optionally (--momentFit) the cut-cell VOLUME rule is replaced by the
    moment-fitting rule: the adaptive point cloud of every cut cell is
    compressed onto the fixed tensor Gauss grid of that cell, so the number of
    volume quadrature points per cut cell no longer depends on the octree
    depth. The SURFACE (Nitsche) rule is never compressed -- its points live on
    {phi==0} and cannot be represented on the element grid -- so all Nitsche
    terms are bit-identical to the non-momentFit run.

    Order caveat (important). Moment fitting reproduces exactly what the nodal
    Lagrange basis of the output grid interpolates exactly, i.e. tensor degree
    n-1 per direction, whereas an n-point Gauss rule is exact to degree 2n-1.
    The output grid has n = quA*p + quB points per direction (quB = 1 here), so
    the inherited default quA = 1 gives n = p+1 and reproduces only degree p --
    NOT enough for the stiffness integrand grad(u).grad(v), whose tensor degree
    reaches 2p. Use --quA 2 (n = 2p+1) for an adequate operator. --quA is
    applied to the moment-fitting option list only, so it also raises the order
    of the underlying Algoim rule (the rule has a single quA key for both);
    with --momentFit off it has no effect at all.

    Diagnostics: --quadStats prints the interior/cut/surface quadrature-point
    split and the CG iteration count and residual (this build has no PARDISO,
    so the solve is CGDiagonal, which does NOT throw on an indefinite matrix).
    --errRefRule computes the error norms with the uncompressed adaptive rule
    while keeping the compressed rule in the assembly, which separates "the
    solution got worse" from "the error functional was mis-integrated".

    Usage:
    ./bin/poisson2_nitsche_immersed_geometry_octree_adaptive_example -f obj/spot.obj -r 2 -g 1e3 --quadDepth 2 --plot
    ./bin/poisson2_nitsche_immersed_geometry_octree_adaptive_example -f obj/spot.obj -r 2 --quadDepth 3 --indicator uniform --momentFit --quA 2

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  main
// =============================================================================
int main(int argc, char * argv[])
{
    std::string filename =
        "obj/spot.obj";
    std::string out = gsFileManager::getCanonicRepresentation(
        gsFileManager::getPath(std::string(__FILE__), true) + "../output_poisson_immersed3d");
    index_t numRefine  = 3;
    index_t numElevate = 0;
    real_t  fill       = 0.9;  // fraction of [0,1]^3 spanned by the largest extent
    real_t  gamma      = 1e3;  // Nitsche penalty parameter (scaled by 1/h)
    index_t quadDepth  = 0;    // octree subdivision depth for cut-cell quadrature
    std::string indicator = "integralChange";
    real_t indicatorTol = 1e-2;
    bool    plot       = false;
    bool    noNitsche  = false;
    bool    momentFit  = false;
    real_t  alpha      = (real_t)0.0;
    // Inherited default of gsAlgoimGenericRule: n = quA*deg + quB points per
    // direction. See the order caveat in the file header.
    real_t  quA        = (real_t)1.0;
    bool    quadStats  = false;
    bool    errRefRule = false;

    gsCmdLine cmd("Poisson equation on an immersed 3D .obj domain via the Finite Cell Method (adaptive octree cut-cell quadrature).");
    cmd.addString("f", "file",            "Wavefront .obj mesh file",            filename);
    cmd.addString("o", "output",          "Output folder for ParaView files",    out);
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform refinement steps",  numRefine);
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps",    numElevate);
    cmd.addReal  ("",  "fill",            "Fill fraction of [0,1]^3 (0..1)",      fill);
    cmd.addReal  ("g", "gamma",           "Nitsche penalty parameter",           gamma);
    cmd.addInt   ("",  "quadDepth",       "Adaptive octree subdivision depth for cut-cell quadrature", quadDepth);
    cmd.addString("",  "indicator",       "Adaptive indicator: uniform, fallback or "
                                          "integralChange",                      indicator);
    cmd.addReal  ("",  "indicatorTol",    "integralChange acceptance tolerance", indicatorTol);
    cmd.addSwitch("plot",      "Create ParaView output",                        plot);
    cmd.addSwitch("noNitsche", "Disable the weak/Nitsche BC (interior forcing only)", noNitsche);
    cmd.addSwitch("momentFit", "Use moment-fitting compression for the cut-cell "
                               "VOLUME rule (the Nitsche surface rule is never "
                               "compressed)",                                   momentFit);
    cmd.addReal  ("",  "alpha",         "Fictitious-domain weight of the outside "
                                        "material in the moment-fitting rule "
                                        "(0 = drop it; only 0 is supported)",   alpha);
    cmd.addReal  ("",  "quA",           "Moment-fitting output order: n = quA*deg + quB "
                                        "points per direction. quA=1 reproduces only "
                                        "degree deg; the stiffness integrand needs "
                                        "quA=2. Applies to --momentFit only.",  quA);
    cmd.addSwitch("quadStats", "Print the interior/cut/surface quadrature-point split "
                               "and the CG solver diagnostics",                 quadStats);
    cmd.addSwitch("errRefRule", "Compute the error norms with the uncompressed adaptive "
                                "rule even when --momentFit assembles with the "
                                "compressed one (diagnostic)",                  errRefRule);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    // n = quA*deg + quB is clamped to one node, so a non-positive --quA would
    // silently degrade the moment-fitting output grid to a single point.
    GISMO_ENSURE(quA > 0, "--quA must be positive (got " << quA << ").");
    // Capability regression of the move to the core gsMomentRule: the FCM
    // blend w <- (1-alpha)w + alpha*w^G needs the analytic full-cell Gauss
    // weight w^G, which the geometry-free compressor does not have.
    GISMO_ENSURE(alpha == (real_t)0,
                 "--alpha (fictitious-domain blend) is not available with the "
                 "core gsMomentRule; use --alpha 0.");
    out = gsFileManager::getCanonicRepresentation(out);

    // -------------------------------------------------------------------------
    // 1. Load the triangle mesh
    // -------------------------------------------------------------------------
    gsSurfMesh mesh;
    if (!gsReadSurfMesh(filename, mesh))
    {
        gsWarn << "Failed to read a triangle mesh from: " << filename << "\n";
        return EXIT_FAILURE;
    }
    const std::string outputStem = gsFileManager::getBasename(filename);
    const std::size_t nVert = mesh.n_vertices();
    const std::size_t nTri  = mesh.n_faces();
    gsInfo << "Loaded mesh '" << filename << "': "
           << nVert << " vertices, " << nTri << " triangles.\n";

    // -------------------------------------------------------------------------
    // 2. Rescale / center the mesh into the unit parametric box [0,1]^3
    //    mapped = (p - center) * scale + 0.5,  scale = fill / maxExtent
    // -------------------------------------------------------------------------
    const real_t scale = gsNormalizeToUnitBox(mesh, fill);
    GISMO_UNUSED(scale);

    // -------------------------------------------------------------------------
    // 3. Analytic level set on the rescaled mesh
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3, 2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(mesh, bbox);

    // -------------------------------------------------------------------------
    // 4. Background box [0,1]^3 (identity map: parameter space == physical space)
    // -------------------------------------------------------------------------
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    const index_t deg = dbasis.maxCwiseDegree();

    // Manufactured solution: f = -Delta(u_exact), imposed weakly on {phi==0}.
    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
    gsFunctionExpr<> f_rhs("3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);

    gsFileManager::mkdir(out);

    gsInfo << "\nDegree of basis: " << deg << "\n";
    gsInfo << "Adaptive cut-cell quadrature octree depth: " << quadDepth << "\n";
    gsInfo << "Solving the Poisson equation on the immersed domain in [0,1]^3.\n";
    if (noNitsche)
        gsInfo << "Weak/Nitsche BC disabled: only the interior forcing is imposed.\n";
    gsInfo << "\nDoFs (dot1=assembled, dot2=solved, dot3=error): ";

    // Number of uniform sample points per direction used as fallback when
    // Algoim returns 0 nodes (re-entrant / multi-crossing cut cells).
    const index_t nFallback = 4;
    const int D = 3;

    gsVector<real_t> l2err(numRefine + 1), h1err(numRefine + 1);

    // ParaView collections spanning all refinement levels: each category is
    // written into its own subfolder and gathered by a single .pvd, where the
    // refinement level r plays the role of a time step.
    std::unique_ptr<gsParaviewCollection> colInterior, colCut, colAll,
                                          colBg, colSol, colExact;
    if (plot)
    {
        colInterior.reset(new gsParaviewCollection(out + "/points_interior/interior"));
        colCut     .reset(new gsParaviewCollection(out + "/points_cutcells/cutcells"));
        colAll     .reset(new gsParaviewCollection(out + "/points_all/all"));
        colBg      .reset(new gsParaviewCollection(out + "/background/background"));
        colSol     .reset(new gsParaviewCollection(out + "/solution/solution"));
        colExact   .reset(new gsParaviewCollection(out + "/exact/exact"));
    }

    // -------------------------------------------------------------------------
    // 5. Refinement loop: classify cells, assemble, solve, compute errors
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<3, real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3, real_t> *>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        real_t hmax = 0;
        for (std::size_t p = 0; p != dbasis.nBases(); ++p)
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

        // --- profiling -----------------------------------------------------
        gsStopwatch swTotal, swPhase;
        double t_classify = 0, t_ruleSetup = 0, t_reg = 0, t_interior = 0,
               t_cut = 0, t_nitsche = 0, t_setFrom = 0, t_solve = 0, t_error = 0;
        // ---------------------------------------------------------------------

        swPhase.restart();
        // Cut-cell classification of the background grid.
        gsImplicitTrimmedDomain<3, real_t> tr_domain(impl_fun, *tbsPtr);
        t_classify = swPhase.stop();

        swPhase.restart();
        // Quadrature rules -- exactly the same objects as the volume example.
        gsGaussRule<real_t>         gauss(gsVector<index_t, 3>::Constant(deg + 1));
        gsOptionList adaptiveOptions = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        adaptiveOptions.setInt("maxDepth", quadDepth);
        adaptiveOptions.setInt("nFallback", nFallback);
        adaptiveOptions.setString("indicator", indicator);
        adaptiveOptions.setReal("indicatorTol", indicatorTol);
        gsAlgoimAdaptiveRule<real_t> volumeRule(impl_fun, *tbsPtr, adaptiveOptions);
        gsOptionList surfaceOptions = adaptiveOptions;
        surfaceOptions.setInt("dim", D);
        // NOTE: the surface rule is deliberately built from adaptiveOptions and
        // is NEVER replaced by moment fitting (its points lie on {phi==0}).
        // Only the options are kept here: the Nitsche loop is OpenMP-parallel
        // and builds one rule per thread from them (the rules are not
        // thread-safe -- see the cut-cell loop).

        // Optional replacement of the cut-cell VOLUME rule only: the same
        // adaptive rule runs underneath (identical maxDepth / nFallback /
        // indicator / indicatorTol), but its point cloud is compressed onto the
        // tensor Gauss grid of the cell.  Rebuilt every refinement, so the
        // statistics printed below refer to this level only.
        // The rule classifies octree children by the midpoint + Lipschitz test
        // on phi -- exactly the test boxClass() below performs, so unlike in the
        // SAT-based volume driver there is no classification mismatch here.
        // The compressor OWNS its adaptive rule (the gsQuadRule hierarchy has
        // no clone(), so ownership transfer is the only non-slicing option).
        // Sets the output grid order AND (single key) the order of the
        // underlying Algoim rule -- see the order caveat in the header.
        gsOptionList momentFitOptions = adaptiveOptions;
        momentFitOptions.setReal("quA", quA);

        // Output points per direction, exactly as the former moment-fitting
        // rule computed them: n = round(quA*maxDegree) + quB, clamped to 1.
        const long nRaw = std::lround(momentFitOptions.askReal("quA", 1.0)
                                      * static_cast<double>(tbsPtr->maxDegree()))
                        + static_cast<long>(momentFitOptions.askInt("quB", 1));
        const index_t mfOrder1d = nRaw > 0 ? static_cast<index_t>(nRaw) : 1;

        // Factory rather than a single instance: the cut-cell loop below is
        // OpenMP-parallel and gsMomentRule is NOT thread-safe (mutable 1D
        // caches, Lagrange scratch buffers and statistics -- see its members),
        // so every thread builds and owns one. Same for gsAlgoimAdaptiveRule
        // (mutable Stats). Construction is cheap: the rules only store the
        // level set by reference plus their option lists.
        auto makeMomentRule = [&]() -> gsMomentRule<real_t>::uPtr
        {
            return gsMomentRule<real_t>::make(
                gsQuadRule<real_t>::uPtr(
                    new gsAlgoimAdaptiveRule<real_t>(impl_fun, *tbsPtr, momentFitOptions)),
                gsVector<index_t>::Constant(D, mfOrder1d));
        };

        gsMomentRule<real_t>::uPtr momentFitRule;
        if (momentFit)
            momentFitRule = makeMomentRule();
        t_ruleSetup = swPhase.stop();

        // Signed-distance box classification (phi is ~1-Lipschitz):
        //   phi(center) < -radius  => box fully inside  {phi<0}
        //   phi(center) >  radius  => box fully outside
        //   otherwise              => box still cut / uncertain
        // radius = half the box diagonal. Conservative: uncertain boxes are
        // treated as cut and refined further.
        auto boxClass = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c) -> int
        {
            gsMatrix<real_t> c(3,1);
            c.col(0) = (real_t)0.5 * (lo_c + hi_c);
            gsMatrix<real_t> v;
            impl_fun.eval_into(c, v);
            const real_t radius = (real_t)0.5 * (hi_c - lo_c).norm();
            if (v(0,0) >  radius) return  1;   // outside
            if (v(0,0) < -radius) return -1;   // inside
            return 0;                          // cut
        };

        // Entry point: adaptive octree quadrature for one background cut cell,
        // returning inside (Gauss) and cut (Algoim) points separately.
        // With \a useMomentFit the whole cell is handled by the moment-fitting
        // rule, which has no mapToSeparated: it emits a single compressed cloud
        // (returned in cutP/cutW), so the inside accumulator stays empty and
        // must be cleared explicitly -- mapTo() never touches it and stale
        // points from the previous cell would otherwise be re-assembled.
        auto cutCellQuadrature = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c,
                                     gsMatrix<real_t> & insP, gsVector<real_t> & insW,
                                     gsMatrix<real_t> & cutP, gsVector<real_t> & cutW,
                                     bool useMomentFit)
        {
            if (useMomentFit)
            {
                insP.resize(D, 0);
                insW.resize(0);
                momentFitRule->mapTo(lo_c, hi_c, cutP, cutW);
                return;
            }
            volumeRule.mapToSeparated(lo_c, hi_c, insP, insW, cutP, cutW, boxClass);
        };

        // ---------------------------------------------------------------
        // Explicit assembly. The geometry map is the identity on [0,1]^3,
        // so parametric gradients equal physical gradients (Jacobian = I).
        // ---------------------------------------------------------------
        const index_t nDofs = tbsPtr->size();
        gsVector<real_t> F_vec(nDofs);
        F_vec.setZero();
        gsSparseEntries<real_t> triplets;

        // Accumulate stiffness (+ load) from one batch of volume points.
        auto addVolume = [&](const gsMatrix<real_t> & pts, const gsVector<real_t> & wts,
                              real_t coeff = (real_t)1, bool withRhs = true)
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

        gsMatrix<real_t> pts_tmp; gsVector<real_t> wts_tmp;
        gsMatrix<real_t> insP, cutP; gsVector<real_t> insW, cutW;

        // Tiny full-background regularization: keeps K nonsingular for DOFs
        // whose support never touches {phi<0}.
        {
            swPhase.restart();
            const real_t alpha = 1e-10;
            for (auto it = tbsPtr->domain()->beginAll(); it != tbsPtr->domain()->endAll(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
                addVolume(pts_tmp, wts_tmp, alpha, /*withRhs=*/false);
            }
            t_reg = swPhase.stop();
        }

        // Total number of volume quadrature points actually used (cost metric).
        // Split into the fully-interior background cells (plain Gauss, never
        // touched by moment fitting) and the cut cells (the only part the
        // moment-fitting rule can compress), so the saving is not diluted by a
        // term that cannot change.
        index_t nQuad = 0, nQuadInterior = 0, nQuadCut = 0, nQuadSurf = 0;

        // Statistics of THIS assembly pass only. The assembly pass now tallies
        // into per-thread rules (reduced into mfs below), so the shared
        // momentFitRule only ever sees the error loop's and the plot branch's
        // calls -- but reset it anyway so that a later read of it cannot pick
        // up a previous refinement level's counts.
        if (momentFit) momentFitRule->resetStats();

        // int_Omega grad(u).grad(v)  and  int_Omega f.v,  Omega = {phi<0}.
        swPhase.restart();
        index_t nInteriorCells = 0;
        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addVolume(pts_tmp, wts_tmp);
            nQuad += pts_tmp.cols();
            nQuadInterior += pts_tmp.cols();
            ++nInteriorCells;
        }
        t_interior = swPhase.stop();

        // Cache of this pass's cut-cell quadrature (insP/insW/cutP/cutW per
        // cell), reused by the error loop below instead of recomputing it: a
        // cut cell's quadrature is the expensive part (Algoim/moment-fitting
        // + gsMeshSignedDist queries), and the error loop otherwise reran the
        // exact same cutCellQuadrature() calls from scratch (measured to cost
        // as much as the assembly pass itself). Only valid when the error
        // loop would ask for the SAME rule as assembly, i.e. not when
        // errRefRule forces the uncompressed rule while momentFit assembled
        // with the compressed one.
        struct CutCellQuad { gsMatrix<real_t> insP, cutP; gsVector<real_t> insW, cutW; };
        std::vector<CutCellQuad> cutCache;
        const bool cutCacheable = !(momentFit && errRefRule);

        swPhase.restart();

        // Materialize the cut-cell boxes: the domain iterator is not
        // random-access, so the OpenMP loop below cannot run on it directly.
        std::vector<std::pair<gsVector<real_t>, gsVector<real_t> > > cutBoxes;
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
            cutBoxes.push_back(std::make_pair(it.lowerCorner(), it.upperCorner()));
        const index_t nCutCells = static_cast<index_t>(cutBoxes.size());

        // Two-phase, because this loop is the bulk of the runtime and its two
        // halves have opposite parallel character:
        //   (1) building a cell's quadrature is a pure function of the cell box
        //       (level set + Algoim / moment fitting), embarrassingly parallel,
        //       and is where essentially all of the time goes -- the emitted
        //       clouds are tiny compared to the octree work behind them;
        //   (2) scattering it into triplets/F_vec touches shared state, is
        //       cheap, and must stay in cell order anyway so that the assembled
        //       matrix -- and hence the CG iterate sequence -- stays
        //       bit-identical to the serial run.
        // Phase 1 therefore fills per-cell slots in parallel, phase 2 replays
        // them serially in cell order.
        std::vector<CutCellQuad> cutQuad(cutBoxes.size());
        gsMomentRule<real_t>::Stats mfs;
#       pragma omp parallel
        {
            // Per-thread rules: both gsMomentRule and gsAlgoimAdaptiveRule keep
            // mutable per-call scratch and statistics, so sharing one instance
            // across threads is a data race. See makeMomentRule() above.
            gsAlgoimAdaptiveRule<real_t> thrVolumeRule(impl_fun, *tbsPtr, adaptiveOptions);
            gsMomentRule<real_t>::uPtr   thrMomentRule;
            if (momentFit) thrMomentRule = makeMomentRule();

#           pragma omp for schedule(dynamic)
            for (index_t c = 0; c < nCutCells; ++c)
            {
                CutCellQuad & q = cutQuad[c];
                if (momentFit)
                {
                    // Moment fitting emits one compressed cloud for the whole
                    // cell (no mapToSeparated), so the inside part stays empty.
                    q.insP.resize(D, 0);
                    q.insW.resize(0);
                    thrMomentRule->mapTo(cutBoxes[c].first, cutBoxes[c].second,
                                         q.cutP, q.cutW);
                }
                else
                    thrVolumeRule.mapToSeparated(cutBoxes[c].first, cutBoxes[c].second,
                                                 q.insP, q.insW, q.cutP, q.cutW, boxClass);
            }

            // Stats are pure sums (and a min for minWeight), so the reduction
            // reproduces the serial tally exactly, thread count irrespective.
            if (momentFit)
#               pragma omp critical
                mfs += thrMomentRule->stats();
        }

        for (index_t c = 0; c < nCutCells; ++c)
        {
            addVolume(cutQuad[c].insP, cutQuad[c].insW);
            addVolume(cutQuad[c].cutP, cutQuad[c].cutW);
            nQuad    += cutQuad[c].insP.cols() + cutQuad[c].cutP.cols();
            nQuadCut += cutQuad[c].insP.cols() + cutQuad[c].cutP.cols();
        }
        t_cut = swPhase.stop();
        if (cutCacheable) cutCache = give(cutQuad);

        // ---------------------------------------------------------------
        // 6. Weak (Nitsche) immersed boundary condition on {phi==0}.
        //    Normal: n = grad(phi)/|grad(phi)| (finite differences via
        //    gsMeshSignedDist). Disable with --noNitsche if this ever proves
        //    fragile -- without it, the interior forcing is still solved
        //    correctly, but u=0 is not imposed on the immersed boundary, so
        //    the solution will not converge to u_exact.
        // ---------------------------------------------------------------
        if (!noNitsche)
        {
            swPhase.restart();

            // Same two-phase split as the cut-cell volume loop above: the
            // surface rule's octree and the finite-difference normals are the
            // expensive, pure-per-cell part and run in parallel; the assembly
            // replay stays serial and in cell order.  phi_d is computed here
            // rather than in phase 2 because gsFunction's central-difference
            // deriv_into() costs 4 level-set evaluations per direction per
            // point -- 12 BVH queries for every surface point -- which belongs
            // on the parallel side.
            struct SurfCellQuad
            {
                gsMatrix<real_t> spts, phi_d;
                gsVector<real_t> swts;
            };
            std::vector<SurfCellQuad> surfQuad(cutBoxes.size());
            gsAlgoimAdaptiveRule<real_t>::Stats surfStats;
#           pragma omp parallel
            {
                gsAlgoimAdaptiveRule<real_t> thrSurfaceRule(impl_fun, *tbsPtr, surfaceOptions);

#               pragma omp for schedule(dynamic)
                for (index_t c = 0; c < nCutCells; ++c)
                {
                    gsMatrix<real_t> sInterior(3,0);
                    gsVector<real_t> sInteriorW;
                    SurfCellQuad & s = surfQuad[c];
                    thrSurfaceRule.mapToSeparated(cutBoxes[c].first, cutBoxes[c].second,
                                                  sInterior, sInteriorW, s.spts, s.swts,
                                                  boxClass);
                    if (s.spts.cols() > 0)
                        impl_fun.deriv_into(s.spts, s.phi_d);
                }

#               pragma omp critical
                surfStats += thrSurfaceRule.stats();
            }

            for (index_t c = 0; c < nCutCells; ++c)
            {
                const gsMatrix<real_t> & spts  = surfQuad[c].spts;
                const gsVector<real_t> & swts  = surfQuad[c].swts;
                const gsMatrix<real_t> & phi_d = surfQuad[c].phi_d;
                if (spts.cols() == 0) continue;

                gsMatrix<index_t> act;
                gsMatrix<real_t>  bv, bd, gDv;
                tbsPtr->active_into(spts, act);
                tbsPtr->eval_into  (spts, bv);
                tbsPtr->deriv_into (spts, bd);
                u_exact.eval_into  (spts, gDv);      // 1 x numSurfPts
                const index_t na = act.rows();
                nQuadSurf += spts.cols();

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

                        // RHS: -(grad(phi_i).n)g_D + (gamma/h)phi_i g_D
                        F_vec(di) += w * (-gi_n * gDv(0,q) + (gamma/hmax) * vi * gDv(0,q));

                        for (index_t j = 0; j < na; ++j)
                        {
                            const index_t dj = act(j,q);
                            const real_t  vj = bv(j,q);
                            real_t gj_n = 0;
                            for (int d = 0; d < D; ++d)
                                gj_n += bd(j*D+d,q) * phi_d(d,q) / nNorm;

                            // -(gi.n)vj - vi(gj.n) + (gamma/h)vi.vj
                            const real_t kij = -gi_n*vj - vi*gj_n + (gamma/hmax)*vi*vj;
                            triplets.add(di, dj, w * kij);
                        }
                    }
                }
            }
            if (surfStats.nZeroSurfaceLeaves > 0)
                gsInfo << " zeroSurfaceLeaves=" << surfStats.nZeroSurfaceLeaves
                       << " unresolvedSurfaceMeasure~"
                       << surfStats.unresolvedSurfaceMeasure;
            t_nitsche = swPhase.stop();
        }

        // Assemble sparse matrix and solve.
        swPhase.restart();
        gsSparseMatrix<real_t> K(nDofs, nDofs);
        K.setFrom(triplets);
        triplets.clear();
        t_setFrom = swPhase.stop();

        gsInfo << nDofs << "." << std::flush;

        gsVector<real_t> solVector(nDofs);
        // Solver health. CGDiagonal does NOT signal failure on an indefinite
        // matrix (which moment fitting can produce via negative weights): it
        // just returns a badly converged vector. Record iterations/residual so
        // a "converged" error norm can never be reported on a junk solve.
        index_t cgIters = -1;
        real_t  cgError = -1;
        bool    solveOk = true;
        {
            swPhase.restart();
#           ifdef GISMO_WITH_PARDISO
                gsSparseSolver<real_t>::PardisoLDLT slvr;
                slvr.compute(K);
                solVector = slvr.solve(F_vec);
#           else
                gsSparseSolver<real_t>::CGDiagonal slvr;
                slvr.compute(K);
                solVector = slvr.solve(F_vec);
                cgIters = static_cast<index_t>(slvr.iterations());
                cgError = slvr.error();
                solveOk = slvr.succeed();
#           endif
            t_solve = swPhase.stop();
        }
        gsInfo << "." << std::flush;

        // Error norms, integrated over {phi<0} with the same quadrature.
        real_t l2sq = 0, h1sq = 0;
        auto addError = [&](const gsMatrix<real_t> & pts, const gsVector<real_t> & wts)
        {
            if (pts.cols() == 0) return;
            gsMatrix<index_t> act;
            gsMatrix<real_t>  bv, bd, uExV, uExD;
            tbsPtr->active_into(pts, act);
            tbsPtr->eval_into  (pts, bv);
            tbsPtr->deriv_into (pts, bd);
            u_exact.eval_into  (pts, uExV);        // 1 x numPts
            u_exact.deriv_into (pts, uExD);        // 3 x numPts
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

        swPhase.restart();
        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addError(pts_tmp, wts_tmp);
        }
        // The error norms are integrated with the SAME rule the assembly used,
        // so with --momentFit the L2/H1 values are themselves computed on the
        // compressed cloud: at low quA the error integrand (u_h-u_exact)^2 is
        // far outside what the output grid interpolates exactly, and the
        // reported error can move (in either direction) for reasons unrelated
        // to solution quality. --errRefRule integrates the norms with the
        // uncompressed adaptive rule instead, which separates the two effects.
        index_t cutIdx = 0;
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it, ++cutIdx)
        {
            if (cutCacheable)
            {
                addError(cutCache[cutIdx].insP, cutCache[cutIdx].insW);
                addError(cutCache[cutIdx].cutP, cutCache[cutIdx].cutW);
            }
            else
            {
                cutCellQuadrature(it.lowerCorner(), it.upperCorner(), insP, insW, cutP, cutW,
                                  momentFit && !errRefRule);
                addError(insP, insW);
                addError(cutP, cutW);
            }
        }
        t_error = swPhase.stop();
        l2err[r] = math::sqrt(l2sq);
        h1err[r] = l2err[r] + math::sqrt(h1sq);
        gsInfo << ". " << std::flush;

        gsInfo << "\nprofile[r=" << r << "]: classify=" << t_classify
               << "s ruleSetup=" << t_ruleSetup
               << "s reg=" << t_reg
               << "s interior(" << nInteriorCells << " cells)=" << t_interior
               << "s cutVol(" << nCutCells << " cells)=" << t_cut
               << "s nitsche=" << t_nitsche
               << "s setFrom=" << t_setFrom
               << "s solve=" << t_solve
               << "s error=" << t_error
               << "s total=" << swTotal.stop() << "s\n";

        // ---------------------------------------------------------------
        // 7. ParaView output, one subfolder + .pvd per category:
        //     - interior / cut-cell / all quadrature points (4 x N clouds).
        //       Inside sub-cells of a cut cell (cheap Gauss) are tagged as
        //       interior; Algoim leaves are tagged as cut.
        //     - the background mesh,
        //     - the numerical and exact solution as point clouds sampled at
        //       the quadrature points inside {phi<0} (4 x N: xyz + value).
        //   The refinement level r is used as the .pvd time step.
        // ---------------------------------------------------------------
        if (plot)
        {
            const std::string rs = std::to_string(r);
            gsMatrix<real_t> pts, phys, quadCut(4, 0), quadIn(4, 0), quadAll(4, 0),
                             solNum(4, 0), solEx(4, 0);
            gsVector<real_t> wts;

            // Append physical points + a scalar (weight or value) as 4 x N.
            auto append = [&](gsMatrix<real_t> & acc, const gsMatrix<real_t> & phys,
                              const gsVector<real_t> & scalar)
            {
                const index_t c = acc.cols();
                acc.conservativeResize(4, c + phys.cols());
                acc.block(0, c, 3, phys.cols()) = phys;
                acc.row(3).segment(c, phys.cols()) = scalar.transpose();
            };

            // Evaluate u_h and u_exact at parametric pts and append them as
            // point clouds (value in row 3).
            auto appendSolution = [&](const gsMatrix<real_t> & pts, const gsMatrix<real_t> & phys)
            {
                gsMatrix<index_t> act;
                gsMatrix<real_t>  bv, uEx;
                tbsPtr->active_into(pts, act);
                tbsPtr->eval_into  (pts, bv);
                u_exact.eval_into  (pts, uEx);
                gsVector<real_t> uh(pts.cols());
                for (index_t q = 0; q < pts.cols(); ++q)
                {
                    real_t v = 0;
                    for (index_t i = 0; i < act.rows(); ++i)
                        v += bv(i,q) * solVector(act(i,q));
                    uh[q] = v;
                }
                append(solNum, phys, uh);
                append(solEx,  phys, uEx.row(0).transpose());
            };

            // Cut cells: inside sub-cells -> interior, Algoim leaves -> cut.
            // This asks for the SAME rule the assembly loop used (momentFit,
            // not momentFit && !errRefRule as the error loop does), so the
            // assembly cache applies verbatim under the same cutCacheable
            // condition -- without it this was a third full recomputation of
            // every cut cell's quadrature, as expensive as assembly itself.
            index_t plotIdx = 0;
            for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it, ++plotIdx)
            {
                if (cutCacheable)
                {
                    insP = cutCache[plotIdx].insP; insW = cutCache[plotIdx].insW;
                    cutP = cutCache[plotIdx].cutP; cutW = cutCache[plotIdx].cutW;
                }
                else
                    cutCellQuadrature(it.lowerCorner(), it.upperCorner(), insP, insW, cutP, cutW,
                                      momentFit);
                if (insP.cols() > 0)
                {
                    mp.patch(0).eval_into(insP, phys);
                    append(quadIn,  phys, insW);
                    append(quadAll, phys, insW);
                    appendSolution(insP, phys);
                }
                if (cutP.cols() > 0)
                {
                    mp.patch(0).eval_into(cutP, phys);
                    append(quadCut, phys, cutW);
                    append(quadAll, phys, cutW);
                    appendSolution(cutP, phys);
                }
            }
            // Interior (fully inside {phi<0}) background Gauss quadrature points.
            for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                append(quadIn,  phys, wts);
                append(quadAll, phys, wts);
                appendSolution(pts, phys);
            }
            if (quadCut.cols() > 0)
            {
                gsWriteParaviewPoints(quadCut, out + "/points_cutcells/cutcells_r" + rs);
                colCut->addPart("cutcells_r" + rs + ".vtp", r);
            }
            if (quadIn.cols() > 0)
            {
                gsWriteParaviewPoints(quadIn, out + "/points_interior/interior_r" + rs);
                colInterior->addPart("interior_r" + rs + ".vtp", r);
            }
            if (quadAll.cols() > 0)
            {
                gsWriteParaviewPoints(quadAll, out + "/points_all/all_r" + rs);
                colAll->addPart("all_r" + rs + ".vtp", r);
            }

            // Background mesh at every refinement level.
            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/background/background_r" + rs);
            colBg->addPart("background_r" + rs + ".vtp", r);

            // Numerical and exact solution as point clouds.
            if (solNum.cols() > 0)
            {
                gsWriteParaviewPoints(solNum, out + "/solution/solution_r" + rs);
                colSol->addPart("solution_r" + rs + ".vtp", r);
                gsWriteParaviewPoints(solEx, out + "/exact/exact_r" + rs);
                colExact->addPart("exact_r" + rs + ".vtp", r);
            }
        }

        gsInfo << " (h=" << hmax << ", gamma/h=" << gamma/hmax
               << ", quadPts=" << nQuad << ")\n";

        // Quadrature-point split of the assembly pass. intQPs (fully interior
        // background cells) and surfQPs (Nitsche points on {phi==0}) are
        // untouched by moment fitting by construction; only cutQPs can shrink.
        if (quadStats || momentFit)
            gsInfo << "quad: intQPs=" << nQuadInterior
                   << " cutQPs=" << nQuadCut
                   << " volQPs=" << nQuad
                   << " surfQPs=" << nQuadSurf
                   << " | cg: ok=" << (solveOk ? 1 : 0)
                   << " iters=" << cgIters
                   << " res=" << std::scientific << std::setprecision(3) << cgError << "\n";

        // -- Moment-fitting compression statistics (this refinement and this
        //    assembly pass only).  outQPs/minW/negW describe the compressed
        //    cells only; cells whose underlying cloud already fits in the
        //    output grid are passed through unchanged (passThrough) and
        //    contribute passQPs, not outQPs.  The honest emitted point count is
        //    emitted = outQPs + passQPs (valid at alpha == 0).
        if (momentFit)
            gsInfo << "momfit: outQPs=" << mfs.nOutputQPs
                   << " passQPs=" << mfs.nPassThroughQPs
                   << " emitted=" << (mfs.nOutputQPs + mfs.nPassThroughQPs)
                   << " underQPs=" << mfs.nUnderlyingQPs
                   << " passThrough=" << mfs.nPassThroughElements
                   << " minW=" << std::scientific << std::setprecision(3) << mfs.minWeight
                   << " negW=" << mfs.nNegativeWeights << "\n";

        // Restore the stream format so the next level's "(h=..., gamma/h=...)"
        // line is printed exactly as in the non-momentFit run.
        if (quadStats || momentFit)
            gsInfo << std::defaultfloat << std::setprecision(6);
    }


    gsInfo << "\nL2 error: " << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err.transpose() << "\n";

    if (numRefine > 0)
    {
        gsInfo << "\nEoC (L2): " << std::fixed << std::setprecision(2)
               << (l2err.head(numRefine).array() / l2err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
        gsInfo << "EoC (H1): " << std::fixed << std::setprecision(2)
               << (h1err.head(numRefine).array() / h1err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
    }

    if (plot)
    {
        // Immersed geometry (the input triangle mesh), written once.
            gsWriteParaview(mesh, out + "/" + outputStem + "_geometry");

        colInterior->save(); colCut->save(); colAll->save();
        colBg->save(); colSol->save(); colExact->save();

        gsInfo << "\nParaView files written to: " << out << "\n";
    }

    return EXIT_SUCCESS;
}
