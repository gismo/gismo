/** @file immersed_face_iterators_example.cpp

    @brief Reports and plots the skeleton/ghost face sets of an immersed
    background grid (a circle/ellipse/sphere level set immersed in a uniform
    tensor-product grid on [0,1]^d), without assembling or solving any PDE.

    Per uniform refinement level this example
      - builds an implicit trimmed domain (gsImplicitTrimmedDomain) on a
        uniform tensor B-spline background grid,
      - enumerates its skeleton faces (between two active elements) and its
        ghost faces (skeleton faces touching at least one cut element)
        through the gsDomain face-iterator API (beginSkeleton()/beginGhost()),
      - checks BOTH face sets against an independent brute-force reference
        rebuilt from element-box neighbour pairs only -- that block never
        calls a face-iterator or a *Faces() count (see the limitation below),
        and
      - reports Sum [[d^k u]]^2 over the skeleton and over the ghost set, for
        a deterministic "random" spline, via gsExprEvaluator::integralSkeleton/
        integralGhost and gsFeSolution::jump().
    With --plot it writes ParaView box sets for the skeleton faces, the
    ghost faces and the interior/cut element classification at the finest
    refinement level -- the deliverable this example exists for; the
    terminal table is secondary.

    What the brute-force check does and does not discriminate: both the
    library face iterators and the brute-force reference read their signs
    from the SAME kd-tree leaf classification (beginInterior()/beginBdr()),
    so agreement rules out a bug in the level-0 sign-grid mapping or in the
    SkeletonFace/GhostFace predicates, but it cannot catch a leaf that the
    classifier itself mis-signed.

    Example command lines:
      ./immersed_face_iterators_example
      ./immersed_face_iterators_example --dim 3 --shape sphere -r 1
      ./immersed_face_iterators_example --shape ellipse --dy 0.6
      ./immersed_face_iterators_example --plot -o /tmp/faces2d
      ./immersed_face_iterators_example --dim 3 --shape sphere -r 1 --plot -o /tmp/faces3d

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

using namespace gismo;

namespace {

// Integer grid coordinate of a point on the uniform n-element grid of [0,1]^d.
// Duplicated (not shared) from optional/gsAlgoim/unittests/gsFaceIterator_test.cpp:
// the brute-force reference built in run<d>() below must stay independent of
// the library's face-iterator API, so this helper intentionally has its own
// copy here rather than living in a header both files include.
index_t gridCoord(const real_t x, const index_t n)
{ return static_cast<index_t>(math::round(x * (real_t)n)); }

// Flat level-0 background-element index of a point's grid coordinates on the
// uniform n-per-direction grid of [0,1]^d, direction 0 running fastest:
// flat = sum_j i_j*n^j. Assumes a uniform grid over [0,1]^d, true here by
// construction (the background basis is refined only through uniformRefine()).
size_t flatIndex(const gsVector<real_t> & p, const index_t n)
{
    size_t flat = 0, stride = 1;
    for (index_t j = 0; j != p.rows(); ++j)
    { flat += (size_t)gridCoord(p[j], n) * stride; stride *= (size_t)n; }
    return flat;
}

// Inverse of flatIndex: multi-index of a flat level-0 element index on the
// uniform n-per-direction grid of dimension d (d<=3 here).
void decodeFlat(size_t flat, const index_t n, const index_t d, index_t idx[3])
{
    for (index_t j = 0; j != d; ++j) { idx[j] = (index_t)(flat % (size_t)n); flat /= (size_t)n; }
}

} // anonymous namespace

/// The embedding box [0,1]^d as a single-patch multi-patch, used as the
/// geometry map for both the trimmed domain's background basis and the
/// PDE-free assembler that provides getSolution()/jump().
template<short_t d> gsMultiPatch<real_t> unitBox();
template<> gsMultiPatch<real_t> unitBox<2>()
{ gsMultiPatch<real_t> mp; mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1,0,0)); return mp; }
template<> gsMultiPatch<real_t> unitBox<3>()
{ gsMultiPatch<real_t> mp; mp.addPatch(gsNurbsCreator<real_t>::BSplineCube(1,0,0,0)); return mp; }

struct Config
{
    index_t     dim        = 2;
    index_t     numRefine  = 3;
    std::string shape      = "circle";
    real_t      dx         = 1.0;
    real_t      dy         = 1.0;
    real_t      dz         = 1.0;
    real_t      radius     = -1.0;
    index_t     degree     = 2;
    bool        plot       = false;
    std::string out        = "output_face_iterators";
};

template<short_t d> int run(const Config & cfg);

int main(int argc, char *argv[])
{
    Config cfg;

    gsCmdLine cmd("Reports and plots the skeleton/ghost face sets of an "
                  "immersed background grid. PDE-free: no assembly, no "
                  "linear solve -- see the module unit tests for matrix-"
                  "level properties.");
    cmd.addInt   ("",     "dim",           "Spatial dimension (2 or 3)", cfg.dim);
    cmd.addInt   ("r",    "uniformRefine", "Number of uniform refinement steps after the base grid", cfg.numRefine);
    cmd.addString("",     "shape",         "Level-set shape: circle | ellipse | sphere", cfg.shape);
    cmd.addReal  ("",     "dx",            "Semi-axis scale factor in x (shape=ellipse only)", cfg.dx);
    cmd.addReal  ("",     "dy",            "Semi-axis scale factor in y (shape=ellipse only)", cfg.dy);
    cmd.addReal  ("",     "dz",            "Semi-axis scale factor in z (shape=ellipse only)", cfg.dz);
    cmd.addReal  ("",     "radius",        "Level-set radius; <=0 selects the auto (grazing) radius", cfg.radius);
    cmd.addInt   ("",     "degree",        "Spline degree of the background basis", cfg.degree);
    cmd.addSwitch("plot", "Write ParaView box sets for the skeleton/ghost faces and elements", cfg.plot);
    cmd.addString("o",    "output",        "Output directory", cfg.out);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (2 != cfg.dim && 3 != cfg.dim)
    { gsWarn << "--dim must be 2 or 3\n"; return EXIT_FAILURE; }

    if ("circle" != cfg.shape && "ellipse" != cfg.shape && "sphere" != cfg.shape)
    { gsWarn << "--shape must be circle, ellipse or sphere\n"; return EXIT_FAILURE; }

    if (("circle" == cfg.shape || "ellipse" == cfg.shape) && 2 != cfg.dim)
    { gsWarn << "--shape " << cfg.shape << " is only valid for --dim 2\n"; return EXIT_FAILURE; }

    if ("sphere" == cfg.shape && 3 != cfg.dim)
    { gsWarn << "--shape sphere is only valid for --dim 3\n"; return EXIT_FAILURE; }

    if (cfg.numRefine < 0)
    { gsWarn << "--uniformRefine must be >= 0\n"; return EXIT_FAILURE; }

    // A negative degree faults inside gsBasis::setDegree below rather than
    // being rejected there, so it must not reach it.
    if (cfg.degree < 1)
    { gsWarn << "--degree must be >= 1\n"; return EXIT_FAILURE; }

    if (2 == cfg.dim) return run<2>(cfg);
    return run<3>(cfg);
}

template<short_t d>
int run(const Config & cfg)
{
    const index_t N0 = 4;   // background-grid elements per direction at level 0

    std::string outPath = cfg.out;
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath) outPath = gsFileManager::getCurrentPath() + "/" + outPath;
    const std::string out = gsFileManager::getCanonicRepresentation(outPath);
    if (cfg.plot) gsFileManager::mkdir(out);

    // ------------------------------------------------------------------
    // Level-set construction. gsImplicitTrimmedDomain holds phi NON-OWNED
    // (memory::make_shared_not_owned), so phi must outlive every trim built
    // below -- it is declared here, once, before the level loop. The radius
    // is likewise fixed once: recomputing it per level would mutate the
    // geometry mid-table and make every column incomparable across rows.
    //
    // The default radius is chosen so the level set grazes a mesh line by
    // delta = 1e-2 * h_finest: xStar is a multiple of the level-0 cell size
    // h0, hence a mesh line at every dyadic refinement level. delta is fixed
    // in ABSOLUTE terms, so the relative overlap delta/h_r = sliver*2^(r-numRefine)
    // GROWS with refinement: it is the intended 1e-2 of a cell at the finest
    // level and a smaller fraction -- i.e. a MORE degenerate cut -- at coarser
    // ones, which is why only the finest level's classification is asserted.
    // This is what makes ghost faces visibly motivated at the default settings.
    // ------------------------------------------------------------------
    const real_t c0 = 0.5 + 1e-3;                     // centre, off every dyadic mesh line
    const bool   autoRadius = (cfg.radius <= 0.0);
    real_t radius = cfg.radius;
    real_t delta = 0.0, h0 = 0.0, hFine = 0.0, xStar = 0.0;
    const real_t sliver = 1e-2;                       // delta / h_finest, auto-radius path only
    if (autoRadius)
    {
        h0    = 1.0 / (real_t)N0;
        hFine = 1.0 / (real_t)(N0 << cfg.numRefine);
        const real_t nominal = 0.3;
        const index_t m = (index_t)std::floor((c0 + nominal) / h0);
        xStar  = (real_t)m * h0;
        radius = xStar - c0 + sliver * hFine;
        GISMO_ENSURE(radius > 0, "Auto-radius computation produced a non-positive radius.");
        delta = sliver * hFine;
    }

    std::ostringstream phiStr;
    phiStr << std::setprecision(17);
    if ("circle" == cfg.shape)
        phiStr << "sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)-" << radius;
    else if ("sphere" == cfg.shape)
        phiStr << "sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2+(z-" << c0 << ")^2)-" << radius;
    else // ellipse
        phiStr << "sqrt(((x-" << c0 << ")/" << cfg.dx << ")^2+((y-" << c0 << ")/" << cfg.dy << ")^2)-" << radius;
    gsFunctionExpr<real_t> phi(phiStr.str(), d);      // NAMED, non-owned by every trim below

    // The finest-level sliver-cell classification check (below) is only a
    // valid probe of THIS geometric construction on the axis-symmetric auto
    // path (shape circle/sphere, unit semi-axes); an ellipse or a
    // user-supplied radius does not graze xStar by construction, so the
    // check is skipped there -- the level-by-level PASS/FAIL table still
    // covers those configurations.
    const bool doSliverCheck = autoRadius
        && ("circle" == cfg.shape || "sphere" == cfg.shape)
        && 1.0 == cfg.dx && 1.0 == cfg.dy && (2 == d || 1.0 == cfg.dz);

    gsInfo << "dim=" << d << "  shape=" << cfg.shape << "  degree=" << cfg.degree
           << "  centre=" << c0 << "  radius=" << (autoRadius ? "auto=" : "user=")
           << std::setprecision(10) << radius << std::setprecision(6);
    if (autoRadius)
        gsInfo << "  delta=" << std::scientific << std::setprecision(4) << delta
               << "  delta/h_finest=" << (delta / hFine) << std::fixed << std::setprecision(6);
    gsInfo << "  N0=" << N0 << "  levels=" << (cfg.numRefine + 1) << "\n\n";

    // Build the mesh hierarchy through gsMultiBasis (the assembler needs it
    // anyway) and derive the trimmed domain's tensor basis from it at every
    // level, so the trimmed-domain breaks and the assembler's space breaks
    // are identical by construction.
    gsMultiPatch<real_t> mp = unitBox<d>();
    gsMultiBasis<real_t> dbasis(mp);
    dbasis.setDegree(static_cast<short_t>(cfg.degree));   // elevate BEFORE refining, or continuity drops to C^0
    dbasis.uniformRefine(N0 - 1);                          // level 0: N0 elements per direction

    gsInfo << " lvl     N    #int    #cut   #skel  #ghost  ghost/skel   bf#skel  bf#ghost   check\n";

    bool allPass = true;

    for (index_t r = 0; r <= cfg.numRefine; ++r)
    {
        const index_t N = N0 << r;
        const real_t  h = 1.0 / (real_t)N;

        const gsTensorBSplineBasis<d,real_t> & tb =
            static_cast<const gsTensorBSplineBasis<d,real_t>&>(dbasis.basis(0));

        // Fresh trimmed domain every level: never hold a face iterator (or
        // the sign-grid pointer behind it) across a rebuild of this object.
        memory::shared_ptr< gsImplicitTrimmedDomain<d,real_t> > trim =
            memory::make_shared(new gsImplicitTrimmedDomain<d,real_t>(phi, tb));
        GISMO_ENSURE(1 == trim->numLevels(), "gsTrimmedDomain requires a single kd-tree level.");

        // ------------------------------------------------------------------
        // Brute-force reference, from element boxes only. This block never
        // calls numSkeletonFaces()/numGhostFaces()/beginSkeleton()/
        // beginGhost()/endSkeleton()/endGhost() or any face iterator: it
        // reads only the interior/cut element classification and derives the
        // face sets itself, so it is independent of the predicates it checks.
        // ------------------------------------------------------------------
        size_t total = 1;
        for (short_t j = 0; j != d; ++j) total *= (size_t)N;
        std::vector<short_t> sgn(total, (short_t)+1);     // default: exterior
        size_t nInterior = 0, nCut = 0;

        gsDomain<real_t>::iterator eInt = trim->template end<InteriorSign>();   // hoisted end
        for (gsDomain<real_t>::iterator it = trim->beginInterior(); it < eInt; ++it)
        { sgn[ flatIndex(it.lowerCorner(), N) ] = -1; ++nInterior; }

        gsDomain<real_t>::iterator eCut = trim->endBdr(boundary::none);        // hoisted end
        for (gsDomain<real_t>::iterator it = trim->beginBdr(boundary::none); it < eCut; ++it)
        { sgn[ flatIndex(it.lowerCorner(), N) ] = 0; ++nCut; }

        std::vector<size_t> stride(d);
        { size_t s = 1; for (short_t j = 0; j != d; ++j) { stride[j] = s; s *= (size_t)N; } }

        std::set< std::pair<short_t,size_t> > bSkel, bGhost;   // canonicalised from the LEFT element
        for (size_t f = 0; f != total; ++f)
        {
            index_t idx[3];
            decodeFlat(f, N, d, idx);
            for (short_t dir = 0; dir != d; ++dir)
            {
                if (idx[dir] + 1 > N-1) continue;               // no neighbour above
                const size_t rgt = f + stride[dir];
                if (sgn[f] <= 0 && sgn[rgt] <= 0)
                {
                    bSkel.insert(std::make_pair(dir, f));
                    if (0 == sgn[f] || 0 == sgn[rgt]) bGhost.insert(std::make_pair(dir, f));
                }
            }
        }

        // Library side (face iterators allowed here).
        std::set< std::pair<short_t,size_t> > lSkel, lGhost;
        gsDomain<real_t>::iterator endSkel = trim->endSkeleton();
        for (gsDomain<real_t>::iterator it = trim->beginSkeleton(); it < endSkel; ++it)
            lSkel.insert(std::make_pair(it.side().direction(), it.leftElementId()));
        gsDomain<real_t>::iterator endGhost = trim->endGhost();
        for (gsDomain<real_t>::iterator it = trim->beginGhost(); it < endGhost; ++it)
            lGhost.insert(std::make_pair(it.side().direction(), it.leftElementId()));

        const size_t nSkel  = trim->numSkeletonFaces();
        const size_t nGhost = trim->numGhostFaces();
        const bool levelPass = (bSkel == lSkel) && (bGhost == lGhost)
            && (nSkel == lSkel.size()) && (nGhost == lGhost.size());
        allPass = allPass && levelPass;

        const real_t ratio = (nSkel > 0) ? (real_t)nGhost / (real_t)nSkel : 0.0;

        gsInfo << std::setw(4) << r << std::setw(6) << N
               << std::setw(8) << nInterior << std::setw(8) << nCut
               << std::setw(8) << nSkel << std::setw(8) << nGhost
               << std::setw(12) << std::fixed << std::setprecision(4) << ratio
               << std::setw(10) << bSkel.size() << std::setw(11) << bGhost.size()
               << std::setw(8) << (levelPass ? "PASS" : "FAIL") << "\n";

        // Sliver-cell classification (printed at every level on the checked
        // path, asserted only at the finest level -- see doSliverCheck above).
        if (doSliverCheck)
        {
            index_t idx[3] = {0,0,0};
            idx[0] = (index_t)math::round(xStar / h);
            idx[1] = (index_t)std::floor(c0 / h);
            if (3 == d) idx[2] = (index_t)std::floor(c0 / h);
            size_t flat = 0, s = 1;
            for (short_t j = 0; j != d; ++j) { flat += (size_t)idx[j]*s; s *= (size_t)N; }

            gsInfo << "        sliver cell idx=(" << idx[0] << "," << idx[1];
            if (3 == d) gsInfo << "," << idx[2];
            gsInfo << ")  sign=" << sgn[flat]
                   << "  delta=" << std::scientific << std::setprecision(4) << delta
                   << "  delta/h_r=" << (delta / h) << std::fixed << std::setprecision(6) << "\n";

            if (r == cfg.numRefine)
                GISMO_ENSURE(0 == sgn[flat],
                    "The near-degenerate cell was not classified as cut at the finest level.");
        }

        // ------------------------------------------------------------------
        // Sum [[d^k u]]^2 over the skeleton and over the ghost face set, for
        // a deterministic "random" spline. getSolution()/jump() exist only
        // through an assembler (a function-backed gsFeVariable does not
        // carry jump()/avg()), which is why initSystem() is called here even
        // though nothing is ever assembled: it only allocates and finalizes
        // the dof mapper that numDofs() and getSolution() require.
        // ------------------------------------------------------------------
        gsBoundaryConditions<real_t> bc; bc.setGeoMap(mp);
        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(dbasis);
        A.setIntegrationDomain(trim);                     // AFTER setIntegrationElements
        auto G = A.getMap(mp);
        auto u = A.getSpace(dbasis);
        u.setup(bc, dirichlet::interpolation, 0);
        A.initSystem();

        gsMatrix<real_t> coefs(A.numDofs(), 1);            // NAMED: getSolution keeps a pointer to this
        for (index_t i = 0; i != coefs.rows(); ++i)
            coefs(i,0) = math::sin(1.7*(real_t)i + 0.3);
        auto sol = A.getSolution(u, coefs);

        gsExprEvaluator<real_t> ev(A);                     // shares A's expression data AND A's domain
        // Uniform knot insertion gives maximal continuity C^{degree-1}: for
        // k < degree the [[d^k u]] jump is IDENTICALLY zero in exact
        // arithmetic. What is printed there is not that zero: the face loop
        // never evaluates both sides exactly ON the face -- gsExprEvaluator.h:
        // 746-757 shifts points() by -faceShift*h and pointsIfc() by
        // +faceShift*h, so the quantity actually computed is the one-sided
        // difference d^k u(x+faceShift*h) - d^k u(x-faceShift*h), a
        // deterministic O(faceShift*h) evaluation offset that scales as
        // (faceShift*h)^2 once squared -- not roundoff. The k == degree row
        // is the first genuinely discontinuous jump, and dnk(.,.,k) for
        // k > degree returns silent exact zeros (the (degree+1)-th
        // derivative of a degree-p polynomial piece vanishes). Capped at 2
        // because that is as far as this example's default degree=2
        // exercises a real signal.
        const index_t kMax = math::min(cfg.degree, (index_t)2);
        for (index_t k = 1; k <= kMax; ++k)
        {
            const real_t Jskel = ev.integralSkeleton(
                  dnk(sol.jump(), G.left(), k) * dnk(sol.jump(), G.left(), k) );
            const real_t Jghost = ev.integralGhost(
                  dnk(sol.jump(), G.left(), k) * dnk(sol.jump(), G.left(), k) );

            // Non-vacuity, gated on a non-empty skeleton: a configuration
            // with no skeleton face at all legitimately integrates to zero.
            if (nSkel > 0)
                GISMO_ENSURE(Jskel > 0, "Skeleton jump integral vanished on a non-empty skeleton.");
            // Ghost is a subset of skeleton with a non-negative integrand.
            GISMO_ENSURE(Jghost <= Jskel * (1 + 1e-10),
                "Ghost jump integral exceeds the skeleton jump integral.");

            gsInfo << "        k=" << k << "   sum [[d^k u]]^2 :  skeleton = "
                   << std::scientific << std::setprecision(3) << Jskel
                   << "   ghost = " << Jghost << std::fixed << "\n";
        }

        // ------------------------------------------------------------------
        // ParaView export at the finest level only: skeleton faces, ghost
        // faces (distinct field value) and the interior/cut element signs.
        // Faces are degenerate boxes (zero extent in side().direction()), so
        // the box-export idiom below handles both without special-casing.
        // ------------------------------------------------------------------
        if (cfg.plot && r == cfg.numRefine)
        {
            auto appendBox = [](gsMatrix<real_t> & boxes, const gsVector<real_t> & lo, const gsVector<real_t> & hi)
            {
                const index_t c = boxes.cols();
                boxes.conservativeResize(lo.rows(), c + 2);
                boxes.col(c)     = lo;
                boxes.col(c + 1) = hi;
            };

            gsMatrix<real_t> boxSkel(d,0), boxGhost(d,0), boxInt(d,0), boxCut(d,0);

            gsDomain<real_t>::iterator endSkelPlot = trim->endSkeleton();
            for (gsDomain<real_t>::iterator it = trim->beginSkeleton(); it < endSkelPlot; ++it)
                appendBox(boxSkel, it.lowerCorner(), it.upperCorner());

            gsDomain<real_t>::iterator endGhostPlot = trim->endGhost();
            for (gsDomain<real_t>::iterator it = trim->beginGhost(); it < endGhostPlot; ++it)
                appendBox(boxGhost, it.lowerCorner(), it.upperCorner());

            gsDomain<real_t>::iterator endIntPlot = trim->template end<InteriorSign>();
            for (gsDomain<real_t>::iterator it = trim->beginInterior(); it < endIntPlot; ++it)
                appendBox(boxInt, it.lowerCorner(), it.upperCorner());

            gsDomain<real_t>::iterator endCutPlot = trim->endBdr(boundary::none);
            for (gsDomain<real_t>::iterator it = trim->beginBdr(boundary::none); it < endCutPlot; ++it)
                appendBox(boxCut, it.lowerCorner(), it.upperCorner());

            if (boxSkel.cols()  > 0) gsWriteParaview(boxSkel,  out + "/skeleton_faces",     real_t(1));
            if (boxGhost.cols() > 0) gsWriteParaview(boxGhost, out + "/ghost_faces",        real_t(2));
            if (boxInt.cols()   > 0) gsWriteParaview(boxInt,   out + "/elements_interior",  real_t(1));
            if (boxCut.cols()   > 0) gsWriteParaview(boxCut,   out + "/elements_cut",       real_t(2));

            gsInfo << "Wrote ParaView output to " << out
                   << " (skeleton_faces, ghost_faces, elements_interior, elements_cut)\n";
        }

        if (r != cfg.numRefine) dbasis.uniformRefine();    // refine at the END of the iteration
    }

    return allPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
