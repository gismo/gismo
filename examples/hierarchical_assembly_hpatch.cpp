/** @file hierarchical_assembly_hpatch.cpp

    @brief Assembly test for gsHPatchRule: mass, stiffness, and biharmonic
    matrices assembled with gsExprAssembler using three quadrature rules and
    compared against element-wise Gauss quadrature.

    The three rules (selected through the gsExprAssembler `quRule` option) are:
      1 = Gauss-Legendre (reference),
      3 = gsPatchRule    (hierarchical basis -> defaults to the finest level),
      4 = gsHPatchRule   (per-atomic-cell hierarchical patch rule).

    For each bilinear form the matrix is assembled with all three rules; the
    error reported is ||M_rule - M_gauss||_inf together with the total number
    of quadrature points and the assembly time.

    Forms:
      MASS       :  u * u^T
      STIFFNESS  :  grad(u) . grad(u)^T
      BIHARMONIC :  lapl(u) * lapl(u)^T   (requires regularity >= 2)

    The geometry is the identity map on the parameter domain, so meas(G)=1 and
    the assembled forms are the parametric mass/stiffness/biharmonic matrices.

    Usage:
        ./hierarchical_assembly_hpatch [-P deg] [-R reg] [-l levels] [-d dim] [--plot]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAssembler/gsHPatchRule.h>
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>

using namespace gismo;

// ============================================================
// Enum for bilinear form type
// ============================================================
enum FormType { MASS, STIFFNESS, BIHARMONIC };
enum MeshPattern { PatternLShape, PatternCenter, PatternDiagonal };

// ============================================================
// Helpers
// ============================================================

template<class T>
T sparseDiffMaxAbs(const gsSparseMatrix<T> & lhs, const gsSparseMatrix<T> & rhs)
{
    GISMO_ASSERT(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols(),
                 "sparseDiffMaxAbs: matrix sizes do not match");

    T maxAbs = T(0);
    const index_t outerCount = lhs.outerSize();
    for (index_t outer = 0; outer < outerCount; ++outer)
    {
        typename gsSparseMatrix<T>::InnerIterator itL(lhs, outer);
        typename gsSparseMatrix<T>::InnerIterator itR(rhs, outer);
        while (itL || itR)
        {
            if (!itR || (itL && itL.index() < itR.index()))
            { maxAbs = std::max(maxAbs, math::abs(itL.value())); ++itL; }
            else if (!itL || itR.index() < itL.index())
            { maxAbs = std::max(maxAbs, math::abs(itR.value())); ++itR; }
            else
            { maxAbs = std::max(maxAbs, math::abs(itL.value() - itR.value())); ++itL; ++itR; }
        }
    }
    return maxAbs;
}

template<class T>
T sparseMaxAbs(const gsSparseMatrix<T> & m)
{
    T mx = T(0);
    for (index_t k = 0; k < m.outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(m, k); it; ++it)
            mx = std::max(mx, math::abs(it.value()));
    return mx;
}

/// Total number of quadrature points produced by the rule selected in \a opts
/// over the domain of \a basis (the rule is built once and mapped per element).
template<class T>
index_t countPoints(const gsBasis<T> & basis, const gsOptionList & opts)
{
    typename gsQuadRule<T>::uPtr rule = gsQuadrature::getPtr(basis, opts);
    index_t n = 0;
    gsMatrix<T> nodes; gsVector<T> weights;
    for (auto & elem : basis.domain()->allElements())
    {
        rule->mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        n += weights.size();
    }
    return n;
}

/// Assemble one bilinear \a form. The DoFs are always defined by \a spaceMB
/// (the original basis); the quadrature rule and element iteration use
/// \a integMB as the integration domain.
template<short_t d, class T>
gsSparseMatrix<T> assembleForm(const gsMultiBasis<T> & spaceMB,
                               const gsMultiBasis<T> & integMB,
                               const gsMultiPatch<T> & mp,
                               FormType                form,
                               const gsOptionList    & opts,
                               T                     & time)
{
    gsExprAssembler<T> A(1,1);
    A.setOptions(opts);
    A.setIntegrationElements(integMB);   // drives quadrature + element iteration

    auto G = A.getMap(mp);
    auto u = A.getSpace(spaceMB);        // drives the DoFs

    gsStopwatch sw;
    A.initSystem();
    switch (form)
    {
        case MASS:       A.assemble( u * u.tr() * meas(G) ); break;
        case STIFFNESS:  A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) ); break;
        case BIHARMONIC: A.assemble( ilapl(u,G) * ilapl(u,G).tr() * meas(G) ); break;
    }
    time = sw.stop();
    return A.matrix();
}

/// Assemble the load vector b_i = ∫ N_i dΩ (integral of each basis function).
/// spaceMB defines the DoFs; integMB drives the quadrature rule and iteration.
template<short_t d, class T>
gsVector<T> assembleRHS(const gsMultiBasis<T> & spaceMB,
                        const gsMultiBasis<T> & integMB,
                        const gsMultiPatch<T> & mp,
                        const gsOptionList    & opts)
{
    gsExprAssembler<T> A(1,1);
    A.setOptions(opts);
    A.setIntegrationElements(integMB);

    auto G = A.getMap(mp);
    auto u = A.getSpace(spaceMB);

    // gsExprAssembler::assemble() expects the sparse system to be initialized,
    // even when we only consume the assembled right-hand side.
    A.initSystem();
    A.assemble( u * meas(G) );   // b_i = ∫ N_i dΩ
    return A.rhs();
}

bool parseMeshPattern(const std::string & pattern, MeshPattern & result)
{
    if (pattern == "lshape")
    {
        result = PatternLShape;
        return true;
    }
    if (pattern == "center")
    {
        result = PatternCenter;
        return true;
    }
    if (pattern == "diagonal")
    {
        result = PatternDiagonal;
        return true;
    }
    return false;
}

template<short_t d>
void appendTargetBox(const typename gsHTensorBasis<d,real_t>::point & low,
                     const typename gsHTensorBasis<d,real_t>::point & upp,
                     index_t targetLevel,
                     std::vector<index_t> & boxes)
{
    boxes.push_back(targetLevel);
    for (short_t dd = 0; dd < d; ++dd)
        boxes.push_back(low[dd]);
    for (short_t dd = 0; dd < d; ++dd)
        boxes.push_back(upp[dd]);
}

template<short_t d>
std::vector<index_t> makePatternBoxes(index_t targetLevel, MeshPattern pattern)
{
    const index_t N = 4 * (1 << targetLevel);
    const index_t width = 4;
    std::vector<index_t> boxes;
    typename gsHTensorBasis<d,real_t>::point low, upp;

    switch (pattern)
    {
        case PatternLShape:
            boxes.reserve(d * (2 * d + 1));
            for (short_t DD = 0; DD < d; ++DD)
            {
                low.setZero();
                for (short_t dd = 0; dd < d; ++dd)
                    upp[dd] = (dd == DD ? N : width);
                appendTargetBox<d>(low, upp, targetLevel, boxes);
            }
            break;

        case PatternCenter:
            boxes.reserve(2 * d + 1);
            for (short_t dd = 0; dd < d; ++dd)
            {
                low[dd] = (N - width) / 2;
                upp[dd] = low[dd] + width;
            }
            appendTargetBox<d>(low, upp, targetLevel, boxes);
            break;

        case PatternDiagonal:
            // Offline construction cost: O(2^targetLevel) boxes in every
            // dimension, still negligible relative to assembly.
            boxes.reserve((N / width) * (2 * d + 1));
            for (index_t start = 0; start < N; start += width)
            {
                for (short_t dd = 0; dd < d; ++dd)
                {
                    low[dd] = start;
                    upp[dd] = start + width;
                }
                appendTargetBox<d>(low, upp, targetLevel, boxes);
            }
            break;
    }

    return boxes;
}

template<short_t d, class T>
typename gsHElementHelper<d,T>::HElementContainer
makeMarkedElementsFromBoxes(const std::vector<index_t> & boxes,
                            index_t targetLevel,
                            const gsHElementHelper<d,T> & helper)
{
    typedef gsHElementHelper<d,T> helper_t;
    typedef typename helper_t::HElementContainer HElementContainer;
    typedef typename helper_t::point_t point_t;

    GISMO_ASSERT(targetLevel > 0, "Target refinement level must be positive.");
    GISMO_ASSERT(boxes.size() % (2 * d + 1) == 0, "Invalid refinement box encoding.");

    HElementContainer marked;
    const index_t parentLevel = targetLevel - 1;
    const index_t stride = 2 * d + 1;

    for (size_t i = 0; i < boxes.size() / stride; ++i)
    {
        const size_t offset = i * stride;
        GISMO_ASSERT(boxes[offset] == targetLevel,
                     "Expected box level " << targetLevel << ", got " << boxes[offset]);

        point_t low, upp;
        for (short_t dd = 0; dd < d; ++dd)
        {
            GISMO_ASSERT(boxes[offset + 1 + dd] % 2 == 0 &&
                         boxes[offset + 1 + d + dd] % 2 == 0,
                         "Pattern boxes must align with parent elements.");
            low[dd] = boxes[offset + 1 + dd] / 2;
            upp[dd] = boxes[offset + 1 + d + dd] / 2;
        }

        typename helper_t::HElementContainer local = helper.explode(low, upp, parentLevel);
        marked.insert(local.begin(), local.end());
    }

    return marked;
}

template<short_t d, class T>
void refineLevel(gsTHBSplineBasis<d,T> & thb,
                 MeshPattern pattern,
                 index_t targetLevel,
                 bool admissible,
                 index_t admissibleJump)
{
    const std::vector<index_t> boxes = makePatternBoxes<d>(targetLevel, pattern);
    if (!admissible)
    {
        thb.refineElements(boxes);
        return;
    }

    gsHElementHelper<d,T> helper(thb);
    typename gsHElementHelper<d,T>::HElementContainer marked =
        makeMarkedElementsFromBoxes<d,T>(boxes, targetLevel, helper);

    // Admissible closure is mesh-construction work only; the extra cost scales
    // with the number of marked elements and their hierarchical neighborhoods.
    marked = helper.markAdmissible(marked, admissibleJump);
    thb.refineElements(helper.toRefBoxes(marked, false));
}

// ============================================================
// Compare the three rules for one bilinear form
// ============================================================

template<short_t d, class T>
void compareRules(const std::string     & formName,
                  const gsMultiBasis<T> & spaceMB,
                  const gsMultiBasis<T> & integMB,
                  const gsMultiPatch<T> & mp,
                  FormType                form,
                  index_t                 quA,
                  index_t                 quB)
{
    auto makeOpts = [&](index_t quRule)
    {
        gsOptionList o = gsExprAssembler<T>::defaultOptions();
        o.setInt ("quRule", quRule);
        o.setReal("quA", static_cast<T>(quA));
        o.setInt ("quB", quB);
        return o;
    };

    gsOptionList optGauss = gsExprAssembler<T>::defaultOptions();
    optGauss.setInt ("quRule", gsQuadrature::GaussLegendre);
    optGauss.setReal("quA", static_cast<T>(1.0));
    optGauss.setInt ("quB", 1);
    // gsOptionList optPatch  = makeOpts(gsQuadrature::PatchRule);
    gsOptionList optHPatch = makeOpts(gsQuadrature::HPatchRule);

    // Gauss and gsPatchRule integrate on the original domain (they do not need
    // the refined integration basis); gsHPatchRule uses integMB.
    T tGauss = 0, tPatch = 0, tHPatch = 0;
    const gsSparseMatrix<T> Mgauss  = assembleForm<d,T>(spaceMB, spaceMB, mp, form, optGauss,  tGauss);
    // const gsSparseMatrix<T> Mpatch  = assembleForm<d,T>(spaceMB, spaceMB, mp, form, optPatch,  tPatch);
    const gsSparseMatrix<T> Mhpatch = assembleForm<d,T>(spaceMB, integMB, mp, form, optHPatch, tHPatch);

    const index_t nGauss  = countPoints(spaceMB.basis(0), optGauss);
    // const index_t nPatch  = countPoints(spaceMB.basis(0), optPatch);
    const index_t nHPatch = countPoints(integMB.basis(0), optHPatch);

    // Reference magnitude for relative error (biharmonic entries are huge,
    // ~h^-4, so absolute error is not a meaningful exactness measure there).
    const T magGauss  = sparseMaxAbs(Mgauss);
    // const T errPatch  = sparseDiffMaxAbs(Mpatch,  Mgauss);
    const T errHPatch = sparseDiffMaxAbs(Mhpatch, Mgauss);

    auto row = [&](const std::string & name, index_t npts, T err, T ratio, T tm)
    {
        gsInfo << std::left
               << std::setw(12) << name
               << std::setw(12) << npts
               << std::setw(14) << std::scientific << std::setprecision(3) << err
               << std::setw(14) << std::scientific << std::setprecision(3)
                                << (magGauss > 0 ? err / magGauss : err)
               << std::setw(9)  << std::fixed << std::setprecision(3) << ratio
               << std::setw(10) << std::fixed << std::setprecision(4) << tm
               << "\n";
    };

    gsInfo << "[" << formName << "]   (||M_gauss||_inf = "
           << std::scientific << std::setprecision(2) << magGauss << ")\n";
    gsInfo << std::left
           << std::setw(12) << "Rule" << std::setw(12) << "nPts"
           << std::setw(14) << "abs.err" << std::setw(14) << "rel.err"
           << std::setw(9)  << "ratio" << std::setw(10) << "Time (s)" << "\n";
    row("Gauss",        nGauss,  T(0),      T(1),                          tGauss);
    // row("gsPatchRule",  nPatch,  errPatch,  (T)nPatch  / nGauss,           tPatch);
    row("gsHPatchRule", nHPatch, errHPatch, (T)nHPatch / nGauss,           tHPatch);

    // --- RHS diagnostic: b_i = ∫ N_i dΩ ---
    // Reveals per-function integration errors more directly than the matrix.
    // For the identity geometry, sum(b) = domain volume = 1 (partition of unity).
    // Truncated THB functions that are piecewise-polynomial within an atomic cell
    // would show non-zero relative error here if the rule misses their kink points.
    {
        const gsVector<T> bGauss  = assembleRHS<d,T>(spaceMB, spaceMB, mp, optGauss);
        // const gsVector<T> bPatch  = assembleRHS<d,T>(spaceMB, spaceMB, mp, optPatch);
        const gsVector<T> bHPatch = assembleRHS<d,T>(spaceMB, integMB, mp, optHPatch);

        const T rhsMag  = bGauss.cwiseAbs().maxCoeff();
        // const T rhsErrP = (bPatch  - bGauss).cwiseAbs().maxCoeff();
        const T rhsErrH = (bHPatch - bGauss).cwiseAbs().maxCoeff();

        gsInfo << "  ∫N_i dΩ — sum(Gauss)=" << std::fixed << std::setprecision(15)
               << bGauss.sum()
               << "  sum(hpatch)=" << bHPatch.sum() << "\n";
        // gsInfo << "  gsPatchRule  rel.err=" << std::scientific << std::setprecision(3)
        //        << (rhsMag > 0 ? rhsErrP / rhsMag : rhsErrP);
        gsInfo << "  gsHPatchRule rel.err="
               << (rhsMag > 0 ? rhsErrH / rhsMag : rhsErrH) << "\n\n";
    }
}

// ============================================================
// Templated main
// ============================================================

template<short_t d, class T>
void run(bool plot,
         index_t numLevels,
         index_t p,
         index_t r,
         bool useIntegrationBasis,
         MeshPattern pattern,
         bool admissible,
         index_t admissibleJump)
{
    // ------------------------------------------------------------------
    // Build THB-spline basis with configurable local refinement
    // ------------------------------------------------------------------
    gsKnotVector<T> kv(0, 1.0, 0, p + 1, 1);
    kv.uniformRefine();
    kv.uniformRefine();

    std::vector<gsKnotVector<T>> kvs;
    for (short_t dd = 0; dd < d; ++dd) kvs.push_back(kv);

    gsTensorBSplineBasis<d,T> tbsb(kvs);
    gsTHBSplineBasis<d,T>     thb(tbsb);

    gsInfo << "Building geometry..." << std::flush;
    {
        thb.uniformRefine();
        for (index_t lvl = 1; lvl <= numLevels; ++lvl)
            refineLevel<d,T>(thb, pattern, lvl, admissible, admissibleJump);
    }
    gsInfo << " done.\n";

    // The DoFs are ALWAYS defined by the original thb. The integration basis
    // (optionally refined via makeIntegrationBasis) only serves as the
    // integration domain for gsHPatchRule.
    gsTHBSplineBasis<d,T> integThb(thb);
    if (useIntegrationBasis)
    {
        integThb = gsHPatchRule<T>::template makeIntegrationBasis<d>(thb);
        gsInfo << "Integration basis: refined from " << thb.size()
               << " to " << integThb.size() << " DOF ("
               << integThb.numElements() << " elements)\n";
    }

    const char * patternName =
        pattern == PatternLShape ? "lshape" :
        pattern == PatternCenter ? "center" : "diagonal";
    gsInfo << "Basis: degree=" << p << "  regularity=" << r
           << "  pattern=" << patternName
           << "  admissible=" << (admissible ? "true" : "false")
           << "  jump=" << admissibleJump
           << "  maxLevel=" << thb.maxLevel() << "  DOF=" << thb.size()
           << "  elements=" << thb.numElements() << "\n\n";

    if (plot && (d == 2 || d == 3))
    {
        gsMesh<T> mesh(thb);
        gsWriteParaview(mesh, "hierarchical_assembly_mesh");
        gsOptionList opts = gsExprAssembler<T>::defaultOptions();
        opts.setInt ("quRule", gsQuadrature::HPatchRule);
        opts.setReal("quA", static_cast<T>(2*p));
        opts.setInt ("quB", r);
        gsMatrix<T> nodes = gsQuadrature::getAllNodes(thb,opts);
        gsInfo<<"Plotting quadrature nodes ("<<nodes.rows()<<" points)... "<<std::flush;
        gsWriteParaviewPoints(nodes, "hierarchical_assembly_nodes");
        gsInfo<<" done.\n";
    }
    if (plot && d==1)
    {
        gsWriteParaview(thb, "hierarchical_assembly_basis", 1000, true);
        gsOptionList opts = gsExprAssembler<T>::defaultOptions();
        opts.setInt ("quRule", gsQuadrature::HPatchRule);
        opts.setReal("quA", static_cast<T>(2*p));
        opts.setInt ("quB", r);
        gsMatrix<T> nodes = gsQuadrature::getAllNodes(thb,opts);
        gsInfo<<"Plotting quadrature nodes ("<<nodes.rows()<<" points)... "<<std::flush;
        gsWriteParaviewPoints(nodes, "hierarchical_assembly_nodes");
        gsInfo<<" done.\n";
    }

    // ------------------------------------------------------------------
    // Identity geometry (so meas(G) = 1): coefs = anchors of coarsest level
    // ------------------------------------------------------------------
    gsTensorBSplineBasis<d,T> & l0 = thb.tensorLevel(0);
    gsMatrix<T> coefs = l0.anchors().transpose();         // (N x d) identity coefs
    gsTensorBSpline<d,T> geom(l0, coefs);

    gsMultiPatch<T> mp;
    mp.addPatch(geom);

    gsMultiBasis<T> spaceMB;                       // DoFs: original thb
    spaceMB.addBasis(thb.clone().release());
    spaceMB.setTopology(mp);

    gsMultiBasis<T> integMB;                        // integration domain for hpatch
    integMB.addBasis(integThb.clone().release());
    integMB.setTopology(mp);

    // ------------------------------------------------------------------
    // Assemble and compare (quA = 2p target degree; quB = target regularity)
    // ------------------------------------------------------------------
    compareRules<d,T>("Mass",      spaceMB, integMB, mp, MASS,      2*p, r);
    compareRules<d,T>("Stiffness", spaceMB, integMB, mp, STIFFNESS, 2*p, r-1);

    if (r >= 2)
        compareRules<d,T>("Biharmonic", spaceMB, integMB, mp, BIHARMONIC, 2*p, r-2);
    else
        gsInfo << "Biharmonic: skipped (requires regularity >= 2; use -P 3 -R 2)\n";

    // gsHPatchRule diagnostics (fallback cells, on the integration domain)
    {
        gsHPatchRule<T> hp(*integThb.domain(), 2*p, r, false);
        gsInfo << "gsHPatchRule: " << hp.numPointsTotal() << " points, "
               << hp.numFallbackCells() << " fallback cells\n";
    }
}

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    index_t p          = 2;
    index_t r          = 1;
    index_t numLevels  = 3;
    index_t dim        = 2;
    index_t admissibleJump = 2;
    bool    plot       = false;
    bool    integBasis = false;
    bool    admissible = false;
    std::string meshPattern = "lshape";

    gsCmdLine cmd("Assembly test: gsHPatchRule vs Gauss and gsPatchRule (gsExprAssembler).");
    cmd.addInt   ("P", "deg",    "Polynomial degree",           p);
    cmd.addInt   ("R", "reg",    "Regularity",                  r);
    cmd.addInt   ("l", "levels", "Number of refinement levels", numLevels);
    cmd.addInt   ("d", "dim",    "Problem dimension (1..4)",    dim);
    cmd.addInt   ("",  "admissibleJump",
                  "Jump parameter m for admissible refinement", admissibleJump);
    cmd.addString("",  "meshPattern",
                  "Mesh pattern: lshape, center, diagonal", meshPattern);
    cmd.addSwitch("",  "admissible",
                  "Apply admissible closure to the marked refinement region", admissible);
    cmd.addSwitch("plot",        "Write ParaView mesh",         plot);
    cmd.addSwitch("integrationBasis",
                  "Refine THB via makeIntegrationBasis before assembly", integBasis);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (r >= p)
    {
        gsInfo << "Error: regularity must be < degree.\n";
        return 1;
    }
    if (admissibleJump < 1)
    {
        gsInfo << "Error: admissibleJump must be >= 1.\n";
        return 1;
    }

    MeshPattern pattern = PatternLShape;
    if (!parseMeshPattern(meshPattern, pattern))
    {
        gsInfo << "Error: unknown meshPattern '" << meshPattern
               << "'. Expected one of: lshape, center, diagonal.\n";
        return 1;
    }

    if      (dim == 1) run<1, real_t>(plot, numLevels, p, r, integBasis, pattern, admissible, admissibleJump);
    else if (dim == 2) run<2, real_t>(plot, numLevels, p, r, integBasis, pattern, admissible, admissibleJump);
    else if (dim == 3) run<3, real_t>(plot, numLevels, p, r, integBasis, pattern, admissible, admissibleJump);
    else if (dim == 4) run<4, real_t>(plot, numLevels, p, r, integBasis, pattern, admissible, admissibleJump);
    else { gsInfo << "Error: dimension must be 1..4.\n"; return 1; }

    return 0;
}
