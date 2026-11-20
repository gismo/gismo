/** @file hierarchical_assembly_hpatch.cpp

    @brief Assembly benchmark for gsHPatchRule: mass, stiffness, and biharmonic
    matrices assembled with gsExprAssembler using gsHPatchRule (reduced patch
    quadrature: one tensor patch rule per kd leaf, or per merged box of leaves
    of one level) and compared against element-wise Gauss-Legendre quadrature.

    For each bilinear form the matrix is assembled with Gauss and with
    gsHPatchRule; the error reported is ||M_rule - M_gauss||_inf together with
    the total number of quadrature points and the timings.

    The Gauss rule is built once, before the assembly. The gsHPatchRule is built
    by gsHPatchRule::makeFactory, whose factory hands out one handle (a copy that
    shares the cell data and has its own lookup position) per call. Both are
    given to gsExprAssembler through its quadrature factory
    (setQuadratureFactory), each behind a forwarding wrapper, so both pay the
    same per-element dispatch.
    Four timings are kept apart (seconds):
      tMesh  : construction of the THB-spline mesh (once per run),
      tBuild : construction of the quadrature rule; for hpatch it is the whole
               makeFactory call, i.e. the clone of the integration basis, the
               partition (--merge) and the rule itself,
      tInit  : initSystem() and computation of the sparsity pattern,
      tAsm   : one assemble() call on the initialised system, minimum over
               --reps repetitions; Gauss and hpatch alternate within the loop.

    Forms:
      MASS       :  u * u^T
      STIFFNESS  :  grad(u) . grad(u)^T
      BIHARMONIC :  lapl(u) * lapl(u)^T   (requires regularity >= 2)

    The geometry is the identity map on the parameter domain, so meas(G)=1 and
    the assembled forms are the parametric mass/stiffness/biharmonic matrices.

    With --csv one row per (form, rule) is written, columns:
      form,dim,p,r,levels,pattern,admissible,jump,nElem,nDof,maxLevel,rule,
      nPts,nFallback,absErr,relErr,tMesh,tBuild,tInit,tAsm,reps,threads,host,cpu,baseRefs,band,
      merge,tPartition,numLeaves,numCells,levelsPairwise,nodesLeaves,nodesPairwise,nodesSweep
    The last eight columns describe the partition of the hpatch rule and are
    empty in the gauss row:
      merge          : the merge mode (--merge),
      tPartition     : wall time of the candidates and the choice, in seconds
                       (0 for none); it is part of tBuild,
      numLeaves      : kd leaves of the domain,
      numCells       : cells of the rule (boxes after merging),
      levelsPairwise : levels on which the pairwise candidate was kept (the
                       non-empty levels for pairwise, 0 for sweep),
      nodesLeaves, nodesPairwise, nodesSweep : model node counts of the kd
                       leaves, the pairwise candidate and the sweep candidate;
                       -1 means not computed. The counts assume that every 1D
                       solve converges, so nPts equals the count of the
                       partition kept only when nFallback is 0.
    nFallback is the number of cells with a direction on Gauss nodes because its
    Newton solve failed.

    The timings are serial by design (the benchmark compares single-thread
    assembly times), so the driver must be run with OMP_NUM_THREADS=1 (every
    run, including --plot); otherwise it exits with code 2.

    Usage:
        OMP_NUM_THREADS=1 ./hierarchical_assembly_hpatch [-P deg] [-R reg] [-l levels]
            [-d dim] [--form mass|stiffness|biharmonic|all] [--reps n] [--csv file]
            [--meshPattern lshape|center|diagonal] [--admissible] [--admissibleJump m]
            [--baseRefs k] [--band scaled|thin] [--merge none|pairwise|sweep|best]
            [--integrationBasis] [--plot]
    The regularity -R defaults to deg-1. --merge selects the leaf merging of the
    hpatch rule (default best, the library default); --plot always shows the
    nodes of the best rule.

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

// Width of the refinement bands: BandScaled keeps the refined region's physical
// size for every --baseRefs >= 1 (the centre box is shifted to an even origin,
// which moves it by one level-targetLevel element at --baseRefs 1); BandThin
// keeps the bands 4 target-level elements wide.
enum BandVariant { BandScaled, BandThin };

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

/// Total number of quadrature points that assemble() evaluates when it maps
/// \a rule onto every element of the domain of \a basis.
template<class T>
index_t countPoints(const gsBasis<T> & basis, const gsQuadRule<T> & rule)
{
    index_t n = 0;
    gsMatrix<T> nodes; gsVector<T> weights;
    for (auto & elem : basis.domain()->allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        n += weights.size();
    }
    return n;
}

/// Quadrature rule that forwards mapTo() to another rule.
///
/// gsExprAssembler asks its quadrature factory for a fresh rule per patch and
/// thread; this wrapper lets the factory hand out rules that were built (and
/// timed) beforehand, behind one layer of virtual dispatch for every arm. It
/// forwards either to a shared stateless rule (Gauss) or to a handle that
/// belongs to this wrapper alone (gsHPatchRule). numNodes() and dim() of the
/// base class are not forwarded.
template<class T>
class gsForwardingQuadRule : public gsQuadRule<T>
{
public:
    explicit gsForwardingQuadRule(memory::shared_ptr<const gsQuadRule<T> > rule)
    : m_rule(give(rule)) { }

    using gsQuadRule<T>::mapTo;

    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    { m_rule->mapTo(lower, upper, nodes, weights); }

    void mapTo(T startVal, T endVal,
               gsMatrix<T>& nodes, gsVector<T>& weights) const override
    { m_rule->mapTo(startVal, endVal, nodes, weights); }

private:
    memory::shared_ptr<const gsQuadRule<T> > m_rule;
};

/// Assembler together with the geometry map and trial space it hands out, set
/// up to integrate with a prebuilt quadrature rule or with a rule factory.
///
/// The DoFs are defined by \a spaceMB (the original basis); element iteration
/// uses \a integMB as the integration domain.
template<class T>
struct AssemblyArm
{
    typedef typename gsExprAssembler<T>::QuadratureRulePtr QuadratureRulePtr;
    typedef typename gsExprAssembler<T>::QuadratureFactory QuadratureFactory;

    gsExprAssembler<T>                      A;
    typename gsExprAssembler<T>::geometryMap G;
    typename gsExprAssembler<T>::space       u;

    AssemblyArm(const gsMultiBasis<T> & spaceMB,
                const gsMultiBasis<T> & integMB,
                const gsMultiPatch<T> & mp,
                memory::shared_ptr<const gsQuadRule<T> > rule)
    : A(1,1), G(A.getMap(mp)), u(A.getSpace(spaceMB))
    {
        A.setIntegrationElements(integMB);
        A.setQuadratureFactory(
            [rule](const gsBasis<T> & basis, const gsOptionList & options,
                   index_t patch, short_t fixedDirection) -> QuadratureRulePtr
            {
                GISMO_UNUSED(basis); GISMO_UNUSED(options); GISMO_UNUSED(patch);
                GISMO_ENSURE(fixedDirection == -1,
                             "The forwarding rule serves volume integration only.");
                return QuadratureRulePtr(new gsForwardingQuadRule<T>(rule));
            });
    }

    /// Every call of \a inner returns a rule that the wrapper takes over.
    AssemblyArm(const gsMultiBasis<T> & spaceMB,
                const gsMultiBasis<T> & integMB,
                const gsMultiPatch<T> & mp,
                QuadratureFactory inner)
    : A(1,1), G(A.getMap(mp)), u(A.getSpace(spaceMB))
    {
        A.setIntegrationElements(integMB);
        A.setQuadratureFactory(
            [inner](const gsBasis<T> & basis, const gsOptionList & options,
                    index_t patch, short_t fixedDirection) -> QuadratureRulePtr
            {
                GISMO_ENSURE(fixedDirection == -1,
                             "The forwarding rule serves volume integration only.");
                QuadratureRulePtr result = inner(basis, options, patch, fixedDirection);
                return QuadratureRulePtr(new gsForwardingQuadRule<T>(
                    memory::shared_ptr<const gsQuadRule<T> >(result.release())));
            });
    }

private:
    AssemblyArm(const AssemblyArm &);
    AssemblyArm & operator=(const AssemblyArm &);
};

/// Timings of one arm, in seconds.
struct ArmTimes
{
    double init;   ///< initSystem() and sparsity pattern
    double assem;  ///< minimum over the repetitions of one assemble() call
};

/// Assemble the bilinear form \a e (same expression in every call) \a reps times
/// on each of the two arms, alternating between them.
///
/// The system is initialised and its pattern computed once per arm (tInit).
/// Before each repetition the values are zeroed in place with clearMatrix(),
/// which keeps the pattern; only the assemble() call is timed. Alternating the
/// arms spreads drifts of the machine (frequency, caches) over both.
template<class T, class ExprA, class ExprB>
void timeAssembly(AssemblyArm<T> & a, ExprA ea,
                  AssemblyArm<T> & b, ExprB eb,
                  index_t reps, ArmTimes & ta, ArmTimes & tb)
{
    gsStopwatch sw;
    a.A.initSystem();
    a.A.computePattern(ea);
    ta.init = sw.stop();
    sw.restart();
    b.A.initSystem();
    b.A.computePattern(eb);
    tb.init = sw.stop();

    ta.assem = tb.assem = std::numeric_limits<double>::max();
    for (index_t k = 0; k < reps; ++k)
    {
        a.A.clearMatrix();
        sw.restart();
        a.A.assemble(ea);
        ta.assem = std::min(ta.assem, sw.stop());

        b.A.clearMatrix();
        sw.restart();
        b.A.assemble(eb);
        tb.assem = std::min(tb.assem, sw.stop());
    }
}

/// Selects the bilinear form and runs timeAssembly() on the two arms.
template<class T>
void timeForm(FormType form, AssemblyArm<T> & a, AssemblyArm<T> & b,
              index_t reps, ArmTimes & ta, ArmTimes & tb)
{
    switch (form)
    {
        case MASS:
            timeAssembly(a, a.u * a.u.tr() * meas(a.G),
                         b, b.u * b.u.tr() * meas(b.G), reps, ta, tb);
            break;
        case STIFFNESS:
            timeAssembly(a, igrad(a.u,a.G) * igrad(a.u,a.G).tr() * meas(a.G),
                         b, igrad(b.u,b.G) * igrad(b.u,b.G).tr() * meas(b.G),
                         reps, ta, tb);
            break;
        case BIHARMONIC:
            timeAssembly(a, ilapl(a.u,a.G) * ilapl(a.u,a.G).tr() * meas(a.G),
                         b, ilapl(b.u,b.G) * ilapl(b.u,b.G).tr() * meas(b.G),
                         reps, ta, tb);
            break;
    }
}

/// Assemble the load vector b_i = ∫ N_i dΩ (integral of each basis function).
/// spaceMB defines the DoFs; integMB drives the quadrature rule and iteration.
/// \a factory is the quadrature factory to integrate with; nullptr selects the
/// rule named by the option "quRule" of \a opts.
template<short_t d, class T>
gsVector<T> assembleRHS(const gsMultiBasis<T> & spaceMB,
                        const gsMultiBasis<T> & integMB,
                        const gsMultiPatch<T> & mp,
                        const gsOptionList    & opts,
                        const typename gsExprAssembler<T>::QuadratureFactory * factory)
{
    gsExprAssembler<T> A(1,1);
    A.setOptions(opts);
    A.setIntegrationElements(integMB);
    if (factory != nullptr)
        A.setQuadratureFactory(*factory);

    auto G = A.getMap(mp);
    auto u = A.getSpace(spaceMB);

    // gsExprAssembler::assemble() expects the sparse system to be initialized,
    // even when we only consume the assembled right-hand side.
    A.initSystem();
    A.assemble( u * meas(G) );   // b_i = ∫ N_i dΩ
    return A.rhs();
}

/// Whether \a name is a merge mode of gsHPatchRule: none, pairwise, sweep, best.
bool isMergeMode(const std::string & name)
{
    return name == "none" || name == "pairwise" || name == "sweep" || name == "best";
}

/// The gsHPatchRule merge mode named \a name (see isMergeMode()).
template<class T>
typename gsHPatchRule<T>::MergeMode mergeModeOf(const std::string & name)
{
    typedef typename gsHPatchRule<T>::MergeMode MergeMode;
    GISMO_ENSURE(isMergeMode(name), "unknown merge mode '" << name << "'");
    return name == "none"     ? MergeMode::none     :
           name == "pairwise" ? MergeMode::pairwise :
           name == "sweep"    ? MergeMode::sweep    : MergeMode::best;
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

bool parseBandVariant(const std::string & band, BandVariant & result)
{
    if (band == "scaled")
    {
        result = BandScaled;
        return true;
    }
    if (band == "thin")
    {
        result = BandThin;
        return true;
    }
    return false;
}

const char * bandName(BandVariant band)
{
    return band == BandScaled ? "scaled" : "thin";
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
std::vector<index_t> makePatternBoxes(index_t targetLevel, MeshPattern pattern,
                                      index_t baseRefs, BandVariant band)
{
    // K base elements per direction; N is the pattern extent in level-targetLevel
    // element indices. Scaled bands are K wide (fixed physical size), thin
    // bands stay 4 elements wide. baseRefs >= 1 makes K, N and both widths even,
    // so every box corner is a multiple of 2 and aligns with parent cells; the
    // centre origin is rounded down to an even index for the same reason.
    const index_t K = index_t(1) << baseRefs;
    const index_t N = K * (index_t(1) << targetLevel);
    const index_t width = (band == BandScaled ? K : std::min<index_t>(4, N));
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
                low[dd] = ((N - width) / 4) * 2;
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
                 index_t baseRefs,
                 BandVariant band,
                 index_t targetLevel,
                 bool admissible,
                 index_t admissibleJump)
{
    const std::vector<index_t> boxes = makePatternBoxes<d>(targetLevel, pattern, baseRefs, band);
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
// Benchmark configuration and CSV output
// ============================================================

/// Command-line settings that are shared by all rows of a run.
struct BenchConfig
{
    index_t       dim;
    index_t       p;
    index_t       r;               ///< basis regularity
    index_t       numLevels;
    MeshPattern   pattern;
    BandVariant   band;            ///< width of the refinement bands
    bool          admissible;
    index_t       admissibleJump;
    bool          plot;
    bool          useIntegrationBasis;
    std::string   form;            ///< mass, stiffness, biharmonic or all
    std::string   merge;           ///< leaf merging of the hpatch rule: none, pairwise, sweep, best
    index_t       reps;
    index_t       baseRefs;        ///< uniform refinements of the base knot vector
    std::ofstream * csv;           ///< open CSV file, or NULL
    // Filled in once the mesh exists
    index_t       nElem;
    index_t       nDof;
    index_t       maxLevel;
    double        tMesh;
};

const char * patternName(MeshPattern pattern)
{
    return pattern == PatternLShape ? "lshape" :
           pattern == PatternCenter ? "center" : "diagonal";
}

/// Header line of the CSV file.
const char * csvHeader()
{
    return "form,dim,p,r,levels,pattern,admissible,jump,nElem,nDof,maxLevel,"
           "rule,nPts,nFallback,absErr,relErr,tMesh,tBuild,tInit,tAsm,"
           "reps,threads,host,cpu,baseRefs,band,"
           "merge,tPartition,numLeaves,numCells,levelsPairwise,"
           "nodesLeaves,nodesPairwise,nodesSweep";
}

/// Partition columns of an hpatch CSV row (see the file documentation).
struct PartitionColumns
{
    PartitionColumns()
    : present(false), seconds(0), numLeaves(0), numCells(0), levelsPairwise(0),
      nodesLeaves(0), nodesPairwise(0), nodesSweep(0) { }

    bool        present;   ///< false: the eight fields are written empty
    std::string merge;
    double      seconds;
    index_t     numLeaves, numCells, levelsPairwise;
    index_t     nodesLeaves, nodesPairwise, nodesSweep;
};

/// Host name of the machine, or "unknown".
std::string hostName()
{
    std::ifstream in("/proc/sys/kernel/hostname");
    std::string line;
    if (!std::getline(in, line) || line.empty())
        return "unknown";
    return line;
}

/// Model name of the first CPU listed in /proc/cpuinfo with commas replaced by
/// blanks (so that it is a single CSV field), or "unknown".
std::string cpuModel()
{
    std::ifstream in("/proc/cpuinfo");
    std::string line;
    while (std::getline(in, line))
    {
        if (line.compare(0, 10, "model name") != 0)
            continue;
        const size_t colon = line.find(':');
        if (colon == std::string::npos)
            break;
        std::string name = line.substr(colon + 1);
        const size_t first = name.find_first_not_of(" \t");
        const size_t last  = name.find_last_not_of(" \t\r\n");
        if (first == std::string::npos)
            break;
        name = name.substr(first, last - first + 1);
        std::replace(name.begin(), name.end(), ',', ' ');
        return name;
    }
    return "unknown";
}

/// Appends one (form, rule) row to the CSV file of \a cfg, if there is one.
void writeCsvRow(const BenchConfig & cfg, const std::string & form,
                 const std::string & rule, index_t nPts, index_t nFallback,
                 double absErr, double relErr, double tBuild, double tInit,
                 double tAsm, const PartitionColumns & part)
{
    if (!cfg.csv)
        return;
    std::ofstream & os = *cfg.csv;
    os << std::scientific << std::setprecision(10);
    os << form << ',' << cfg.dim << ',' << cfg.p << ',' << cfg.r << ','
       << cfg.numLevels << ',' << patternName(cfg.pattern) << ','
       << (cfg.admissible ? 1 : 0) << ',' << cfg.admissibleJump << ','
       << cfg.nElem << ',' << cfg.nDof << ',' << cfg.maxLevel << ','
       << rule << ',' << nPts << ',' << nFallback << ','
       << absErr << ',' << relErr << ','
       << cfg.tMesh << ',' << tBuild << ',' << tInit << ',' << tAsm << ','
       << cfg.reps << ',' << omp_get_max_threads() << ','
       << hostName() << ',' << cpuModel() << ',' << cfg.baseRefs << ',' << bandName(cfg.band);
    if (part.present)
        os << ',' << part.merge << ',' << part.seconds << ',' << part.numLeaves << ','
           << part.numCells << ',' << part.levelsPairwise << ',' << part.nodesLeaves << ','
           << part.nodesPairwise << ',' << part.nodesSweep;
    else
        os << ",,,,,,,,";
    os << '\n';
}

// ============================================================
// Compare the two rules for one bilinear form
// ============================================================

template<short_t d, class T>
void compareRules(const BenchConfig     & cfg,
                  const std::string     & formName,
                  const gsMultiBasis<T> & spaceMB,
                  const gsMultiBasis<T> & integMB,
                  const gsMultiPatch<T> & mp,
                  FormType                form,
                  index_t                 quA,
                  index_t                 quB)
{
    gsOptionList optGauss = gsExprAssembler<T>::defaultOptions();
    optGauss.setInt ("quRule", gsQuadrature::GaussLegendre);
    optGauss.setReal("quA", static_cast<T>(1.0));
    optGauss.setInt ("quB", 1);
    // Options of the RHS diagnostic of the hpatch rule, which integrates with
    // the leaf factory below instead of an option-selected rule.
    const gsOptionList optDefault = gsExprAssembler<T>::defaultOptions();

    // Prebuilt rules. Gauss integrates on the original domain, gsHPatchRule on
    // the integration domain; only the construction is timed.
    gsStopwatch sw;
    double tBuildGauss = 0, tBuildHPatch = 0;
    sw.restart();
    typename gsQuadRule<T>::uPtr gaussUnique = gsQuadrature::getPtr(spaceMB.basis(0), optGauss);
    tBuildGauss = sw.stop();
    memory::shared_ptr<const gsQuadRule<T> > gaussRule(gaussUnique.release());

    sw.restart();
    const typename gsExprAssembler<T>::QuadratureFactory hpatchFactory =
        gsHPatchRule<T>::makeFactory(integMB, quA, quB, false, nullptr,
                                     mergeModeOf<T>(cfg.merge));
    tBuildHPatch = sw.stop();

    // A handle of the factory's rule for counting points and for the
    // diagnostics; outside every timed region.
    const typename gsQuadRule<T>::uPtr handle =
        hpatchFactory(integMB.basis(0), optDefault, 0, -1);
    const gsHPatchRule<T> & hpatchRule = dynamic_cast<const gsHPatchRule<T> &>(*handle);
    const typename gsHPatchRule<T>::PartitionInfo & info = hpatchRule.partitionInfo();

    AssemblyArm<T> armGauss (spaceMB, spaceMB, mp, gaussRule);
    AssemblyArm<T> armHPatch(spaceMB, integMB,  mp, hpatchFactory);

    ArmTimes tGauss, tHPatch;
    timeForm(form, armGauss, armHPatch, cfg.reps, tGauss, tHPatch);
    const gsSparseMatrix<T> Mgauss  = armGauss.A.matrix();
    const gsSparseMatrix<T> Mhpatch = armHPatch.A.matrix();

    const index_t nGauss  = countPoints(spaceMB.basis(0), *gaussRule);
    const index_t nHPatch = countPoints(integMB.basis(0), hpatchRule);
    // A merged cell hands each of its nodes to exactly one element.
    GISMO_ENSURE(nHPatch == hpatchRule.numPointsTotal(),
                 "The element-wise point count (" << nHPatch
                 << ") differs from the point count of the rule ("
                 << hpatchRule.numPointsTotal() << ").");

    // Reference magnitude for relative error (biharmonic entries are huge,
    // ~h^-4, so absolute error is not a meaningful exactness measure there).
    const T magGauss  = sparseMaxAbs(Mgauss);
    const T errHPatch = sparseDiffMaxAbs(Mhpatch, Mgauss);
    const T relHPatch = magGauss > 0 ? errHPatch / magGauss : errHPatch;

    auto row = [&](const std::string & name, index_t npts, T err, T ratio,
                   double tAsm, double tBuild, double tInit)
    {
        gsInfo << std::left
               << std::setw(12) << name
               << std::setw(12) << npts
               << std::setw(14) << std::scientific << std::setprecision(3) << err
               << std::setw(14) << std::scientific << std::setprecision(3)
                                << (magGauss > 0 ? err / magGauss : err)
               << std::setw(9)  << std::fixed << std::setprecision(3) << ratio
               << std::setw(10) << std::fixed << std::setprecision(4) << tAsm
               << std::setw(10) << std::fixed << std::setprecision(4) << tBuild
               << std::setw(10) << std::fixed << std::setprecision(4) << tInit
               << "\n";
    };

    gsInfo << "[" << formName << "]   (||M_gauss||_inf = "
           << std::scientific << std::setprecision(2) << magGauss << ")\n";
    gsInfo << std::left
           << std::setw(12) << "Rule" << std::setw(12) << "nPts"
           << std::setw(14) << "abs.err" << std::setw(14) << "rel.err"
           << std::setw(9)  << "ratio" << std::setw(10) << "Time (s)"
           << std::setw(10) << "tBuild" << std::setw(10) << "tInit" << "\n";
    row("Gauss",        nGauss,  T(0),      T(1),                 tGauss.assem,  tBuildGauss,  tGauss.init);
    row("gsHPatchRule", nHPatch, errHPatch, (T)nHPatch / nGauss,  tHPatch.assem, tBuildHPatch, tHPatch.init);
    gsInfo << "gsHPatchRule: merge=" << cfg.merge
           << " points=" << hpatchRule.numPointsTotal()
           << " fallbackCells=" << hpatchRule.numFallbackCells()
           << " leaves=" << info.numLeaves << " cells=" << info.numCells
           << " nodesLeaves=" << info.nodesLeaves
           << " nodesPairwise=" << info.nodesPairwise
           << " nodesSweep=" << info.nodesSweep
           << " levelsPairwise=" << info.levelsPairwise
           << " tPartition=" << info.seconds << "\n";

    PartitionColumns part;
    part.present        = true;
    part.merge          = cfg.merge;
    part.seconds        = info.seconds;
    part.numLeaves      = info.numLeaves;
    part.numCells       = info.numCells;
    part.levelsPairwise = info.levelsPairwise;
    part.nodesLeaves    = info.nodesLeaves;
    part.nodesPairwise  = info.nodesPairwise;
    part.nodesSweep     = info.nodesSweep;

    const std::string csvForm = form == MASS ? "mass" :
                                form == STIFFNESS ? "stiffness" : "biharmonic";
    writeCsvRow(cfg, csvForm, "gauss", nGauss, 0, 0.0, 0.0,
                tBuildGauss, tGauss.init, tGauss.assem, PartitionColumns());
    writeCsvRow(cfg, csvForm, "hpatch", nHPatch, hpatchRule.numFallbackCells(),
                errHPatch, relHPatch,
                tBuildHPatch, tHPatch.init, tHPatch.assem, part);

    // --- RHS diagnostic: b_i = ∫ N_i dΩ ---
    // Reveals per-function integration errors more directly than the matrix.
    // For the identity geometry, sum(b) = domain volume = 1 (partition of unity).
    // Truncated THB functions that are piecewise-polynomial within an atomic cell
    // would show non-zero relative error here if the rule misses their kink points.
    {
        const gsVector<T> bGauss  = assembleRHS<d,T>(spaceMB, spaceMB, mp, optGauss, nullptr);
        const gsVector<T> bHPatch = assembleRHS<d,T>(spaceMB, integMB, mp, optDefault, &hpatchFactory);

        const T rhsMag  = bGauss.cwiseAbs().maxCoeff();
        const T rhsErrH = (bHPatch - bGauss).cwiseAbs().maxCoeff();

        gsInfo << "  ∫N_i dΩ — sum(Gauss)=" << std::fixed << std::setprecision(15)
               << bGauss.sum()
               << "  sum(hpatch)=" << bHPatch.sum() << "\n";
        gsInfo << "  gsHPatchRule rel.err="
               << (rhsMag > 0 ? rhsErrH / rhsMag : rhsErrH) << "\n\n";
    }
}

// ============================================================
// Templated main
// ============================================================

template<short_t d, class T>
void run(BenchConfig & cfg)
{
    const index_t p = cfg.p, r = cfg.r;
    const bool plot = cfg.plot;

    // ------------------------------------------------------------------
    // Build THB-spline basis with configurable local refinement
    // ------------------------------------------------------------------
    gsKnotVector<T> kv(0, 1.0, 0, p + 1, 1);
    for (index_t i = 0; i < cfg.baseRefs; ++i)
        kv.uniformRefine();

    std::vector<gsKnotVector<T>> kvs;
    for (short_t dd = 0; dd < d; ++dd) kvs.push_back(kv);

    gsTensorBSplineBasis<d,T> tbsb(kvs);
    gsTHBSplineBasis<d,T>     thb(tbsb);

    gsInfo << "Building geometry..." << std::flush;
    gsStopwatch meshTimer;
    {
        thb.uniformRefine();
        for (index_t lvl = 1; lvl <= cfg.numLevels; ++lvl)
            refineLevel<d,T>(thb, cfg.pattern, cfg.baseRefs, cfg.band, lvl, cfg.admissible, cfg.admissibleJump);
    }
    cfg.tMesh     = meshTimer.stop();
    cfg.nElem     = thb.numElements();
    cfg.nDof      = thb.size();
    cfg.maxLevel  = thb.maxLevel();
    gsInfo << " done.\n";

    // The DoFs are ALWAYS defined by the original thb. The integration basis
    // (optionally refined via makeIntegrationBasis) only serves as the
    // integration domain for gsHPatchRule.
    gsTHBSplineBasis<d,T> integThb(thb);
    if (cfg.useIntegrationBasis)
    {
        integThb = gsHPatchRule<T>::template makeIntegrationBasis<d>(thb);
        gsInfo << "Integration basis: refined from " << thb.size()
               << " to " << integThb.size() << " DOF ("
               << integThb.numElements() << " elements)\n";
    }

    gsInfo << "Basis: degree=" << p << "  regularity=" << r
           << "  pattern=" << patternName(cfg.pattern)
           << "  admissible=" << (cfg.admissible ? "true" : "false")
           << "  jump=" << cfg.admissibleJump
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
    const bool all = cfg.form == "all";
    if (all || cfg.form == "mass")
        compareRules<d,T>(cfg, "Mass",      spaceMB, integMB, mp, MASS,      2*p, r);
    if (all || cfg.form == "stiffness")
        compareRules<d,T>(cfg, "Stiffness", spaceMB, integMB, mp, STIFFNESS, 2*p, r-1);

    if (cfg.form == "biharmonic" || (all && r >= 2))
        compareRules<d,T>(cfg, "Biharmonic", spaceMB, integMB, mp, BIHARMONIC, 2*p, r-2);
    else if (all)
        gsInfo << "Biharmonic: skipped (requires regularity >= 2; use -P 3 -R 2)\n";
}

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    index_t p          = 2;
    index_t r          = -1;
    index_t numLevels  = 3;
    index_t dim        = 2;
    index_t admissibleJump = 2;
    index_t reps       = 5;
    index_t baseRefs   = 2;
    bool    plot       = false;
    bool    integBasis = false;
    bool    admissible = false;
    std::string meshPattern = "lshape";
    std::string band = "scaled";
    std::string form = "all";
    std::string merge = "best";
    std::string csvFile;

    gsCmdLine cmd("Assembly benchmark: gsHPatchRule vs Gauss (gsExprAssembler).");
    cmd.addInt   ("P", "deg",    "Polynomial degree",           p);
    cmd.addInt   ("R", "reg",    "Regularity (-1: degree - 1)", r);
    cmd.addInt   ("l", "levels", "Number of refinement levels", numLevels);
    cmd.addInt   ("d", "dim",    "Problem dimension (1..4)",    dim);
    cmd.addInt   ("",  "admissibleJump",
                  "Jump parameter m for admissible refinement", admissibleJump);
    cmd.addString("",  "meshPattern",
                  "Mesh pattern: lshape, center, diagonal", meshPattern);
    cmd.addString("",  "form",
                  "Bilinear form: mass, stiffness, biharmonic, all", form);
    cmd.addString("",  "merge",
                  "Leaf merging of the leaf rule: none, pairwise, sweep, best", merge);
    cmd.addInt   ("",  "reps",
                  "Repetitions of assemble(); the minimum time is reported", reps);
    cmd.addString("",  "band",
                  "Refinement band variant: scaled (band width grows with --baseRefs) or thin (4 elements wide)", band);
    cmd.addInt   ("",  "baseRefs",
                  "Uniform refinements of the base knot vector "
                  "(2: 4 elements per direction before the level-0 refinement)", baseRefs);
    cmd.addString("",  "csv",
                  "Write one CSV row per (form, rule) to this file", csvFile);
    cmd.addSwitch("",  "admissible",
                  "Apply admissible closure to the marked refinement region", admissible);
    cmd.addSwitch("plot",        "Write ParaView mesh and the quadrature nodes of the best rule "
                                 "(whatever --merge says)", plot);
    cmd.addSwitch("integrationBasis",
                  "Refine THB via makeIntegrationBasis before assembly", integBasis);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // The timings compare single-thread assembly; they are serial by design.
    if (omp_get_max_threads() > 1)
    {
        gsWarn << "This benchmark requires OMP_NUM_THREADS=1 (found "
               << omp_get_max_threads() << " threads).\n";
        return 2;
    }

    if (r == -1)
        r = p - 1;
    if (r < 0 || r >= p)
    {
        gsInfo << "Error: regularity must be in [0, degree).\n";
        return EXIT_FAILURE;
    }
    if (admissibleJump < 1)
    {
        gsInfo << "Error: admissibleJump must be >= 1.\n";
        return EXIT_FAILURE;
    }
    if (reps < 1)
    {
        gsInfo << "Error: reps must be >= 1.\n";
        return EXIT_FAILURE;
    }
    if (baseRefs < 1)
    {
        gsInfo << "Error: baseRefs must be >= 1 (pattern boxes must align with parent elements).\n";
        return 2;
    }
    if (form != "mass" && form != "stiffness" && form != "biharmonic" && form != "all")
    {
        gsInfo << "Error: unknown form '" << form
               << "'. Expected one of: mass, stiffness, biharmonic, all.\n";
        return EXIT_FAILURE;
    }
    if (!isMergeMode(merge))
    {
        gsInfo << "Error: unknown merge '" << merge
               << "'. Expected one of: none, pairwise, sweep, best.\n";
        return EXIT_FAILURE;
    }
    if (form == "biharmonic" && r < 2)
    {
        gsInfo << "Error: the biharmonic form requires regularity >= 2.\n";
        return EXIT_FAILURE;
    }

    BandVariant bandVariant = BandScaled;
    if (!parseBandVariant(band, bandVariant))
    {
        gsInfo << "Error: unknown band '" << band
               << "'. Expected one of: scaled, thin.\n";
        return 2;
    }

    MeshPattern pattern = PatternLShape;
    if (!parseMeshPattern(meshPattern, pattern))
    {
        gsInfo << "Error: unknown meshPattern '" << meshPattern
               << "'. Expected one of: lshape, center, diagonal.\n";
        return EXIT_FAILURE;
    }

    std::ofstream csv;
    if (!csvFile.empty())
    {
        csv.open(csvFile.c_str(), std::ios::out | std::ios::trunc);
        if (!csv)
        {
            gsInfo << "Error: cannot open '" << csvFile << "' for writing.\n";
            return EXIT_FAILURE;
        }
        csv << csvHeader() << '\n';
    }

    BenchConfig cfg;
    cfg.dim = dim;
    cfg.p = p;
    cfg.r = r;
    cfg.numLevels = numLevels;
    cfg.pattern = pattern;
    cfg.admissible = admissible;
    cfg.admissibleJump = admissibleJump;
    cfg.plot = plot;
    cfg.useIntegrationBasis = integBasis;
    cfg.form = form;
    cfg.merge = merge;
    cfg.reps = reps;
    cfg.baseRefs = baseRefs;
    cfg.band = bandVariant;
    cfg.csv = csvFile.empty() ? NULL : &csv;
    cfg.nElem = cfg.nDof = cfg.maxLevel = 0;
    cfg.tMesh = 0;

    if      (dim == 1) run<1, real_t>(cfg);
    else if (dim == 2) run<2, real_t>(cfg);
    else if (dim == 3) run<3, real_t>(cfg);
    else if (dim == 4) run<4, real_t>(cfg);
    else { gsInfo << "Error: dimension must be 1..4.\n"; return EXIT_FAILURE; }

    return EXIT_SUCCESS;
}
