/** @file hierarchical_rule_cost.cpp

    @brief Rule-generation cost and node-count benchmark for gsHPatchRule.

    This driver deliberately performs NO assembly. It builds a hierarchical
    mesh, constructs the quadrature rules on it, and reports

      - the wall-clock cost of CONSTRUCTING the rule (the Newton solves), and
      - the number of quadrature nodes the rule produces,

    for a sweep of polynomial degrees. Because nothing is assembled, the
    numbers are unaffected by the per-thread rule construction inside the
    OpenMP region of gsExprAssembler (see gsExprAssembler.h, the comment above
    the thread-local QuRule): this driver is the place to sweep the degree far
    beyond the range at which a full assembly benchmark is affordable.

    Node counts are reported for three rules:

      Gauss   : element-wise Gauss-Legendre, (p+1)^d nodes per element.
      HPatch  : gsHPatchRule, the per-leaf rule of this paper.
      WQ*     : an ANALYTIC MODEL of hierarchical weighted quadrature, namely
                2 nodes per element per direction, i.e. 2^d per element. This
                is the interior count stated in section 4 of

                  C. Giannelli, T. Kanduc, M. Martinelli, G. Sangalli, M. Tani,
                  "Weighted quadrature for hierarchical B-splines",
                  CMAME 400 (2022) 115465,

                verbatim: "we have 2 quadrature points on each interior element
                of level n, or p+1 on the elements that touch the boundary".
                Both terms are evaluated per element, and the fraction of
                boundary-touching elements is reported next to the count.
                It is NOT an implementation of their method. Three caveats apply
                and must be repeated wherever this column is used:
                  (i)   their points live on the level-nu grid and coarse
                        functions inherit fine points, so the effective count
                        NEAR LEVEL TRANSITIONS exceeds 2^d per element;
                  (ii)  on coarse meshes and at high degree the boundary term
                        dominates (p+1 = 7 against 2 at p = 6), so the interior
                        asymptotic nElem * 2^d quoted in the paper is only
                        meaningful where bdryFrc is small;
                  (iii) their rule is WEIGHTED, so a node carries one weight per
                        test function, not one weight.

    Only truncated (THB) bases are swept. The rule depends on the hierarchical
    DOMAIN, not on the basis, so the node counts of an HB and a THB basis on
    the same mesh are identical by construction; HB-vs-THB is an EXACTNESS
    question, which is the subject of hierarchical_assembly_hpatch.cpp.

    Forms select the target space of Corollary 2.1 via the number of
    derivatives s: mass (s=0, regularity p-1), stiffness (s=1, p-2),
    biharmonic (s=2, p-3). A degree is skipped for a form when p-1-s < 0.

    Usage:
        ./hierarchical_rule_cost [--degrees 2,3,4] [-d dim] [-l levels]
                                 [--form mass|stiffness|biharmonic|all]
                                 [--meshPattern lshape|center|diagonal]
                                 [--admissible] [--admissibleJump m]
                                 [--baseRefs k] [--band scaled|thin] [--reps n] [--csv out.csv]

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

#include <fstream>
#include <limits>
#include <cstdlib>
#include <sstream>
#if !defined(_WIN32) && !defined(_WIN64)
#include <unistd.h>
#endif

using namespace gismo;

enum MeshPattern { PatternLShape, PatternCenter, PatternDiagonal };

// Width of the refinement bands: BandScaled keeps the refined region's physical
// size for every --baseRefs >= 1 (the centre box is shifted to an even origin,
// which moves it by one level-targetLevel element at --baseRefs 1); BandThin
// keeps the bands 4 target-level elements wide.
enum BandVariant { BandScaled, BandThin };

// ============================================================
// Mesh helpers (shared verbatim with hierarchical_assembly_hpatch.cpp)
// ============================================================

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
// Node counting and rule-construction timing
// ============================================================

/// Total number of nodes the rule selected in \a opts produces over the
/// domain of \a basis, together with the wall-clock time spent CONSTRUCTING
/// the rule (\a tBuild, seconds, best of \a reps) and the time spent mapping
/// it to every element (\a tMap, seconds, best of \a reps).
template<class T>
index_t countNodes(const gsBasis<T>    & basis,
                   const gsOptionList  & opts,
                   index_t               reps,
                   T                   & tBuild,
                   T                   & tMap)
{
    gsStopwatch clock;

    tBuild = std::numeric_limits<T>::max();
    for (index_t k = 0; k < reps; ++k)
    {
        clock.restart();
        typename gsQuadRule<T>::uPtr warm = gsQuadrature::getPtr(basis, opts);
        const T t = clock.stop();
        GISMO_ENSURE(warm != nullptr, "quadrature rule construction returned null");
        tBuild = math::min(tBuild, t);
    }

    typename gsQuadRule<T>::uPtr rule = gsQuadrature::getPtr(basis, opts);

    index_t n = 0;
    gsMatrix<T> nodes;
    gsVector<T> weights;
    tMap = std::numeric_limits<T>::max();
    for (index_t k = 0; k < reps; ++k)
    {
        index_t nPass = 0;
        clock.restart();
        for (auto & elem : basis.domain()->allElements())
        {
            rule->mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
            nPass += weights.size();
        }
        const T t = clock.stop();
        GISMO_ENSURE(k == 0 || nPass == n, "node count changed between mapping passes");
        n = nPass;
        tMap = math::min(tMap, t);
    }

    return n;
}

/// Run-wide metadata repeated in every CSV row.
struct RunInfo
{
    std::string pattern;   ///< mesh pattern name
    std::string band;      ///< band variant name
    index_t     threads;   ///< OpenMP thread count
    std::string host;      ///< host name, CSV-safe
    std::string cpu;       ///< CPU model, CSV-safe
};

/// Trims surrounding whitespace and replaces commas by semicolons so that the
/// string occupies exactly one CSV field.
std::string csvSafe(std::string str)
{
    for (size_t k = 0; k != str.size(); ++k)
        if (str[k] == ',') str[k] = ';';
    const std::string ws = " \t\r\n";
    const size_t first = str.find_first_not_of(ws);
    if (first == std::string::npos) return std::string();
    return str.substr(first, str.find_last_not_of(ws) - first + 1);
}

// ============================================================
// One (degree, form) measurement
// ============================================================

template<short_t d, class T>
void measure(index_t        p,
             index_t        s,
             const std::string & formName,
             index_t        numLevels,
             MeshPattern    pattern,
             BandVariant    bandVariant,
             bool           admissible,
             index_t        admissibleJump,
             index_t        reps,
             index_t        baseRefs,
             const RunInfo & info,
             std::ostream * csv)
{
    const index_t reg = p - 1 - s;   // target regularity at a simple knot
    if (reg < 0) return;             // form not representable at this degree

    // ---- mesh ----
    gsStopwatch clock;
    gsKnotVector<T> kv(0, 1.0, 0, p + 1, 1);
    for (index_t i = 0; i < baseRefs; ++i)
        kv.uniformRefine();

    std::vector<gsKnotVector<T>> kvs;
    for (short_t dd = 0; dd < d; ++dd) kvs.push_back(kv);

    gsTensorBSplineBasis<d,T> tbsb(kvs);
    gsTHBSplineBasis<d,T>     thb(tbsb);

    clock.restart();
    thb.uniformRefine();
    for (index_t lvl = 1; lvl <= numLevels; ++lvl)
        refineLevel<d,T>(thb, pattern, baseRefs, bandVariant, lvl, admissible, admissibleJump);
    const T tMesh = clock.stop();

    const index_t nElem = thb.numElements();
    GISMO_ENSURE(nElem > 0, "empty hierarchical mesh: nothing would be measured");

    // ---- rules ----
    gsOptionList optGauss = gsExprAssembler<T>::defaultOptions();
    optGauss.setInt ("quRule", gsQuadrature::GaussLegendre);
    optGauss.setReal("quA", static_cast<T>(1.0));
    optGauss.setInt ("quB", 1);

    gsOptionList optHPatch = gsExprAssembler<T>::defaultOptions();
    optHPatch.setInt ("quRule", gsQuadrature::HPatchRule);
    optHPatch.setReal("quA", static_cast<T>(2 * p));
    optHPatch.setInt ("quB", reg);

    T tBuildG = 0, tMapG = 0, tBuildH = 0, tMapH = 0;
    const index_t nGauss  = countNodes<T>(thb, optGauss,  reps, tBuildG, tMapG);
    const index_t nHPatch = countNodes<T>(thb, optHPatch, reps, tBuildH, tMapH);

    // Fallback diagnostics come from the concrete rule, not the base class.
    gsHPatchRule<T> hp(*thb.domain(), 2 * p, reg, false);
    const index_t nFallback = hp.numFallbackCells();

    // Node-count model of hierarchical weighted quadrature; see the file header
    // for what this model does and does not represent. Per direction an element
    // gets 2 nodes in the interior and p+1 if it is the first or last knot span
    // in that direction, so the boundary correction is NOT negligible on coarse
    // meshes or at high degree (p+1 = 7 against 2 at p = 6). We therefore also
    // report the fraction of elements touching the boundary: where that fraction
    // is large, the model is dominated by its boundary term and the interior
    // asymptotic nElem * 2^d that the paper quotes does not apply.
    index_t nWQmodel  = 0;
    index_t nBoundary = 0;
    {
        const gsMatrix<T> supp = thb.support();
        const T tol = 1e-12;
        for (auto & elem : thb.domain()->allElements())
        {
            const gsVector<T> lo = elem.lowerCorner();
            const gsVector<T> up = elem.upperCorner();
            index_t nodes  = 1;
            bool    onBdry = false;
            for (short_t k = 0; k != d; ++k)
            {
                const bool bdry = (lo[k] <= supp(k,0) + tol) || (up[k] >= supp(k,1) - tol);
                nodes *= bdry ? (p + 1) : 2;
                onBdry = onBdry || bdry;
            }
            nWQmodel += nodes;
            if (onBdry) ++nBoundary;
        }
    }
    const T bdryFrac = static_cast<T>(nBoundary) / static_cast<T>(nElem);

    const T ratioH  = static_cast<T>(nHPatch)  / static_cast<T>(nGauss);
    const T ratioWQ = static_cast<T>(nWQmodel) / static_cast<T>(nGauss);

    gsInfo << std::left
           << std::setw(12) << formName
           << std::setw(4)  << p
           << std::setw(9)  << nElem
           << std::setw(9)  << thb.size()
           << std::setw(11) << nGauss
           << std::setw(11) << nHPatch
           << std::setw(11) << nWQmodel
           << std::setw(9)  << std::fixed << std::setprecision(3) << ratioH
           << std::setw(9)  << std::fixed << std::setprecision(3) << ratioWQ
           << std::setw(12) << std::scientific << std::setprecision(3) << tBuildH
           << std::setw(12) << std::scientific << std::setprecision(3) << tMapH
           << std::setw(8)  << nFallback
           << std::setw(9)  << std::fixed << std::setprecision(3) << bdryFrac
           << "\n";

    if (csv != nullptr)
        (*csv) << formName << ',' << d << ',' << p << ',' << numLevels << ','
               << (admissible ? 1 : 0) << ',' << admissibleJump << ','
               << nElem << ',' << thb.size() << ',' << thb.maxLevel() << ','
               << nGauss << ',' << nHPatch << ',' << nWQmodel << ','
               << std::scientific << std::setprecision(9)
               << ratioH << ',' << ratioWQ << ','
               << tMesh << ',' << tBuildG << ',' << tMapG << ','
               << tBuildH << ',' << tMapH << ',' << nFallback << ','
               << nBoundary << ',' << bdryFrac << ','
               << (p - 1) << ',' << info.pattern << ',' << reps << ','
               << info.threads << ',' << info.host << ',' << info.cpu << ',' << baseRefs << ',' << info.band << '\n';
}

template<short_t d, class T>
void run(const std::vector<index_t> & degrees,
         const std::string & form,
         index_t numLevels,
         MeshPattern pattern,
         BandVariant bandVariant,
         bool admissible,
         index_t admissibleJump,
         index_t reps,
         index_t baseRefs,
         const RunInfo & info,
         std::ostream * csv)
{
    gsInfo << "dim=" << d << "  levels=" << numLevels
           << "  admissible=" << (admissible ? "true" : "false")
           << "  jump=" << admissibleJump
           << "  reps=" << reps
           << "  threads=" << info.threads << "\n"
           << "host=" << info.host << "  cpu=" << info.cpu << "\n\n";

    gsInfo << std::left
           << std::setw(12) << "form"   << std::setw(4)  << "p"
           << std::setw(9)  << "nElem"  << std::setw(9)  << "nDof"
           << std::setw(11) << "nGauss" << std::setw(11) << "nHPatch"
           << std::setw(11) << "nWQ*"   << std::setw(9)  << "H/G"
           << std::setw(9)  << "WQ*/G"  << std::setw(12) << "tBuild(s)"
           << std::setw(12) << "tMap(s)"<< std::setw(8)  << "nFallb"
           << std::setw(9)  << "bdryFrc"
           << "\n";

    for (size_t k = 0; k != degrees.size(); ++k)
    {
        const index_t p = degrees[k];
        if (form == "mass"       || form == "all")
            measure<d,T>(p, 0, "mass",       numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
        if (form == "stiffness"  || form == "all")
            measure<d,T>(p, 1, "stiffness",  numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
        if (form == "biharmonic" || form == "all")
            measure<d,T>(p, 2, "biharmonic", numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
    }
    gsInfo << "\n";
}

int main(int argc, char * argv[])
{
    index_t     dim            = 2;
    index_t     numLevels      = 3;
    index_t     admissibleJump = 2;
    index_t     reps           = 3;
    index_t     baseRefs       = 2;
    bool        admissible     = false;
    std::string degreeList     = "2,3,4";
    std::string form           = "all";
    std::string meshPattern    = "lshape";
    std::string band           = "scaled";
    std::string csvFile        = "";

    gsCmdLine cmd("Rule-generation cost and node counts for gsHPatchRule (no assembly).");
    cmd.addString("",  "degrees", "Comma-separated list of polynomial degrees", degreeList);
    cmd.addInt   ("d", "dim",     "Problem dimension (1..4)",                   dim);
    cmd.addInt   ("l", "levels",  "Number of refinement levels",                numLevels);
    cmd.addString("",  "form",    "mass, stiffness, biharmonic or all",         form);
    cmd.addString("",  "meshPattern", "Mesh pattern: lshape, center, diagonal", meshPattern);
    cmd.addSwitch("",  "admissible", "Apply admissible closure to the marked region", admissible);
    cmd.addInt   ("",  "admissibleJump", "Jump parameter m for admissible refinement", admissibleJump);
    cmd.addInt   ("",  "reps",    "Repetitions for the construction and mapping timings (best of)", reps);
    cmd.addString("",  "band",
                  "Refinement band variant: scaled (band width grows with --baseRefs) or thin (4 elements wide)", band);
    cmd.addInt   ("",  "baseRefs", "Uniform refinements of the base knot vector (2: 4 elements per direction before the level-0 refinement)", baseRefs);
    cmd.addString("",  "csv",     "Write one row per measurement to this file", csvFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (form != "mass" && form != "stiffness" && form != "biharmonic" && form != "all")
    {
        gsInfo << "Error: unknown form '" << form
               << "'. Expected one of: mass, stiffness, biharmonic, all.\n";
        return 1;
    }
    if (reps < 1)
    {
        gsInfo << "Error: reps must be >= 1.\n";
        return 1;
    }
    if (admissibleJump < 1)
    {
        gsInfo << "Error: admissibleJump must be >= 1.\n";
        return 1;
    }
    if (baseRefs < 1)
    {
        gsInfo << "Error: baseRefs must be >= 1 (pattern boxes must align with parent elements).\n";
        return 2;
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
        return 1;
    }

    std::vector<index_t> degrees;
    {
        std::stringstream ss(degreeList);
        std::string tok;
        while (std::getline(ss, tok, ','))
            if (!tok.empty()) degrees.push_back(std::atoi(tok.c_str()));
    }
    GISMO_ENSURE(!degrees.empty(),
                 "no degrees parsed from --degrees '" << degreeList
                 << "': the sweep would report nothing");
    for (size_t k = 0; k != degrees.size(); ++k)
        GISMO_ENSURE(degrees[k] >= 1, "degree must be >= 1, got " << degrees[k]);

    std::ofstream csvStream;
    std::ostream * csv = nullptr;
    if (!csvFile.empty())
    {
        csvStream.open(csvFile.c_str());
        GISMO_ENSURE(csvStream.is_open(), "cannot open csv file '" << csvFile << "'");
        csvStream << "form,dim,p,levels,admissible,jump,nElem,nDof,maxLevel,"
                     "nGauss,nHPatch,nWQmodel,ratioH,ratioWQ,"
                     "tMesh,tBuildGauss,tMapGauss,tBuildHPatch,tMapHPatch,nFallback,"
                     "nBoundaryElem,boundaryFraction,"
                     "r,pattern,reps,threads,host,cpu,baseRefs,band\n";
        csv = &csvStream;
    }

    RunInfo info;
    info.pattern = meshPattern;
    info.band    = band;
    info.threads = omp_get_max_threads();
    {
#if !defined(_WIN32) && !defined(_WIN64)
        char buf[256];
        const bool ok = (gethostname(buf, sizeof(buf)) == 0);
        buf[sizeof(buf) - 1] = '\0';
        info.host = csvSafe(ok ? std::string(buf) : std::string("unknown"));
#else
        const char * h = std::getenv("COMPUTERNAME");
        info.host = csvSafe(h ? std::string(h) : std::string("unknown"));
#endif
    }
    info.cpu = csvSafe(gsSysInfo::getCpuInfo());

    if      (dim == 1) run<1, real_t>(degrees, form, numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
    else if (dim == 2) run<2, real_t>(degrees, form, numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
    else if (dim == 3) run<3, real_t>(degrees, form, numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
    else if (dim == 4) run<4, real_t>(degrees, form, numLevels, pattern, bandVariant, admissible, admissibleJump, reps, baseRefs, info, csv);
    else { gsInfo << "Error: dimension must be 1..4.\n"; return 1; }

    if (csv != nullptr) csvStream.close();
    return 0;
}
