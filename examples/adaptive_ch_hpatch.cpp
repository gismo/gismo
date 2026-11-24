/** @file adaptive_ch_hpatch.cpp

    @brief Adaptive Cahn-Hilliard simulation on a THB-spline space, with the
    element-wise Gauss-Legendre rule of the module or gsHPatchRule (one reduced
    patch rule per kd leaf, or per merged box of leaves) as the rule of the solve.

    One run per call. The unit square carries a THB space of degree p on
    n0 x n0 base elements, refined in bands around the interface of the initial
    state (levels 1..L, see refineBands). The state is
        c(x) = max_i tanh( -(|x - x_i| - r_i) / (sqrt(2) eps) ),   eps = 2 h_L,
    with h_L the cell size of level L: c = +1 inside the bubbles and -1 outside.
    The model is the quartic potential in the AroundZero convention with
    constant mobility and lambda = eps^2, with the strong Clamped condition
    (dc/dn = 0) on all four sides. Time stepping is generalized-alpha with
    Newton, through gsCahnHilliardTimeStepper.

    After each solved step the mesh is adapted once, in place, on the one THB
    basis that serves as space and integration basis:
      - indicator   eta_K = 1 - |mean_K c|, evaluated with a Gauss rule;
      - refine      every element with eta_K > tauR and level l < L is
                    refined directly to level L;
      - coarsen     a parent cell is coarsened when all four children are
                    active elements of one level and each has eta < tauC.
    Because tauC < tauR the two sets are disjoint. The state and its time
    derivative are L2-projected onto the new space, the assemblers are
    re-initialised, the leaf quadrature factory is rebuilt (one cache serves the
    whole run) and the mass matrix is re-assembled.

    --rule gauss|leaf selects the rule of the solve. --shadow (only with
    --rule gauss) additionally assembles the residual and the Jacobian with a
    leaf-rule assembler on the same state after each step, and times both arms
    in alternation (warm-up call, then --reps repetitions, minima reported).

    The timings are serial by design, so the driver must be run with
    OMP_NUM_THREADS=1; otherwise it exits with code 2.
    Exit codes: 0 done, 1 invalid arguments or unopenable CSV file, 2 more than
    one thread, 3 Newton failure.
    It needs the gsCahnHilliard module (GISMO_OPTIONAL=gsCahnHilliard); without
    it the program only prints a message.

    Output: CSV (--csv file; the same lines are echoed to stdout). The file is
    a run header of lines "# key: value", the column line, one row per step
    (flushed after every row) and a last line "# status: completed" or
    "# status: newton-failed at step k". Reals are written in scientific format
    with 10 digits, integers plainly; an empty field means "not computed". Times
    are wall seconds. Columns:
      step, t (time at the end of the step), converged, newtonIt, resRatio,
      nElem, nDof, nBasis, maxLevel (the mesh the step was solved on),
      nAboveTauR (elements with eta > tauR before adapting), nRef (elements
      refined), nCrs (elements removed by coarsening), meshChanged, nElemNew,
      nUnresolved (pre-adapt elements with eta > tauR whose centre is not at
      level L after adapting),
      tWall, tStep, tRes, tJac, tFact, tSolve, tMass (stepper buckets; tMass
      includes the mass re-assembly after remeshing), tIndicator, tMark,
      tAdapt, tInit, tProj, tBuild (rebuild of the leaf factory for the next
      step), tDiag, tShadow,
      energy, mass, cmin, cmax, certified (on the solved mesh),
      ptsLeaf, numLeaves, numCells, nFallback, cacheHits, cacheSolves,
      cacheFailures (the leaf rule used by this step's solve or shadow arm; the
      cache counters are cumulative; empty in a pure Gauss run),
      ptsGauss, tResGauss, tJacGauss, tAsmGauss, tResLeaf, tJacLeaf, tAsmLeaf,
      R_step = tAsmLeaf/tAsmGauss, relDiffRes, relDiffJac (shadow runs only).
    tShadow lies inside tWall and outside every other bucket.

    Usage:
        OMP_NUM_THREADS=1 ./adaptive_ch_hpatch [-P deg] [--n0 n] [-l levels] [--beta b]
            [--ic bubbles3|bubbles2] [--dt dt] [-N steps] [--rhoinf r] [--tol t]
            [--maxit n] [--taur tr] [--tauc tc] [--rule gauss|leaf] [--shadow]
            [--reps n] [--csv file] [--caseId name] [--commit sha] [--modcommit sha]
    The bubbles are (0.34,0.5,r 0.15), (0.66,0.5,r 0.15) and, for bubbles3,
    (0.5,0.2,r 0.08); the two large ones merge.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Hugo Verhelst
*/

#include <gismo.h>

using namespace gismo;

#ifdef gsCahnHilliard_ENABLED

#include <gsAssembler/gsHPatchRule.h>
#include <gsCahnHilliard/gsCahnHilliardAssembler.h>
#include <gsCahnHilliard/gsCahnHilliardTimeStepper.h>
#include <gsDomain/gsHDomainIterator.h>

typedef gsCahnHilliardAssembler<real_t, CH::Convention::AroundZero,
                                CH::Mobility::Constant, CH::Potential::Quartic> CHAssembler;
typedef gsCahnHilliardTimeStepper<real_t, CH::Convention::AroundZero,
                                  CH::Mobility::Constant, CH::Potential::Quartic> CHStepper;

// ============================================================
// Interface geometry and initial state
// ============================================================

/// Union of disjoint circles: signed distance (negative inside) and the exact
/// minimum of its modulus over an axis-aligned box.
class Interface
{
public:
    Interface(const std::vector<real_t> & cx, const std::vector<real_t> & cy,
              const std::vector<real_t> & r)
    : m_cx(cx), m_cy(cy), m_r(r) { }

    /// Signed distance at (x,y); for several circles the smallest in modulus.
    real_t distance(real_t x, real_t y) const
    {
        real_t best = std::numeric_limits<real_t>::max();
        for (size_t c = 0; c < m_r.size(); ++c)
        {
            const real_t d = std::hypot(x - m_cx[c], y - m_cy[c]) - m_r[c];
            if (math::abs(d) < math::abs(best))
                best = d;
        }
        return best;
    }

    /// min over the box [x0,x1]x[y0,y1] of |d|, computed exactly (no sampling).
    /// Exact because the circles are disjoint, so the band {|d| < delta} is the
    /// union of the bands of the single circles.
    real_t minAbsDistance(real_t x0, real_t x1, real_t y0, real_t y1) const
    {
        real_t best = std::numeric_limits<real_t>::max();
        for (size_t c = 0; c < m_r.size(); ++c)
        {
            // dn / df: smallest / largest distance from the centre to the box.
            const real_t dxn = std::max(std::max(x0 - m_cx[c], m_cx[c] - x1), real_t(0));
            const real_t dyn = std::max(std::max(y0 - m_cy[c], m_cy[c] - y1), real_t(0));
            const real_t dn  = std::hypot(dxn, dyn);
            const real_t dxf = std::max(math::abs(x0 - m_cx[c]), math::abs(x1 - m_cx[c]));
            const real_t dyf = std::max(math::abs(y0 - m_cy[c]), math::abs(y1 - m_cy[c]));
            const real_t df  = std::hypot(dxf, dyf);
            real_t m = 0;
            if (m_r[c] < dn)
                m = dn - m_r[c];
            else if (m_r[c] > df)
                m = m_r[c] - df;
            best = std::min(best, m);
        }
        return best;
    }

private:
    std::vector<real_t> m_cx, m_cy, m_r;
};

/// The diffuse-interface state c = max_i tanh(-(|x - x_i| - r_i) / (sqrt(2) eps)):
/// +1 inside the bubbles, -1 outside. For the quartic potential (c^2-1)^2/4 and
/// lambda = eps^2 the tanh is the 1D equilibrium profile across a flat interface.
class BubbleState : public gsFunction<real_t>
{
public:
    BubbleState(const Interface & interface, real_t eps)
    : m_interface(interface), m_eps(eps) { }

    short_t domainDim() const { return 2; }
    short_t targetDim() const { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const
    {
        result.resize(1, u.cols());
        for (index_t i = 0; i < u.cols(); ++i)
            result(0,i) = std::tanh(-m_interface.distance(u(0,i), u(1,i)) / (math::sqrt(real_t(2)) * m_eps));
    }

private:
    const Interface & m_interface;
    real_t m_eps;
};

/// Refines \a thb in bands around the interface: for l = 0, ..., L-1 every
/// level-l cell with min |d| < beta h_l becomes one box of level l+1.
/// A level-l cell meeting the band of width beta h_l lies in its level-(l-1)
/// parent, which meets the wider band of level l-1 and has been refined, so
/// every marked cell is active at level l. Cost: O((n0 2^l)^2) per level.
void refineBands(gsTHBSplineBasis<2,real_t> & thb, const Interface & interface,
                 index_t n0, index_t numLevels, real_t beta)
{
    for (index_t l = 0; l < numLevels; ++l)
    {
        const index_t nl = n0 << l;
        const real_t  h  = real_t(1) / nl;
        const real_t  delta = beta * h;
        std::vector<index_t> boxes;
        for (index_t i = 0; i < nl; ++i)
            for (index_t j = 0; j < nl; ++j)
                if (interface.minAbsDistance(i*h, (i+1)*h, j*h, (j+1)*h) < delta)
                {
                    // box format: level, lower corner, upper corner, in the
                    // index space of that level
                    boxes.push_back(l + 1);
                    boxes.push_back(2*i);     boxes.push_back(2*j);
                    boxes.push_back(2*i + 2); boxes.push_back(2*j + 2);
                }
        if (!boxes.empty())
            thb.refineElements(boxes);
    }
}

// ============================================================
// Mesh data and marking
// ============================================================

/// An active element: its level, the index (i,j) of its cell in the element
/// index space of that level, the indicator eta and the cell centre.
struct ElemData
{
    index_t level, i, j;
    real_t  eta, cx, cy;
};

/// Maximum element level of \a basis.
index_t maxElementLevel(const gsBasis<real_t> & basis)
{
    index_t m = 0;
    const auto dom = basis.domain();
    for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
        m = std::max(m, static_cast<index_t>(static_cast<const gsHDomainIterator<real_t,2> *>(it.get())->getLevel()));
    return m;
}

/// Collects level, cell index, indicator and centre of every active element of
/// the hierarchical \a basis, indexed by the element id. The knots are uniform
/// and dyadic, so the cell index at level l is the lower corner times n0 2^l.
std::vector<ElemData> gatherElements(const gsBasis<real_t> & basis, index_t n0,
                                     const std::vector<real_t> & eta, index_t maxLev)
{
    std::vector<ElemData> elems(basis.numElements());
    const auto dom = basis.domain();
    for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
    {
        ElemData d;
        d.level = static_cast<const gsHDomainIterator<real_t,2> *>(it.get())->getLevel();
        GISMO_ENSURE(d.level <= maxLev, "Element level " << d.level << " exceeds L = " << maxLev << ".");
        const gsVector<real_t> lo = it.lowerCorner();
        const gsVector<real_t> up = it.upperCorner();
        const real_t cellsPerDir = real_t(n0) * real_t(index_t(1) << d.level);
        d.i   = static_cast<index_t>(std::lround(lo(0) * cellsPerDir));
        d.j   = static_cast<index_t>(std::lround(lo(1) * cellsPerDir));
        d.eta = eta[it.id()];
        d.cx  = real_t(0.5) * (lo(0) + up(0));
        d.cy  = real_t(0.5) * (lo(1) + up(1));
        elems[it.id()] = d;
    }
    return elems;
}

// ============================================================
// Helpers
// ============================================================

/// Total number of quadrature points that assemble() evaluates when it maps
/// \a rule onto every element of the domain of \a basis.
index_t countPoints(const gsBasis<real_t> & basis, const gsQuadRule<real_t> & rule)
{
    index_t n = 0;
    gsMatrix<real_t> nodes; gsVector<real_t> weights;
    for (auto & elem : basis.domain()->allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        n += weights.size();
    }
    return n;
}

/// The two timed calls. Both take the state as a gsMultiPatch, so that no
/// constructSolution() is part of the timing.
void callResidual(CHAssembler & assembler, const gsMultiPatch<real_t> & state)
{ assembler.assembleResidual(state); }

void callJacobian(CHAssembler & assembler, const gsMultiPatch<real_t> & state)
{ assembler.assembleJacobian(state); }

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

/// Minimum over the repetitions of the residual time, the Jacobian time and
/// their sum within one repetition.
struct ArmTimes
{
    ArmTimes()
    : tRes(std::numeric_limits<double>::max()), tJac(std::numeric_limits<double>::max()),
      tAsm(std::numeric_limits<double>::max()) { }

    void update(double res, double jac)
    {
        tRes = std::min(tRes, res);
        tJac = std::min(tJac, jac);
        tAsm = std::min(tAsm, res + jac);
    }

    double tRes, tJac, tAsm;
};

/// Statistics of one leaf rule (one build of the factory) and the cumulative
/// counters of the cache at the time of that build.
struct LeafStats
{
    LeafStats() : ptsLeaf(0), numLeaves(0), numCells(0), nFallback(0),
                  cacheHits(0), cacheSolves(0), cacheFailures(0) { }
    index_t ptsLeaf, numLeaves, numCells, nFallback, cacheHits, cacheSolves, cacheFailures;
};

/// Statistics of the rule that \a factory serves for patch 0 of \a basis. The
/// handle shares the data of the rule the assembler uses.
LeafStats leafStats(const gsExprAssembler<real_t>::QuadratureFactory & factory,
                    const gsBasis<real_t> & basis, const gsOptionList & options,
                    const gsHPatchRuleCache<real_t> & cache)
{
    const gsQuadRule<real_t>::uPtr handle = factory(basis, options, 0, -1);
    const gsHPatchRule<real_t> & rule = dynamic_cast<const gsHPatchRule<real_t> &>(*handle);
    LeafStats s;
    s.numLeaves     = rule.partitionInfo().numLeaves;
    s.numCells      = rule.partitionInfo().numCells;
    s.nFallback     = rule.numFallbackCells();
    s.ptsLeaf       = countPoints(basis, rule);
    s.cacheHits     = cache.numHits();
    s.cacheSolves   = cache.numSolves();
    s.cacheFailures = cache.numFailures();
    GISMO_ENSURE(s.ptsLeaf == rule.numPointsTotal(),
                 "The element-wise point count (" << s.ptsLeaf
                 << ") differs from the point count of the rule ("
                 << rule.numPointsTotal() << ").");
    return s;
}

/// Splits a comma-separated list.
std::vector<std::string> splitList(const std::string & s)
{
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string item;
    while (std::getline(in, item, ','))
        out.push_back(item);
    return out;
}

const char * csvColumns()
{
    return "step,t,converged,newtonIt,resRatio,nElem,nDof,nBasis,maxLevel,nAboveTauR,nRef,nCrs,"
           "meshChanged,nElemNew,nUnresolved,tWall,tStep,tRes,tJac,tFact,tSolve,tMass,tIndicator,"
           "tMark,tAdapt,tInit,tProj,tBuild,tDiag,tShadow,energy,mass,cmin,cmax,certified,ptsLeaf,"
           "numLeaves,numCells,nFallback,cacheHits,cacheSolves,cacheFailures,ptsGauss,tResGauss,"
           "tJacGauss,tAsmGauss,tResLeaf,tJacLeaf,tAsmLeaf,R_step,relDiffRes,relDiffJac";
}

/// One CSV row: fields are set by column name; an unset field is written empty.
class CsvRow
{
public:
    CsvRow() : m_cols(splitList(csvColumns())), m_vals(m_cols.size()) { }

    template<class V>
    void set(const std::string & name, V value)
    {
        std::ostringstream os;
        write(os, value, std::is_floating_point<V>());
        m_vals[indexOf(name)] = os.str();
    }

    std::string str() const
    {
        std::string s;
        for (size_t i = 0; i < m_vals.size(); ++i)
            s += (i ? "," : "") + m_vals[i];
        return s;
    }

private:
    template<class V> static void write(std::ostream & os, V v, std::true_type)
    { os << std::scientific << std::setprecision(10) << v; }
    template<class V> static void write(std::ostream & os, V v, std::false_type)
    { os << v; }

    size_t indexOf(const std::string & name) const
    {
        for (size_t i = 0; i < m_cols.size(); ++i)
            if (m_cols[i] == name)
                return i;
        GISMO_ERROR("Unknown CSV column '" << name << "'.");
    }

    std::vector<std::string> m_cols;
    std::vector<std::string> m_vals;
};

/// Writes a line to stdout and, when open, to the CSV file (flushed).
struct Output
{
    void line(const std::string & s)
    {
        gsInfo << s << std::endl;
        if (csv.is_open())
            csv << s << std::endl;
    }
    std::ofstream csv;
};

template<class V>
std::string headerLine(const std::string & key, const V & value)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(10) << "# " << key << ": ";
    if (std::is_floating_point<V>::value)
        os << value;
    else
    {
        os.unsetf(std::ios::floatfield);
        os << value;
    }
    return os.str();
}

// ============================================================
// Main
// ============================================================

int main(int argc, char * argv[])
{
    index_t p = 3;
    index_t n0 = 16;
    index_t numLevels = 2;
    real_t  beta = 2;
    std::string ic = "bubbles3";
    real_t  dt = 1e-5;
    index_t numSteps = 10;
    real_t  rhoInf = 0.5;
    real_t  tol = 1e-4;
    index_t maxIt = 50;
    real_t  tauR = 0.05;
    real_t  tauC = 0.005;
    std::string rule = "gauss";
    bool    shadow = false;
    index_t reps = 5;
    std::string csvFile;
    std::string caseId = "adhoc";
    std::string commit = "unknown";
    std::string modCommit = "unknown";

    gsCmdLine cmd("Adaptive Cahn-Hilliard simulation on a THB space: Gauss or "
                  "gsHPatchRule as the rule of the solve, with per-step timing buckets (CSV).");
    cmd.addInt   ("P", "deg",      "Polynomial degree of the space basis (>= 2)", p);
    cmd.addInt   ("",  "n0",       "Elements per direction of the base mesh", n0);
    cmd.addInt   ("l", "levels",   "Maximum refinement level L (>= 1)", numLevels);
    cmd.addReal  ("",  "beta",     "Half-width of the initial refined bands in cells of the refined level", beta);
    cmd.addString("",  "ic",       "Initial state: bubbles3 or bubbles2", ic);
    cmd.addReal  ("",  "dt",       "Time step", dt);
    cmd.addInt   ("N", "steps",    "Number of time steps", numSteps);
    cmd.addReal  ("",  "rhoinf",   "RhoInf of the generalized-alpha scheme", rhoInf);
    cmd.addReal  ("",  "tol",      "Relative Newton tolerance", tol);
    cmd.addInt   ("",  "maxit",    "Maximum number of Newton iterations (>= 2)", maxIt);
    cmd.addReal  ("",  "taur",     "Refine elements with indicator eta > taur", tauR);
    cmd.addReal  ("",  "tauc",     "Coarsen parents whose four children have eta < tauc (0 < tauc < taur)", tauC);
    cmd.addString("",  "rule",     "Quadrature rule of the solve: gauss or leaf", rule);
    cmd.addSwitch("",  "shadow",   "Also time a leaf-rule assembler against Gauss after each step (needs --rule gauss)", shadow);
    cmd.addInt   ("",  "reps",     "Repetitions of the shadow timing; the minimum is reported", reps);
    cmd.addString("",  "csv",      "CSV output path (run header, one row per step, status line); empty: stdout only", csvFile);
    cmd.addString("",  "caseId",   "Label of this run in the CSV header", caseId);
    cmd.addString("",  "commit",   "Root commit, echoed into the header", commit);
    cmd.addString("",  "modcommit","Module commit, echoed into the header", modCommit);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // The timings are serial by design.
    if (omp_get_max_threads() > 1)
    {
        gsWarn << "This benchmark requires OMP_NUM_THREADS=1 (found "
               << omp_get_max_threads() << " threads).\n";
        return 2;
    }

    if (p < 2)            { gsInfo << "Error: the degree must be >= 2.\n"; return EXIT_FAILURE; }
    if (n0 < 1)           { gsInfo << "Error: n0 must be >= 1.\n"; return EXIT_FAILURE; }
    if (numLevels < 1)    { gsInfo << "Error: levels must be >= 1.\n"; return EXIT_FAILURE; }
    if (beta <= 0)        { gsInfo << "Error: beta must be > 0.\n"; return EXIT_FAILURE; }
    if (ic != "bubbles3" && ic != "bubbles2")
    {
        gsInfo << "Error: unknown ic '" << ic << "'. Expected bubbles3 or bubbles2.\n";
        return EXIT_FAILURE;
    }
    if (dt <= 0)          { gsInfo << "Error: dt must be > 0.\n"; return EXIT_FAILURE; }
    if (numSteps < 1)     { gsInfo << "Error: steps must be >= 1.\n"; return EXIT_FAILURE; }
    if (tol <= 0)         { gsInfo << "Error: tol must be > 0.\n"; return EXIT_FAILURE; }
    if (maxIt < 2)        { gsInfo << "Error: maxit must be >= 2.\n"; return EXIT_FAILURE; }
    if (!(0 < tauC && tauC < tauR && tauR < 1))
    {
        gsInfo << "Error: the thresholds must satisfy 0 < tauc < taur < 1.\n";
        return EXIT_FAILURE;
    }
    if (rule != "gauss" && rule != "leaf")
    {
        gsInfo << "Error: unknown rule '" << rule << "'. Expected gauss or leaf.\n";
        return EXIT_FAILURE;
    }
    if (shadow && rule != "gauss")
    {
        gsInfo << "Error: --shadow needs --rule gauss.\n";
        return EXIT_FAILURE;
    }
    if (reps < 1)         { gsInfo << "Error: reps must be >= 1.\n"; return EXIT_FAILURE; }

    Output out;
    if (!csvFile.empty())
    {
        out.csv.open(csvFile.c_str(), std::ios::out | std::ios::trunc);
        if (!out.csv)
        {
            gsInfo << "Error: cannot open '" << csvFile << "' for writing.\n";
            return EXIT_FAILURE;
        }
    }

    // The leaf rule is built for the solve (--rule leaf) or for the shadow arm.
    const bool leafRun = (rule == "leaf");
    const bool useLeaf = leafRun || shadow;

    // ---- Interface and mesh -------------------------------------------------
    std::vector<real_t> cx = {0.34, 0.66}, cy = {0.5, 0.5}, cr = {0.15, 0.15};
    if (ic == "bubbles3")
    {
        cx.push_back(0.5); cy.push_back(0.2); cr.push_back(0.08);
    }
    const Interface interface(cx, cy, cr);

    // Interface width: twice the cell size of the finest level.
    const real_t hL  = real_t(1) / (real_t(n0) * std::pow(real_t(2), real_t(numLevels)));
    const real_t eps = 2 * hL;
    const real_t lambda = eps * eps;
    const index_t gaussPtsPerElem = (2*p + 1) * (2*p + 1);

    // Geometry: the bilinear unit square, one element (affine map).
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare());
    mp.addAutoBoundaries();

    // Space basis: degree p, maximal smoothness, n0 x n0 elements, then the
    // initial band refinement.
    gsKnotVector<real_t> kv(0, 1, n0 - 1, p + 1, 1);
    gsTensorBSplineBasis<2,real_t> tbsb(kv, kv);
    gsTHBSplineBasis<2,real_t> thb(tbsb);
    refineBands(thb, interface, n0, numLevels, beta);

    // The one multi-basis of the run. Every assembler borrows it as space and
    // integration basis, so it is declared before them and never reassigned;
    // mesh adaptation refines and coarsens mb.basis(0) in place.
    gsMultiBasis<real_t> mb;
    mb.addBasis(thb.clone().release());
    mb.setTopology(mp);

    // Strong Clamped condition (dc/dn = 0) on all four sides.
    gsBoundaryConditions<real_t> bc;
    for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::clamped, 0);
    bc.setGeoMap(mp);

    // ---- Assemblers ---------------------------------------------------------
    // The cache of the leaf rule outlives every factory built in the run and
    // is not thread-safe (the run is serial).
    gsHPatchRuleCache<real_t> cache;

    CHAssembler asmMain(mp, mb, bc);
    memory::unique_ptr<CHAssembler> asmShadow;
    if (shadow)
        asmShadow.reset(new CHAssembler(mp, mb, bc));

    std::vector<CHAssembler *> assemblers(1, &asmMain);
    if (shadow)
        assemblers.push_back(asmShadow.get());
    for (CHAssembler * a : assemblers)
    {
        a->options().setReal  ("Lambda", lambda);
        a->options().setSwitch("AssembleWeakBCs", false);
        // The space basis is also the integration basis.
        a->options().setSwitch("ExprAssembler.SameElement", true);
        // (2p+1)^2 Gauss points per element: exact for the residual and
        // Jacobian integrands of the quartic model on an affine geometry.
        a->options().setReal  ("ExprAssembler.quA", 2.0);
        a->options().setInt   ("ExprAssembler.quB", 1);
        a->initialize();
    }

    const index_t nElem0 = mb.basis(0).numElements();
    {
        const gsQuadRule<real_t>::uPtr gaussRule =
            gsQuadrature::getPtr(mb.basis(0), asmMain.options().getGroup("ExprAssembler"));
        GISMO_ENSURE(gaussRule->numNodes() == gaussPtsPerElem,
                     "The Gauss rule has " << gaussRule->numNodes() << " nodes per element, not (2p+1)^2 = "
                     << gaussPtsPerElem << ".");
        GISMO_ENSURE(countPoints(mb.basis(0), *gaussRule) == nElem0 * gaussPtsPerElem,
                     "The Gauss point count differs from nElem*(2p+1)^2.");
    }

    // ---- Leaf rule ------------------------------------------------------------
    // One rule for every integral of the residual and Jacobian: degree 4p,
    // regularity p-3 (-1 at p=2 means no continuity assumed). The factory
    // clones mb, so it is rebuilt after every mesh change.
    CHAssembler * leafAsm = shadow ? asmShadow.get() : &asmMain;
    gsStopwatch sw;
    LeafStats leaf;
    double setupBuild = 0;
    auto buildLeaf = [&]() -> double
    {
        sw.restart();
        gsExprAssembler<real_t>::QuadratureFactory factory =
            gsHPatchRule<real_t>::makeFactory(mb, 4*p, p - 3, false, &cache,
                                              gsHPatchRule<real_t>::MergeMode::best);
        leafAsm->setQuadratureFactory(factory);
        const double t = sw.stop();
        leaf = leafStats(factory, mb.basis(0), leafAsm->options().getGroup("ExprAssembler"), cache);
        return t;
    };
    if (useLeaf)
        setupBuild = buildLeaf();

    // ---- Stepper -----------------------------------------------------------------
    // The stepper borrows asmMain, so the factory installed on it reaches the
    // mass matrix and every Newton assembly.
    CHStepper stepper(asmMain);
    stepper.options().setReal("TimeStep", dt);
    stepper.options().setReal("RhoInf", rhoInf);
    stepper.options().setInt ("Method", 1);
    stepper.options().setInt ("MaxIterations", maxIt);
    stepper.options().setReal("ToleranceRel", tol);
    stepper.initialize();
    const double setupMass = stepper.timings().mass;
    stepper.resetTimings();

    // ---- Initial state --------------------------------------------------------------
    // Quasi-interpolant of the bubble profile, brought into the Clamped space
    // by keeping only the free DoFs; dc/dt = 0 at t = 0.
    const BubbleState fun(interface, eps);
    gsMatrix<real_t> coefs;
    gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0), fun, coefs);
    gsMultiPatch<real_t> C0;
    C0.addPatch(mb.basis(0).makeGeometry(give(coefs)));
    gsMatrix<real_t> Cold, dCold, Cnew, dCnew;
    asmMain.constructSolution(C0, Cold);
    dCold.setZero(asmMain.numDofs(), 1);
    gsMultiPatch<real_t> Cmp, dCmp;
    asmMain.constructSolution(Cold, Cmp);
    asmMain.constructSolution(dCold, dCmp);

    // ---- Header -------------------------------------------------------------------------
#ifdef GISMO_WITH_PARDISO
    const std::string solverName = "PardisoLU";
#else
    const std::string solverName = "SparseLU";
#endif
    out.line(headerLine("driver", std::string("adaptive_ch_hpatch")));
    out.line(headerLine("caseId", caseId));
    out.line(headerLine("rule", rule));
    out.line(headerLine("shadow", index_t(shadow)));
    out.line(headerLine("p", p));
    out.line(headerLine("n0", n0));
    out.line(headerLine("L", numLevels));
    out.line(headerLine("beta", beta));
    out.line(headerLine("ic", ic));
    out.line(headerLine("eps", eps));
    out.line(headerLine("lambda", lambda));
    out.line(headerLine("dt", dt));
    out.line(headerLine("N", numSteps));
    out.line(headerLine("rhoInf", rhoInf));
    out.line(headerLine("tol", tol));
    out.line(headerLine("maxIt", maxIt));
    out.line(headerLine("tauR", tauR));
    out.line(headerLine("tauC", tauC));
    out.line(headerLine("reps", reps));
    out.line(headerLine("quA", real_t(2)));
    out.line(headerLine("quB", index_t(1)));
    out.line(headerLine("gaussPtsPerElem", gaussPtsPerElem));
    out.line(headerLine("leafDegree", 4*p));
    out.line(headerLine("leafRegularity", p - 3));
    out.line(headerLine("merge", std::string("best")));
    out.line(headerLine("solver", solverName));
    out.line(headerLine("threads", index_t(omp_get_max_threads())));
    out.line(headerLine("host", hostName()));
    out.line(headerLine("cpu", cpuModel()));
    out.line(headerLine("commit", commit));
    out.line(headerLine("modcommit", modCommit));
    out.line(headerLine("setup_nElem", nElem0));
    out.line(headerLine("setup_nDof", asmMain.numDofs()));
    out.line(headerLine("setup_tBuild", setupBuild));
    out.line(headerLine("setup_tMass", setupMass));
    out.line(csvColumns());

    // ---- Time loop --------------------------------------------------------------------------
    // Per step: solve, diagnostics, optional shadow timing, indicator,
    // marking, adaptation, and (when the mesh changed) re-initialisation,
    // transfer, factory rebuild and mass re-assembly. The stepper's timers
    // accumulate, so the buckets of a step are differences of its totals.
    const std::vector<std::string> timeCols = {"tIndicator", "tMark", "tAdapt", "tInit",
                                               "tProj", "tBuild", "tDiag", "tShadow"};
    gsStopwatch swWall;
    real_t t = 0;
    for (index_t step = 1; step <= numSteps; ++step)
    {
        swWall.restart();
        const gsCahnHilliardStepTimings<real_t> T0 = stepper.timings();
        const index_t nElem  = mb.basis(0).numElements();
        const index_t nDof   = asmMain.numDofs();
        const index_t nBasis = mb.basis(0).size();
        GISMO_ENSURE(Cold.rows() == nDof, "The state has " << Cold.rows() << " rows, not nDof = " << nDof << ".");

        // The rule of this step's solve; the rebuild after adapting overwrites
        // `leaf` for the next step.
        const LeafStats leafUsed = leaf;
        CsvRow row;
        auto setLeafColumns = [&]()
        {
            row.set("ptsLeaf", leafUsed.ptsLeaf);       row.set("numLeaves", leafUsed.numLeaves);
            row.set("numCells", leafUsed.numCells);     row.set("nFallback", leafUsed.nFallback);
            row.set("cacheHits", leafUsed.cacheHits);   row.set("cacheSolves", leafUsed.cacheSolves);
            row.set("cacheFailures", leafUsed.cacheFailures);
        };
        auto setStepperColumns = [&](const gsCahnHilliardStepTimings<real_t> & T1)
        {
            row.set("tRes",   T1.residual  - T0.residual);
            row.set("tJac",   T1.jacobian  - T0.jacobian);
            row.set("tFact",  T1.factorize - T0.factorize);
            row.set("tSolve", T1.solve     - T0.solve);
            row.set("tMass",  T1.mass      - T0.mass);
            row.set("newtonIt", T1.newtonIterations - T0.newtonIterations);
        };

        // 1. Solve. Cold/Cnew must be distinct objects.
        sw.restart();
        const gsCahnHilliardStepStatus<real_t> status = stepper.step(Cold, dCold, Cnew, dCnew);
        const double tStep = sw.stop();
        t += dt;
        const gsCahnHilliardStepTimings<real_t> Tafter = stepper.timings();

        if (!status.converged)
        {
            row.set("step", step);
            row.set("t", t);
            row.set("converged", index_t(0));
            row.set("resRatio", status.residualRatio);
            row.set("nElem", nElem);
            row.set("nDof", nDof);
            row.set("nBasis", nBasis);
            row.set("maxLevel", maxElementLevel(mb.basis(0)));
            row.set("tStep", tStep);
            row.set("tWall", swWall.stop());
            setStepperColumns(Tafter);
            if (useLeaf)
                setLeafColumns();
            out.line(row.str());
            out.line(headerLine("status", "newton-failed at step " + std::to_string(step)));
            return 3;
        }
        Cold.swap(Cnew);
        dCold.swap(dCnew);

        for (const std::string & c : timeCols)
            row.set(c, 0.0);

        // 2. The end-of-step state as a function on the solved mesh.
        double tDiag = 0;
        sw.restart();
        asmMain.constructSolution(Cold, Cmp);
        asmMain.constructSolution(dCold, dCmp);
        tDiag += sw.stop();

        // 3. Shadow timing: Gauss arm = asmMain, leaf arm = asmShadow, same
        //    state, same protocol as the stage-A benchmark.
        if (shadow)
        {
            sw.restart();
            const gsQuadRule<real_t>::uPtr gaussRule =
                gsQuadrature::getPtr(mb.basis(0), asmMain.options().getGroup("ExprAssembler"));
            const index_t ptsGauss = countPoints(mb.basis(0), *gaussRule);
            GISMO_ENSURE(ptsGauss == nElem * gaussPtsPerElem,
                         "The Gauss point count (" << ptsGauss << ") differs from nElem*(2p+1)^2 = "
                         << nElem * gaussPtsPerElem << ".");

            callResidual(asmMain, Cmp);      callJacobian(asmMain, Cmp);
            callResidual(*asmShadow, Cmp);   callJacobian(*asmShadow, Cmp);
            ArmTimes tG, tL;
            gsStopwatch swArm;
            double res, jac;
            for (index_t rep = 0; rep < reps; ++rep)
            {
                swArm.restart(); callResidual(asmMain, Cmp);    res = swArm.stop();
                swArm.restart(); callJacobian(asmMain, Cmp);    jac = swArm.stop();
                tG.update(res, jac);

                swArm.restart(); callResidual(*asmShadow, Cmp); res = swArm.stop();
                swArm.restart(); callJacobian(*asmShadow, Cmp); jac = swArm.stop();
                tL.update(res, jac);
            }
            // Each arm holds the output of its last assembly.
            const gsMatrix<real_t> rG = asmMain.rhs();
            const gsMatrix<real_t> rL = asmShadow->rhs();
            const gsSparseMatrix<real_t> JG = asmMain.matrix();
            const gsSparseMatrix<real_t> dJ = asmShadow->matrix() - JG;
            row.set("ptsGauss", ptsGauss);
            row.set("tResGauss", tG.tRes); row.set("tJacGauss", tG.tJac); row.set("tAsmGauss", tG.tAsm);
            row.set("tResLeaf",  tL.tRes); row.set("tJacLeaf",  tL.tJac); row.set("tAsmLeaf",  tL.tAsm);
            row.set("R_step", tL.tAsm / tG.tAsm);
            row.set("relDiffRes", real_t((rL - rG).norm() / rG.norm()));
            row.set("relDiffJac", real_t(dJ.norm() / JG.norm()));
            row.set("tShadow", sw.stop());
        }

        // 4. Diagnostics on the converged state.
        sw.restart();
        const real_t energy = asmMain.computeEnergy(Cmp);
        const real_t mass   = asmMain.computeMass(Cmp);
        const gsCahnHilliardBounds<real_t> bounds = asmMain.checkConcentrationBounds(Cmp);
        tDiag += sw.stop();

        // 5. Indicator eta_K = 1 - |mean_K c| through a standalone Gauss
        //    evaluator (it never uses the quadrature factory), rebuilt each
        //    step because it stores the element set of mb at this call.
        sw.restart();
        gsExprEvaluator<real_t> ev;
        ev.options().setReal("quA", 2.0);
        ev.options().setInt ("quB", 1);
        ev.setIntegrationElements(mb);
        auto G = ev.getMap(mp);
        auto c = ev.getVariable(Cmp);
        ev.integralElWise(meas(G) * c);
        const std::vector<real_t> cInt = ev.elementwise();
        ev.integralElWise(meas(G));
        const std::vector<real_t> areas = ev.elementwise();
        GISMO_ENSURE(static_cast<index_t>(cInt.size()) == nElem && areas.size() == cInt.size(),
                     "The evaluator returned " << cInt.size() << " element values for " << nElem << " elements.");
        std::vector<real_t> eta(nElem);
        for (index_t e = 0; e < nElem; ++e)
            eta[e] = 1 - math::abs(cInt[e] / areas[e]);
        const std::vector<ElemData> elems = gatherElements(mb.basis(0), n0, eta, numLevels);
        index_t maxLevel = 0;
        for (const ElemData & d : elems)
            maxLevel = std::max(maxLevel, d.level);
        const double tIndicator = sw.stop();

        // 6. Marking by thresholds. A refined element goes straight to level
        //    L. A parent is coarsened only when all four children are active
        //    elements of one level with eta < tauC (unrefineElements clears the
        //    whole parent box, so a partial family must not be marked).
        sw.restart();
        std::vector<index_t> refBoxes, crsBoxes;
        index_t nAbove = 0, nRef = 0, nCrs = 0;
        std::map<std::array<index_t,3>, real_t> etaOf;
        for (const ElemData & d : elems)
            etaOf[{{d.level, d.i, d.j}}] = d.eta;
        for (const ElemData & d : elems)
        {
            if (d.eta > tauR)
            {
                ++nAbove;
                if (d.level < numLevels)
                {
                    const index_t s = numLevels - d.level;
                    refBoxes.push_back(numLevels);
                    refBoxes.push_back(d.i << s);       refBoxes.push_back(d.j << s);
                    refBoxes.push_back((d.i + 1) << s); refBoxes.push_back((d.j + 1) << s);
                    ++nRef;
                }
            }
            if (d.level >= 1 && d.i % 2 == 0 && d.j % 2 == 0)
            {
                bool family = true;
                for (index_t a = 0; a < 2 && family; ++a)
                    for (index_t b = 0; b < 2 && family; ++b)
                    {
                        const auto child = etaOf.find({{d.level, d.i + a, d.j + b}});
                        family = (child != etaOf.end() && child->second < tauC);
                    }
                if (family)
                {
                    crsBoxes.push_back(d.level - 1);
                    crsBoxes.push_back(d.i / 2);       crsBoxes.push_back(d.j / 2);
                    crsBoxes.push_back(d.i / 2 + 1);   crsBoxes.push_back(d.j / 2 + 1);
                    nCrs += 4;
                }
            }
        }
        const double tMark = sw.stop();
        const bool meshChanged = (nRef > 0 || nCrs > 0);

        // 7. Adapt the basis in place: refinement first, then coarsening.
        sw.restart();
        if (!refBoxes.empty())
            mb.basis(0).refineElements(refBoxes);
        if (!crsBoxes.empty())
            mb.basis(0).unrefineElements(crsBoxes);
        const double tAdapt = sw.stop();

        double tInit = 0, tProj = 0, tBuild = 0;
        const index_t levelsOld = maxLevel;
        if (meshChanged)
        {
            // 8. Re-initialise the assemblers on the new mesh.
            sw.restart();
            for (CHAssembler * a : assemblers)
                a->initialize();
            tInit = sw.stop();

            // 9. L2 projection of c and dc/dt. Both the old field and the new
            //    space are piecewise polynomial on dyadic meshes of level <=
            //    levEv, so the finest tensor level integrates the product
            //    exactly. The Clamped DoF matching is restored by the
            //    constructSolution round trip.
            sw.restart();
            const index_t levEv = std::max(levelsOld, maxElementLevel(mb.basis(0)));
            gsMultiBasis<real_t> evbasis;
            evbasis.addBasis(static_cast<const gsHTensorBasis<2,real_t> &>(mb.basis(0)).tensorLevel(levEv).clone().release());
            const gsOptionList quOpts = asmMain.options().getGroup("ExprAssembler");
            gsMatrix<real_t> coefsC, coefsdC;
            gsL2Projection<real_t>::project(mb, evbasis, mp, Cmp.patch(0),  coefsC,  gsBoundaryConditions<real_t>(), quOpts);
            gsL2Projection<real_t>::project(mb, evbasis, mp, dCmp.patch(0), coefsdC, gsBoundaryConditions<real_t>(), quOpts);
            gsMultiPatch<real_t> Cproj, dCproj;
            Cproj.addPatch(mb.basis(0).makeGeometry(give(coefsC)));
            dCproj.addPatch(mb.basis(0).makeGeometry(give(coefsdC)));
            asmMain.constructSolution(Cproj, Cold);
            asmMain.constructSolution(dCproj, dCold);
            asmMain.constructSolution(Cold, Cmp);
            asmMain.constructSolution(dCold, dCmp);
            tProj = sw.stop();

            // 10. New leaf factory on the new mesh; it serves the next step.
            if (useLeaf)
                tBuild = buildLeaf();

            // 11. The mass matrix is re-assembled with the new rule and mesh;
            //     its time lands in the stepper's mass bucket.
            stepper.assembleMassMatrix();
        }
        const index_t nElemNew = mb.basis(0).numElements();

        // 12. Validity of the adaptation: every pre-adapt element with
        //     eta > tauR must now be at level L at its centre.
        sw.restart();
        index_t nUnresolved = 0;
        const gsHTensorBasis<2,real_t> & hb = static_cast<const gsHTensorBasis<2,real_t> &>(mb.basis(0));
        for (const ElemData & d : elems)
            if (d.eta > tauR)
            {
                gsMatrix<real_t> ctr(2,1);
                ctr << d.cx, d.cy;
                if (hb.getLevelAtPoint(ctr) != numLevels)
                    ++nUnresolved;
            }
        tDiag += sw.stop();

        const gsCahnHilliardStepTimings<real_t> T1 = stepper.timings();
        GISMO_ENSURE(T1.steps - T0.steps == 1, "The stepper counted " << T1.steps - T0.steps << " steps in one step.");
        const double tWall = swWall.stop();

        row.set("step", step);
        row.set("t", t);
        row.set("converged", index_t(1));
        row.set("resRatio", status.residualRatio);
        row.set("nElem", nElem);
        row.set("nDof", nDof);
        row.set("nBasis", nBasis);
        row.set("maxLevel", maxLevel);
        row.set("nAboveTauR", nAbove);
        row.set("nRef", nRef);
        row.set("nCrs", nCrs);
        row.set("meshChanged", index_t(meshChanged));
        row.set("nElemNew", nElemNew);
        row.set("nUnresolved", nUnresolved);
        row.set("tWall", tWall);
        row.set("tStep", tStep);
        setStepperColumns(T1);
        row.set("tIndicator", tIndicator);
        row.set("tMark", tMark);
        row.set("tAdapt", tAdapt);
        row.set("tInit", tInit);
        row.set("tProj", tProj);
        row.set("tBuild", tBuild);
        row.set("tDiag", tDiag);
        row.set("energy", energy);
        row.set("mass", mass);
        row.set("cmin", bounds.cmin);
        row.set("cmax", bounds.cmax);
        row.set("certified", index_t(bounds.certified));
        if (useLeaf)
            setLeafColumns();
        // tShadow was set in step 3 (it is not part of tWall's other buckets).
        out.line(row.str());

        gsInfo << "step " << step << ": t=" << t << ", newton=" << T1.newtonIterations - T0.newtonIterations
               << ", nElem " << nElem << " -> " << nElemNew << " (ref " << nRef << ", crs " << nCrs
               << "), nUnresolved=" << nUnresolved << ", wall=" << tWall << " s\n";
    }

    out.line(headerLine("status", std::string("completed")));
    return EXIT_SUCCESS;
}

#else

int main()
{
    gsInfo << "adaptive_ch_hpatch needs the gsCahnHilliard module (GISMO_OPTIONAL=gsCahnHilliard).\n";
    return EXIT_SUCCESS;
}

#endif
