/** @file hierarchical_ch_hpatch.cpp

    @brief Assembly benchmark for the Cahn-Hilliard residual and Jacobian on an
    interface-adapted 2D THB-spline mesh: element-wise Gauss-Legendre quadrature
    (the default of gsCahnHilliardAssembler) against gsHPatchRule (one reduced
    patch rule per kd leaf, or per merged box of leaves).

    One case is run per call. The unit square is refined in bands around a
    diffuse interface (a union of circles, or a line): level l+1 is added on every
    level-l cell that lies closer than beta*h_l to the interface, for
    l = 0, ..., L-1, with h_l the cell size of level l. The state is
        c(x) = tanh( d(x) / (sqrt(2) * eps) ),   eps = 2 h_L,
    with d the signed distance to the interface, quasi-interpolated onto the THB
    space and projected onto the space of the strong Clamped boundary condition
    (dc/dn = 0 on all four sides, imposed by DoF matching).

    Two gsCahnHilliardAssembler objects (AroundZero convention, constant
    mobility, quartic potential) share the same THB basis. The Gauss arm uses the
    options of the module; the leaf arm gets a gsHPatchRule factory of degree 4p
    and regularity p-3, installed with setQuadratureFactory(). Per call,
        residual :  psi''(c) grad(w).grad(c) + lambda lapl(w) lapl(c),
        Jacobian :  the matching bilinear form,
    each one assemble() on an initialised assembler (assembleResidual and
    assembleJacobian with the gsMultiPatch state). Times in seconds:
      tBuild : the makeFactory call of the leaf arm (partition, cache, rule),
      tRes, tJac : one residual / one Jacobian assembly, minimum over --reps,
      tAsm   : minimum over --reps of (tRes + tJac) of the same repetition.
    The arms alternate within the repetition loop, Gauss first, and the state is
    built outside every timed region. R = tAsmLeaf / tAsmGauss. The driver also
    checks that both arms assemble the same residual and Jacobian
    (relDiffRes, relDiffJac, relative 2-norm and Frobenius norm); both are
    expected to be at rounding level, since the rule is exact for the integrands
    on every element of the affine geometry.

    The geometry is the bilinear unit square (one element, affine map), separate
    from the degree-p space basis.

    At p=2 the target regularity p-3 = -1 makes every direction of every cell tie
    with Gauss, so the leaf arm is plain (2p+1)-point Gauss and not a reduced
    rule: numCells < numLeaves is then no saving, and relDiffRes and relDiffJac
    are exactly 0.

    The CSV row has the columns
      caseId,family,k,theta,p,L,n0,beta,eps,lambda,nElem,nDof,nBasis,maxLevel,
      numLeaves,numCells,ptsGauss,ptsLeaf,tBuild,tResGauss,tJacGauss,tAsmGauss,
      tResLeaf,tJacLeaf,tAsmLeaf,R,relDiffRes,relDiffJac,nFallback,reps,threads,
      host,cpu
    where k is empty for the line family and theta is empty for circles; nDof is
    the number of free DoFs after the Clamped matching; nBasis is the size of the
    THB basis; numLeaves and numCells are the kd leaves of the domain and the
    cells of the leaf rule (boxes after merging); ptsGauss and ptsLeaf are the
    quadrature points assemble() evaluates per integral; nFallback is the number
    of cells with a direction on Gauss nodes because a 1D Newton solve failed;
    lambda is the value used.

    The timings are serial by design, so the driver must be run with
    OMP_NUM_THREADS=1; otherwise it exits with code 2.
    It needs the gsCahnHilliard module (GISMO_OPTIONAL=gsCahnHilliard); without
    it the program only prints a message.

    Usage:
        OMP_NUM_THREADS=1 ./hierarchical_ch_hpatch [--family circles|line] [-k 1|2|3 | --ncircles 1|2|3]
            [--theta deg] [-P deg] [-l levels] [--n0 n] [--beta b] [--lambda l]
            [--reps n] [--caseId name] [--csv file]
    Circles: k=1 one circle (centre (0.5,0.5), r 0.3); k=2 two circles of
    radius 0.15; k=3 three circles of radius 0.15. Line: through (0.5,0.5) at
    angle theta (degrees, counter-clockwise from the x-axis). --lambda -1 means
    lambda = eps^2.

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

typedef gsCahnHilliardAssembler<real_t, CH::Convention::AroundZero,
                                CH::Mobility::Constant, CH::Potential::Quartic> CHAssembler;

// ============================================================
// Interface geometry
// ============================================================

/// Signed distance to the interface (negative inside) and the exact minimum of
/// its modulus over an axis-aligned box.
class Interface
{
public:
    /// Circles: centres cx, cy and radii r. The circles are disjoint, so the
    /// band {|d| < delta} is the union of the bands of the single circles.
    Interface(const std::vector<real_t> & cx, const std::vector<real_t> & cy,
              const std::vector<real_t> & r)
    : m_cx(cx), m_cy(cy), m_r(r), m_isLine(false), m_nx(0), m_ny(0) { }

    /// Line through (0.5,0.5) with the normal (-sin theta, cos theta).
    explicit Interface(real_t thetaDeg)
    : m_isLine(true)
    {
        const real_t th = thetaDeg * EIGEN_PI / 180;
        m_nx = -std::sin(th);
        m_ny =  std::cos(th);
    }

    /// Signed distance at (x,y); for several circles the smallest in modulus.
    real_t distance(real_t x, real_t y) const
    {
        if (m_isLine)
            return (x - 0.5) * m_nx + (y - 0.5) * m_ny;
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
    real_t minAbsDistance(real_t x0, real_t x1, real_t y0, real_t y1) const
    {
        if (m_isLine)
        {
            // d is affine: the minimum of |d| is 0 if the corner values change
            // sign (or vanish), else the smallest corner modulus.
            const real_t v[4] = { distance(x0,y0), distance(x1,y0),
                                  distance(x0,y1), distance(x1,y1) };
            const bool allPos = v[0] > 0 && v[1] > 0 && v[2] > 0 && v[3] > 0;
            const bool allNeg = v[0] < 0 && v[1] < 0 && v[2] < 0 && v[3] < 0;
            if (!allPos && !allNeg)
                return 0;
            real_t m = math::abs(v[0]);
            for (index_t i = 1; i < 4; ++i)
                m = std::min(m, math::abs(v[i]));
            return m;
        }
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
    bool   m_isLine;
    real_t m_nx, m_ny;
};

/// The diffuse-interface state c = tanh(d / (sqrt(2) eps)) on the unit square.
/// For the quartic potential (c^2-1)^2/4 and lambda = eps^2 this is the 1D
/// equilibrium profile across a flat interface.
class DiffuseState : public gsFunction<real_t>
{
public:
    DiffuseState(const Interface & interface, real_t eps)
    : m_interface(interface), m_eps(eps) { }

    short_t domainDim() const { return 2; }
    short_t targetDim() const { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const
    {
        result.resize(1, u.cols());
        for (index_t i = 0; i < u.cols(); ++i)
            result(0,i) = std::tanh(m_interface.distance(u(0,i), u(1,i)) / (math::sqrt(real_t(2)) * m_eps));
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
/// constructSolution() is part of the timing; switching to the gsMatrix
/// overloads (free DoF vector) only changes these two functions.
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

const char * csvHeader()
{
    return "caseId,family,k,theta,p,L,n0,beta,eps,lambda,nElem,nDof,nBasis,maxLevel,"
           "numLeaves,numCells,ptsGauss,ptsLeaf,tBuild,tResGauss,tJacGauss,tAsmGauss,"
           "tResLeaf,tJacLeaf,tAsmLeaf,R,relDiffRes,relDiffJac,nFallback,reps,threads,"
           "host,cpu";
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

// ============================================================
// Main
// ============================================================

int main(int argc, char * argv[])
{
    std::string family = "circles";
    std::string caseId = "adhoc";
    std::string csvFile;
    index_t k = 1;
    index_t p = 2;
    index_t numLevels = 2;
    index_t n0 = 16;
    index_t reps = 5;
    real_t  theta = 0;
    real_t  beta = 2;
    real_t  lambdaIn = -1;

    gsCmdLine cmd("Cahn-Hilliard residual and Jacobian assembly on band-refined THB meshes: "
                  "gsHPatchRule vs Gauss.");
    cmd.addString("",  "family", "Interface family: circles or line", family);
    cmd.addInt   ("k", "ncircles", "Number of circles (1, 2 or 3)", k);
    cmd.addReal  ("",  "theta",  "Angle of the line in degrees", theta);
    cmd.addInt   ("P", "deg",    "Polynomial degree of the space basis (>= 2)", p);
    cmd.addInt   ("l", "levels", "Number of refinement levels", numLevels);
    cmd.addInt   ("",  "n0",     "Elements per direction of the base mesh", n0);
    cmd.addReal  ("",  "beta",   "Band half-width in cells of the refined level", beta);
    cmd.addReal  ("",  "lambda", "Interface parameter lambda (-1: eps^2)", lambdaIn);
    cmd.addInt   ("",  "reps",   "Repetitions of the timed assemblies; the minimum is reported", reps);
    cmd.addString("",  "caseId", "Label of this case in the CSV row", caseId);
    cmd.addString("",  "csv",    "Write the header and the row of this case to this file", csvFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // The timings compare single-thread assembly; they are serial by design.
    if (omp_get_max_threads() > 1)
    {
        gsWarn << "This benchmark requires OMP_NUM_THREADS=1 (found "
               << omp_get_max_threads() << " threads).\n";
        return 2;
    }

    if (p < 2)            { gsInfo << "Error: the degree must be >= 2.\n"; return EXIT_FAILURE; }
    if (numLevels < 0)    { gsInfo << "Error: levels must be >= 0.\n"; return EXIT_FAILURE; }
    if (n0 < 1)           { gsInfo << "Error: n0 must be >= 1.\n"; return EXIT_FAILURE; }
    if (beta <= 0)        { gsInfo << "Error: beta must be > 0.\n"; return EXIT_FAILURE; }
    if (reps < 1)         { gsInfo << "Error: reps must be >= 1.\n"; return EXIT_FAILURE; }
    if (family != "circles" && family != "line")
    {
        gsInfo << "Error: unknown family '" << family << "'. Expected circles or line.\n";
        return EXIT_FAILURE;
    }
    if (family == "circles" && (k < 1 || k > 3))
    {
        gsInfo << "Error: k must be 1, 2 or 3.\n";
        return EXIT_FAILURE;
    }
    if (lambdaIn != -1 && lambdaIn <= 0)
    {
        gsInfo << "Error: lambda must be positive (or -1 for eps^2).\n";
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

    // ---- Interface and mesh -------------------------------------------------
    std::vector<real_t> cx, cy, cr;
    if (k == 1)      { cx = {0.5};           cy = {0.5};           cr = {0.3};  }
    else if (k == 2) { cx = {0.3, 0.7};      cy = {0.5, 0.5};      cr = {0.15, 0.15}; }
    else             { cx = {0.3, 0.7, 0.5}; cy = {0.3, 0.3, 0.72}; cr = {0.15, 0.15, 0.15}; }
    const Interface interface = (family == "line") ? Interface(theta) : Interface(cx, cy, cr);

    // Interface width: twice the cell size of the finest level.
    const real_t hL  = real_t(1) / (real_t(n0) * std::pow(real_t(2), real_t(numLevels)));
    const real_t eps = 2 * hL;
    const real_t lambda = (lambdaIn == -1) ? eps * eps : lambdaIn;

    // Geometry: the bilinear unit square, one element (affine map).
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare());
    mp.addAutoBoundaries();

    // Space basis: degree p, maximal smoothness, n0 x n0 elements, then the
    // band refinement.
    gsKnotVector<real_t> kv(0, 1, n0 - 1, p + 1, 1);
    gsTensorBSplineBasis<2,real_t> tbsb(kv, kv);
    gsTHBSplineBasis<2,real_t> thb(tbsb);
    refineBands(thb, interface, n0, numLevels, beta);

    // The assemblers borrow this multi-basis as space and integration basis, so
    // it is declared before them and never reassigned.
    gsMultiBasis<real_t> mb;
    mb.addBasis(thb.clone().release());
    mb.setTopology(mp);
    const index_t nElem    = mb.basis(0).numElements();
    const index_t nBasis   = mb.basis(0).size();
    const index_t maxLevel = thb.maxLevel();

    // Strong Clamped condition (dc/dn = 0) on all four sides.
    gsBoundaryConditions<real_t> bc;
    for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::clamped, 0);
    bc.setGeoMap(mp);

    // ---- Assemblers ---------------------------------------------------------
    // The cache of the leaf rule must outlive the factory and all its copies.
    gsHPatchRuleCache<real_t> cache;

    CHAssembler gauss(mp, mb, bc);
    CHAssembler leaf (mp, mb, bc);
    for (CHAssembler * a : { &gauss, &leaf })
    {
        a->options().setReal  ("Lambda", lambda);
        a->options().setSwitch("AssembleWeakBCs", false);
        // No separate integration basis: the space basis is also used for
        // the integration.
        a->options().setSwitch("ExprAssembler.SameElement", true);
        a->initialize();
    }
    const index_t nDof = gauss.numDofs();

    // ---- State --------------------------------------------------------------
    // Quasi-interpolant of the tanh profile, projected onto the Clamped space by
    // keeping only the free DoFs. Built once; both arms assemble at this state.
    const DiffuseState fun(interface, eps);
    gsMatrix<real_t> coefs;
    gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0), fun, coefs);
    gsMultiPatch<real_t> C0;
    C0.addPatch(mb.basis(0).makeGeometry(give(coefs)));
    gsMatrix<real_t> Cvec;
    gauss.constructSolution(C0, Cvec);
    gsMultiPatch<real_t> C;
    gauss.constructSolution(Cvec, C);

    // ---- Leaf rule ----------------------------------------------------------
    // One rule for every integral of the residual and Jacobian: degree 4p,
    // regularity p-3 (-1 at p=2 means no continuity assumed).
    gsStopwatch sw;
    gsExprAssembler<real_t>::QuadratureFactory factory =
        gsHPatchRule<real_t>::makeFactory(mb, 4*p, p - 3, false, &cache,
                                          gsHPatchRule<real_t>::MergeMode::best);
    const double tBuild = sw.stop();
    leaf.setQuadratureFactory(factory);

    // A handle of the factory's rule: it shares the data of the rule the
    // assembler uses and serves the point counts and diagnostics.
    const gsOptionList optLeaf = leaf.options().getGroup("ExprAssembler");
    const gsQuadRule<real_t>::uPtr handle = factory(mb.basis(0), optLeaf, 0, -1);
    const gsHPatchRule<real_t> & hpatchRule = dynamic_cast<const gsHPatchRule<real_t> &>(*handle);
    const gsHPatchRule<real_t>::PartitionInfo & info = hpatchRule.partitionInfo();
    const index_t numLeaves  = info.numLeaves;
    const index_t numCells   = info.numCells;
    const index_t nFallback  = hpatchRule.numFallbackCells();
    const index_t ptsLeaf    = countPoints(mb.basis(0), hpatchRule);
    GISMO_ENSURE(ptsLeaf == hpatchRule.numPointsTotal(),
                 "The element-wise point count (" << ptsLeaf
                 << ") differs from the point count of the rule ("
                 << hpatchRule.numPointsTotal() << ").");

    // ---- Gauss rule of the module and its point count ------------------------
    const gsOptionList optGauss = gauss.options().getGroup("ExprAssembler");
    const gsQuadRule<real_t>::uPtr gaussRule = gsQuadrature::getPtr(mb.basis(0), optGauss);
    const index_t ptsGauss = countPoints(mb.basis(0), *gaussRule);
    GISMO_ENSURE(ptsGauss == nElem * (2*p + 1) * (2*p + 1),
                 "The Gauss point count (" << ptsGauss << ") differs from nElem*(2p+1)^2 = "
                 << nElem * (2*p + 1) * (2*p + 1) << ".");

    // ---- Timing -------------------------------------------------------------
    // Warm-up: first calls allocate and fill caches.
    callResidual(gauss, C); callJacobian(gauss, C);
    callResidual(leaf,  C); callJacobian(leaf,  C);

    ArmTimes tG, tL;
    double res, jac;
    for (index_t rep = 0; rep < reps; ++rep)
    {
        sw.restart(); callResidual(gauss, C); res = sw.stop();
        sw.restart(); callJacobian(gauss, C); jac = sw.stop();
        tG.update(res, jac);

        sw.restart(); callResidual(leaf, C);  res = sw.stop();
        sw.restart(); callJacobian(leaf, C);  jac = sw.stop();
        tL.update(res, jac);
    }
    const double R = tL.tAsm / tG.tAsm;

    // ---- Exactness ------------------------------------------------------------
    // Each arm holds the output of its last assembly.
    const gsMatrix<real_t> rG = gauss.rhs();
    const gsMatrix<real_t> rL = leaf.rhs();
    const gsSparseMatrix<real_t> JG = gauss.matrix();
    const gsSparseMatrix<real_t> JL = leaf.matrix();
    const gsSparseMatrix<real_t> dJ = JL - JG;
    const real_t relDiffRes = (rL - rG).norm() / rG.norm();
    const real_t relDiffJac = dJ.norm() / JG.norm();
    const bool valid = relDiffRes <= 1e-10 && relDiffJac <= 1e-10;

    // ---- Output ---------------------------------------------------------------
    gsInfo << "Case " << caseId << ": " << family;
    if (family == "circles") gsInfo << " k=" << k; else gsInfo << " theta=" << theta;
    gsInfo << ", p=" << p << ", L=" << numLevels << ", n0=" << n0 << ", beta=" << beta
           << ", eps=" << eps << ", lambda=" << lambda << "\n";
    gsInfo << "Mesh: nElem=" << nElem << ", nDof=" << nDof << ", nBasis=" << nBasis
           << ", maxLevel=" << maxLevel << "\n";
    gsInfo << "Leaf rule: leaves=" << numLeaves << ", cells=" << numCells
           << ", nFallback=" << nFallback << ", cache solves=" << cache.numSolves()
           << ", failures=" << cache.numFailures() << "\n";
    gsInfo << "Points: Gauss=" << ptsGauss << " (" << gaussRule->numNodes()
           << " per element), leaf=" << ptsLeaf << ", leaf/Gauss="
           << real_t(ptsLeaf) / ptsGauss << "\n";
    if (p == 2 && ptsLeaf != ptsGauss)
        gsInfo << "FINDING: p=2 point counts differ: leaf=" << ptsLeaf
               << ", Gauss=" << ptsGauss << "\n";
    gsInfo << "tBuild=" << tBuild << " s\n";
    gsInfo << "Gauss: tRes=" << tG.tRes << ", tJac=" << tG.tJac << ", tAsm=" << tG.tAsm << " s\n";
    gsInfo << "Leaf : tRes=" << tL.tRes << ", tJac=" << tL.tJac << ", tAsm=" << tL.tAsm << " s\n";
    gsInfo << "R = tAsmLeaf/tAsmGauss = " << R << "\n";
    gsInfo << "relDiffRes=" << relDiffRes << ", relDiffJac=" << relDiffJac
           << (valid ? "  VALID" : "  INVALID") << "\n";

    std::ostringstream row;
    row << std::scientific << std::setprecision(10);
    row << caseId << ',' << family << ',';
    if (family == "circles") row << k;
    row << ',';
    if (family == "line") row << theta;
    row << ',' << p << ',' << numLevels << ',' << n0 << ',' << beta << ',' << eps << ','
        << lambda << ',' << nElem << ',' << nDof << ',' << nBasis << ',' << maxLevel << ','
        << numLeaves << ',' << numCells << ',' << ptsGauss << ',' << ptsLeaf << ','
        << tBuild << ',' << tG.tRes << ',' << tG.tJac << ',' << tG.tAsm << ','
        << tL.tRes << ',' << tL.tJac << ',' << tL.tAsm << ',' << R << ','
        << relDiffRes << ',' << relDiffJac << ',' << nFallback << ',' << reps << ','
        << omp_get_max_threads() << ',' << hostName() << ',' << cpuModel();

    gsInfo << csvHeader() << '\n' << row.str() << '\n';
    if (csv.is_open())
        csv << row.str() << '\n';

    return EXIT_SUCCESS;
}

#else

int main()
{
    gsInfo << "hierarchical_ch_hpatch needs the gsCahnHilliard module (GISMO_OPTIONAL=gsCahnHilliard).\n";
    return EXIT_SUCCESS;
}

#endif
