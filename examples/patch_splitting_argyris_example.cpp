/** @file patch_splitting_argyris_example.cpp
    @brief Comprehensive Benchmark, Verification, and End-to-End ArgyrisSplines
           Biharmonic PDE Integration Suite for gsPatchSplitter<T>.

    Tests Mode A (existing gsMultiPatch splitting) and Mode B (gsPlanarDomain
    decomposition) across stretched rectangles, single-patch degenerate disks,
    L-shapes, T-shapes, smooth disks, plates with circular holes, and previously
    stalled concave/star/heptagon domains, followed by AS-G1 fitting and 4th-order
    biharmonic PDE convergence verification using gsAsG1Basis.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#include <gismo.h>
#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <gsModeling/gsASG1LocalFitter.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>
#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsReadFile.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace gismo;
using namespace gismo::expr;

namespace {

using T = real_t;
using Pt2 = gsPoint2D<T>;

gsBSpline<T>* makeLineBSpline(const Pt2& a, const Pt2& b, int deg = 3)
{
    int n = deg + 1;
    gsKnotVector<T> kv(0.0, 1.0, 0, n, 1, deg);
    gsMatrix<T> coefs(n, 2);
    for (int i = 0; i < n; ++i) {
        T t = static_cast<T>(i) / static_cast<T>(deg);
        coefs(i, 0) = (1 - t) * a.x + t * b.x;
        coefs(i, 1) = (1 - t) * a.y + t * b.y;
    }
    return new gsBSpline<T>(kv, coefs);
}

gsPlanarDomain<T>* makePolygonDomain(const std::vector<Pt2>& pts, int deg = 3)
{
    std::vector<gsCurve<T>*> curves;
    int n = static_cast<int>(pts.size());
    curves.reserve(n);
    for (int i = 0; i < n; ++i)
        curves.push_back(makeLineBSpline(pts[i], pts[(i + 1) % n], deg));
    return new gsPlanarDomain<T>(new gsCurveLoop<T>(curves));
}

gsMultiPatch<T> makeStretchedRectMultiPatch(T Lx, T Ly, int deg = 3)
{
    gsKnotVector<T> kv(0.0, 1.0, 0, deg + 1, 1, deg);
    gsMatrix<T> coefs((deg + 1) * (deg + 1), 2);
    for (int j = 0; j <= deg; ++j) {
        T v = static_cast<T>(j) / deg;
        for (int i = 0; i <= deg; ++i) {
            T u = static_cast<T>(i) / deg;
            coefs(i + j * (deg + 1), 0) = Lx * u;
            coefs(i + j * (deg + 1), 1) = Ly * v;
        }
    }
    gsMultiPatch<T> mp;
    mp.addPatch(typename gsGeometry<T>::uPtr(new gsTensorBSpline<2, T>(kv, kv, coefs)));
    mp.computeTopology();
    return mp;
}

// Single-patch cubic B-spline disk/ellipse with 4 degenerate 180-deg corners
gsMultiPatch<T> makeSinglePatchUnitDisk(int deg = 3)
{
    auto dom = std::unique_ptr<gsPlanarDomain<T>>(gsPlanarDomainGenerator<T>::preset("circle", 1.0, deg));
    const auto& outer = dom->outer();
    // outer has 4 quadrant B-splines in CCW order: south, east, north(rev), west(rev)
    gsBSpline<T> s0(*dynamic_cast<const gsBSpline<T>*>(&outer.curve(0)));
    gsBSpline<T> s1(*dynamic_cast<const gsBSpline<T>*>(&outer.curve(1)));
    gsBSpline<T> s2(*dynamic_cast<const gsBSpline<T>*>(&outer.curve(2)));
    gsBSpline<T> s3(*dynamic_cast<const gsBSpline<T>*>(&outer.curve(3)));
    s0.knots().affineTransformTo(0.0, 1.0);
    s1.knots().affineTransformTo(0.0, 1.0);
    s2.knots().affineTransformTo(0.0, 1.0);
    s3.knots().affineTransformTo(0.0, 1.0);
    s2.reverse(); // north: runs u: 0->1
    s3.reverse(); // west: runs v: 0->1

    auto harmonize = [](gsBSpline<T>& a, gsBSpline<T>& b) {
        gsKnotVector<T> ku = a.knots().knotUnion(b.knots());
        std::vector<T> ia, ib;
        ku.symDifference(a.knots(), ia);
        ku.symDifference(b.knots(), ib);
        if (!ia.empty()) a.insertKnots(ia.begin(), ia.end());
        if (!ib.empty()) b.insertKnots(ib.begin(), ib.end());
    };
    harmonize(s0, s2);
    harmonize(s3, s1);

    const int nu = s0.basis().size(), nv = s3.basis().size();
    gsMatrix<T> gru = s0.knots().greville(), grv = s3.knots().greville();
    const gsMatrix<T>& cS = s0.coefs();
    const gsMatrix<T>& cE = s1.coefs();
    const gsMatrix<T>& cN = s2.coefs();
    const gsMatrix<T>& cW = s3.coefs();
    gsMatrix<T> coefs(nu * nv, 2);
    for (int j = 0; j < nv; ++j) {
        T eta = grv(0, j);
        for (int i = 0; i < nu; ++i) {
            T xi = gru(0, i);
            for (int d = 0; d < 2; ++d) {
                T sn = (1 - eta) * cS(i, d) + eta * cN(i, d);
                T we = (1 - xi) * cW(j, d) + xi * cE(j, d);
                T cc = (1 - xi) * (1 - eta) * cS(0, d) + xi * (1 - eta) * cS(nu - 1, d)
                     + (1 - xi) * eta * cN(0, d) + xi * eta * cN(nu - 1, d);
                coefs(i + j * nu, d) = sn + we - cc;
            }
        }
    }

    gsMultiPatch<T> mp;
    mp.addPatch(typename gsGeometry<T>::uPtr(new gsTensorBSpline<2, T>(s0.knots(), s3.knots(), coefs)));
    mp.computeTopology();
    return mp;
}

gsPlanarDomain<T>* makeTShapeDomain(int deg = 3)
{
    // Symmetric T-shape: stem [-0.5, 0.5] x [-1.5, 0.5], top bar [-1.5, 1.5] x [0.5, 1.5]
    std::vector<Pt2> pts = {
        {-0.5, -1.5},
        { 0.5, -1.5},
        { 0.5,  0.5},
        { 1.5,  0.5},
        { 1.5,  1.5},
        {-1.5,  1.5},
        {-1.5,  0.5},
        {-0.5,  0.5}
    };
    return makePolygonDomain(pts, deg);
}

gsPlanarDomain<T>* makeLShapeDomain(int deg = 3)
{
    std::vector<Pt2> pts = {
        {-1.0, -1.0},
        { 1.0, -1.0},
        { 1.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  1.0},
        {-1.0,  1.0}
    };
    return makePolygonDomain(pts, deg);
}

T computeCurveArcLength(const gsGeometry<T>& crv, int nS = 128)
{
    gsMatrix<T> pr = crv.parameterRange();
    T u0 = pr(0, 0), u1 = pr(0, 1);
    gsMatrix<T> u(1, nS + 1), pts;
    for (int i = 0; i <= nS; ++i)
        u(0, i) = u0 + (u1 - u0) * (static_cast<T>(i) / nS);
    crv.eval_into(u, pts);
    T len = 0;
    for (int i = 0; i < nS; ++i)
        len += std::hypot(pts(0, i + 1) - pts(0, i), pts(1, i + 1) - pts(1, i));
    return len;
}

T computeDomainPerimeter(const gsPlanarDomain<T>& dom)
{
    T perim = 0;
    for (int l = 0; l < dom.numLoops(); ++l) {
        const auto& lp = dom.loop(l);
        for (int c = 0; c < lp.size(); ++c)
            perim += computeCurveArcLength(lp.curve(c));
    }
    return perim;
}

T computeMultiPatchBoundaryPerimeter(const gsMultiPatch<T>& mp)
{
    T perim = 0;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it) {
        auto bnd = mp.patch(it->patch).boundary(it->side());
        perim += computeCurveArcLength(*bnd);
    }
    return perim;
}

struct BenchMetrics {
    std::string name;
    int nPatches = 0;
    int nInterfaces = 0;
    int nEV = 0;
    T minAR = 1, maxAR = 1;
    T minAng = 0, maxAng = 180;
    T minDetJ = 0;
    T perimErr = 0;
    bool passed = false;
};

bool verifySplitterResult(
    const std::string& name,
    const gsPatchSplitter<T>& splitter,
    const gsMultiPatch<T>& mpIn,
    T expectedPerimeter,
    const gsPlanarDomain<T>* optDom,
    BenchMetrics& out)
{
    out.name = name;
    gsMultiPatch<T> mp(mpIn);

    if (mp.nBoxes() < 1) {
        std::cerr << "  [FAIL] " << name << ": mp.nBoxes() < 1\n";
        return false;
    }

    if (!mp.computeTopology(1e-4, false, false)) {
        std::cerr << "  [FAIL] " << name << ": computeTopology failed\n";
        return false;
    }

    // Verify 0 hanging nodes: every side is either a boundary or a 1-to-1 interface
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        for (int s = 1; s <= 4; ++s) {
            patchSide ps(p, static_cast<boxSide>(s));
            bool isBnd = mp.isBoundary(ps);
            bool isIfc = mp.isInterface(ps);
            if (isBnd == isIfc) {
                std::cerr << "  [FAIL] " << name << ": patch " << p << " side " << s
                          << " isBnd=" << isBnd << " isIfc=" << isIfc << " (hanging node / topology mismatch)\n";
                return false;
            }
        }
    }

    // Perimeter conservation check
    T mpPerim = computeMultiPatchBoundaryPerimeter(mp);
    out.perimErr = std::abs(mpPerim - expectedPerimeter);
    if (out.perimErr > 2e-2 * std::max<T>(1.0, expectedPerimeter)) {
        std::cerr << "  [FAIL] " << name << ": boundary perimeter mismatch ("
                  << mpPerim << " vs expected " << expectedPerimeter << ")\n";
        return false;
    }

    if (optDom) {
        T bdrDist = 0;
        if (!gsPatchSplitter<T>::verifyBoundaryMatch(*optDom, mp, 2e-3, &bdrDist)) {
            std::cerr << "  [FAIL] " << name << ": boundary Hausdorff mismatch (" << bdrDist << ")\n";
            return false;
        }
    }

    // Dense 20x20 sampling per patch: det(J) > 0 everywhere
    T minDet20 = std::numeric_limits<T>::max();
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const gsMatrix<T> supp = mp.patch(p).support();
        const T u0 = supp(0, 0), u1 = supp(0, 1);
        const T v0 = supp(1, 0), v1 = supp(1, 1);
        gsMatrix<T> uv(2, 20 * 20), J;
        int idx = 0;
        for (int j = 0; j < 20; ++j) {
            T v = v0 + (v1 - v0) * (static_cast<T>(j) / 19);
            for (int i = 0; i < 20; ++i) {
                T u = u0 + (u1 - u0) * (static_cast<T>(i) / 19);
                uv(0, idx) = u;
                uv(1, idx) = v;
                ++idx;
            }
        }
        mp.patch(p).jacobian_into(uv, J);
        for (index_t k = 0; k < uv.cols(); ++k) {
            T detJ = J(0, 2 * k) * J(1, 2 * k + 1) - J(0, 2 * k + 1) * J(1, 2 * k);
            if (detJ < minDet20) minDet20 = detJ;
        }
    }
    if (minDet20 <= T(1e-9)) {
        std::cerr << "  [FAIL] " << name << ": min det(J) on 20x20 grid = " << minDet20 << " <= 0\n";
        return false;
    }

    // Verify ArgyrisSplines degree & knot span requirements
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        if (mp.basis(p).degree(0) < 3 || mp.basis(p).degree(1) < 3) {
            std::cerr << "  [FAIL] " << name << ": patch " << p << " degree < 3\n";
            return false;
        }
        if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p))) {
            if (tb->knots(0).unique().size() - 1 < 3 || tb->knots(1).unique().size() - 1 < 3) {
                std::cerr << "  [FAIL] " << name << ": patch " << p << " has < 3 knot spans\n";
                return false;
            }
        }
    }

    auto qr = splitter.analyzeQuality(mp);
    out.nPatches = qr.nPatches;
    out.nInterfaces = qr.nInterfaces;
    out.nEV = qr.nExtraordinaryVertices;
    out.minAR = qr.minAspectRatio;
    out.maxAR = qr.maxAspectRatio;
    out.minAng = qr.minAngleDeg;
    out.maxAng = qr.maxAngleDeg;
    out.minDetJ = minDet20;
    out.passed = qr.isValid && !qr.hasHangingNodes && qr.isConforming;
    return out.passed;
}

// End-to-End ArgyrisSplines Biharmonic PDE Test: Delta^2 u = f, u = sin(pi x) cos(pi y)
bool runArgyrisBiharmonicTest(
    const std::string& label,
    const gsMultiPatch<T>& inputMp,
    T& finalL2Err,
    T& finalH2Err,
    T& lastL2Rate,
    T& lastH2Rate)
{
    gsMultiPatch<T> target(inputMp);
    gsASG1LocalFitter<T>::unifyPatchOrientations(target);
    target.computeTopology(1e-4, false, false);

    // Ensure AS-G1 compatibility via gsASG1LocalFitter (p=3, r=1, k=2)
    gsASG1LocalFitter<T> fitter(target, 3, 1, 2);
    fitter.setMode(gsASG1LocalFitter<T>::FittingMode::Global);
    fitter.fit();
    gsMultiPatch<T> mp = fitter.result();
    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(1e-4, false, false);

    if (!fitter.checkASG1(1e-3)) {
        std::cerr << "  [FAIL] Argyris PDE (" << label << "): AS-G1 fitting residual too large ("
                  << fitter.maxAlphaResidual() << ", " << fitter.maxBetaResidual() << ")\n";
        return false;
    }

    gsFunctionExpr<T> exact_u("sin(pi*x)*cos(pi*y)", 2);
    gsFunctionExpr<T> exact_grad("pi*cos(pi*x)*cos(pi*y)", "-pi*sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<T> exact_hess("-pi^2*sin(pi*x)*cos(pi*y)", "-pi^2*cos(pi*x)*sin(pi*y)",
                                 "-pi^2*cos(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*cos(pi*y)", 2);
    gsFunctionExpr<T> rhs_f("4*pi^4*sin(pi*x)*cos(pi*y)", 2);

    gsConstantFunction<T> zero;
    gsBoundaryConditions<T> bc;
    const std::string bcLabel = "ValuesAndDerivatives";
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
        bc.add(it->patch, it->side(), bcLabel, zero);

    gsMatrix<T> gd = computeGluingData(mp, 1e-3, 0);
    std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, bcLabel);
    std::vector<cornerCondition> cc;
    gsMatrix<T> cornerDataForPatches = gsMatrix<T>::Zero(mp.nPatches(), 8);
    for (size_t i = 0; i < vertices.size(); ++i) {
        gsVector<T> normal = getOuterNormalDerivative(mp, vertices[i]);
        for (size_t j = 0; j < vertices[i].size(); ++j)
            cornerDataForPatches.block(vertices[i][j].patch, 2 * (vertices[i][j].m_index - 1), 1, 2) = normal.transpose();
        cc.push_back(cornerCondition{
            vertices[i][0],
            normal.norm() < T(1e-6) ? cornerConditionType::all : cornerConditionType::valuesNormals
        });
    }
    for (size_t i = 0; i < mp.nPatches(); ++i)
        for (index_t j = 0; j < 4; ++j)
            if (cornerDataForPatches.block(i, 2 * j, 1, 2).norm() < T(1e-6)) {
                cornerDataForPatches(i, 2 * j) = T(1);
                cornerDataForPatches(i, 2 * j + 1) = T(0);
            }

    std::cout << "\n--- ArgyrisSplines Biharmonic Convergence Test: " << label
              << " (P=" << mp.nPatches() << ", IF=" << mp.nInterfaces() << ") ---\n";
    std::cout << std::setw(5) << "ref" << std::setw(10) << "N_free"
              << std::setw(14) << "L2 Error" << std::setw(10) << "L2 Rate"
              << std::setw(14) << "H2 Error" << std::setw(10) << "H2 Rate" << "\n";

    T prevH = 0, prevL2 = 0, prevH2 = 0;
    lastL2Rate = 0;
    lastH2Rate = 0;

    for (index_t ref = 0; ref <= 2; ++ref) {
        gsMultiBasis<T> dbasis(mp);
        for (index_t i = 0; i < ref; ++i)
            dbasis.uniformRefine(1, 2);

        T h_max = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
            h_max = std::max(h_max, dbasis[i].getMaxCellLength());

        std::vector<gsArgyrisEmbedding<T>> argBasis;
        gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
        for (size_t i = 0; i < mp.nPatches(); ++i) {
            argBasis.push_back(deriveArgyrisBasisEmbedding(
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(dbasis[i]),
                gsMatrix<T>(gd.row(i)),
                gsMatrix<T>(cornerDataForPatches.row(i)),
                mp.patch(i)));
            argBasisGlobalCreator.set(i, i, argBasis[i].matrix);
        }
        gsSparseMatrix<T> argBasisGlobal = argBasisGlobalCreator;
        gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc, cc);
        gsSparseMatrix<T> T_global = argBasisGlobal * asEmbeddingMatrix<T>(mapper.size(), mapper.asVector()).transpose();

        const index_t nFree = mapper.freeSize(), nBnd = mapper.boundarySize();
        gsSparseMatrix<T> T_free = T_global.leftCols(nFree);
        gsSparseMatrix<T> T_bnd  = T_global.rightCols(nBnd);

        gsMatrix<T> sol_bnd;
        if (nBnd > 0) {
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());
            auto G_map   = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);
            auto u_c     = A.getCoeff(exact_u, G_map);
            auto g_c     = A.getCoeff(exact_grad, G_map);
            A.initSystem();
            A.assembleBdr(bc.get(bcLabel),
                          u_space * u_space.tr() * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (igrad(u_space, G_map) * unv(G_map)).tr() * meas(G_map),
                          u_space * u_c * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (g_c.tr() * unv(G_map)) * meas(G_map));
            gsSparseMatrix<T> M_bnd = T_bnd.transpose() * A.matrix() * T_bnd;
            gsMatrix<T>       F_bnd = T_bnd.transpose() * A.rhs();
            makeSparseCholeskySolver(M_bnd)->apply(F_bnd, sol_bnd);
        } else {
            sol_bnd.setZero(nBnd, 1);
        }

        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());
        auto G_map   = A.getMap(mp);
        auto u_space = A.getSpace(dbasis);
        auto f_c     = A.getCoeff(rhs_f, G_map);
        A.initSystem();
        A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                   u_space * f_c * meas(G_map));

        gsSparseMatrix<T> K_free = T_free.transpose() * A.matrix() * T_free;
        gsMatrix<T>       F_free = T_free.transpose() * (A.rhs() - A.matrix() * (T_bnd * sol_bnd));
        gsMatrix<T> sol_free;
        makeSparseCholeskySolver(K_free)->apply(F_free, sol_free);

        gsMatrix<T> solVec = T_free * sol_free + T_bnd * sol_bnd;
        gsMultiPatch<T> sol_mp;
        index_t off = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i) {
            index_t nC = dbasis[i].size();
            gsMatrix<T> c = solVec.middleRows(off, nC);
            sol_mp.addPatch(dbasis[i].makeGeometry(give(c)));
            off += nC;
        }

        gsExprEvaluator<T> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G_ev = ev.getMap(mp);
        auto u_num = ev.getVariable(sol_mp);
        auto u_ex  = ev.getVariable(exact_u, G_ev);
        auto h_ex  = reshape(ev.getVariable(exact_hess, G_ev), 2, 2);

        T l2Err = std::sqrt(ev.integral((u_num - u_ex).sqNorm() * meas(G_ev)));
        T h2Err = std::sqrt(ev.integral((ihess(u_num, G_ev) - h_ex).sqNorm() * meas(G_ev)));

        T l2Rate = 0, h2Rate = 0;
        if (ref > 0 && prevL2 > 0 && l2Err > 0 && prevH > h_max) {
            T lr = std::log(prevH / h_max);
            l2Rate = std::log(prevL2 / l2Err) / lr;
            h2Rate = std::log(prevH2 / h2Err) / lr;
            lastL2Rate = l2Rate;
            lastH2Rate = h2Rate;
        }

        std::cout << std::setw(5) << ref << std::setw(10) << nFree
                  << std::scientific << std::setprecision(3) << std::setw(14) << l2Err
                  << std::fixed << std::setprecision(2) << std::setw(10) << l2Rate
                  << std::scientific << std::setprecision(3) << std::setw(14) << h2Err
                  << std::fixed << std::setprecision(2) << std::setw(10) << h2Rate << "\n";

        prevH = h_max;
        prevL2 = l2Err;
        prevH2 = h2Err;
        finalL2Err = l2Err;
        finalH2Err = h2Err;
    }

    return (finalL2Err < T(1e-3)) && (lastL2Rate > T(2.5)) && (lastH2Rate > T(1.3));
}

} // namespace

int main(int argc, char* argv[])
{
    bool plot = false;
    std::string filesFrom = "";
    std::string outDir = "";
    gsCmdLine cmd("Benchmark & ArgyrisSplines Integration Suite for gsPatchSplitter.");
    cmd.addSwitch("plot", "Export ParaView .pvd/.vts files for each benchmark.", plot);
    cmd.addString("f", "files-from", "Run batch on XML paths listed in file (one path per line).", filesFrom);
    cmd.addString("o", "out-dir", "Output directory for PVD files and summary.csv in batch mode.", outDir);
    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    if (!filesFrom.empty()) {
        std::ifstream fin(filesFrom);
        if (!fin) {
            std::cerr << "Cannot open files list: " << filesFrom << "\n";
            return 1;
        }
        if (outDir.empty()) outDir = ".";
        mkdir(outDir.c_str(), 0755);

        auto sanitizeName = [](const std::string& p) {
            std::string s = p;
            const std::string prefix = "filedata/domain2d/";
            if (s.compare(0, prefix.size(), prefix) == 0)
                s = s.substr(prefix.size());
            if (s.size() > 4 && s.substr(s.size() - 4) == ".xml")
                s = s.substr(0, s.size() - 4);
            for (char& c : s)
                if (c == '/' || c == '\\' || c == ' ') c = '_';
            return s;
        };

        std::vector<BenchMetrics> batchSummary;
        int nTotal = 0, nPass = 0, nSkip = 0;
        std::string line;
        while (std::getline(fin, line)) {
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
                line.pop_back();
            if (line.empty() || line[0] == '#') continue;

            std::string label = sanitizeName(line);
            std::string pvdPrefix = outDir + "/" + label;
            gsFileData<T> fd(line);

            if (fd.has<gsPlanarDomain<T>>()) {
                ++nTotal;
                memory::unique_ptr<gsPlanarDomain<T>> dom = fd.getFirst<gsPlanarDomain<T>>();
                T perim = computeDomainPerimeter(*dom);
                gsPatchSplitter<T> splitter(*dom);
                gsMultiPatch<T> res = splitter.compute();
                BenchMetrics bm;
                bool ok = verifySplitterResult(label, splitter, res, perim, dom.get(), bm);
                if (ok) ++nPass;
                batchSummary.push_back(bm);
                if (res.nPatches() > 0)
                    gsWriteParaview(res, pvdPrefix, 200, true);
            } else if (fd.has<gsMultiPatch<T>>()) {
                ++nTotal;
                memory::unique_ptr<gsMultiPatch<T>> inMp = fd.getFirst<gsMultiPatch<T>>();
                inMp->computeTopology();
                T perim = computeMultiPatchBoundaryPerimeter(*inMp);
                gsPatchSplitter<T> splitter(*inMp);
                gsMultiPatch<T> res = splitter.compute();
                BenchMetrics bm;
                bool ok = verifySplitterResult(label + "[ModeA]", splitter, res, perim, nullptr, bm);
                if (ok) ++nPass;
                batchSummary.push_back(bm);
                if (res.nPatches() > 0)
                    gsWriteParaview(res, pvdPrefix, 200, true);
            } else if (fd.has<gsGeometry<T>>()) {
                ++nTotal;
                memory::unique_ptr<gsGeometry<T>> g = fd.getFirst<gsGeometry<T>>();
                gsMultiPatch<T> inMp;
                inMp.addPatch(give(g));
                inMp.computeTopology();
                T perim = computeMultiPatchBoundaryPerimeter(inMp);
                gsPatchSplitter<T> splitter(inMp);
                gsMultiPatch<T> res = splitter.compute();
                BenchMetrics bm;
                bool ok = verifySplitterResult(label + "[Geom]", splitter, res, perim, nullptr, bm);
                if (ok) ++nPass;
                batchSummary.push_back(bm);
                if (res.nPatches() > 0)
                    gsWriteParaview(res, pvdPrefix, 200, true);
            } else {
                ++nSkip;
            }
        }

        std::ofstream csv(outDir + "/summary.csv");
        csv << "domain,patches,interfaces,ev,min_ar,max_ar,min_ang_deg,max_ang_deg,min_det_j,status\n";
        std::cout << "\n==================================================================================================================\n";
        std::cout << "  filedata/domain2d Batch Results (" << nPass << " / " << nTotal << " PASSED, " << nSkip << " skipped)\n";
        std::cout << "==================================================================================================================\n";
        for (const auto& m : batchSummary) {
            std::cout << std::left << std::setw(54) << m.name
                      << " P=" << std::setw(4) << m.nPatches
                      << " IF=" << std::setw(4) << m.nInterfaces
                      << " EV=" << std::setw(3) << m.nEV
                      << " AR=[" << std::fixed << std::setprecision(2) << m.minAR << "," << m.maxAR << "]"
                      << " Ang=[" << std::setprecision(1) << m.minAng << "," << m.maxAng << "]"
                      << " minJ=" << std::scientific << std::setprecision(2) << m.minDetJ
                      << " " << (m.passed ? "PASS" : "FAIL") << "\n";
            csv << m.name << "," << m.nPatches << "," << m.nInterfaces << "," << m.nEV << ","
                << m.minAR << "," << m.maxAR << "," << m.minAng << "," << m.maxAng << ","
                << m.minDetJ << "," << (m.passed ? "PASS" : "FAIL") << "\n";
        }
        return (nPass == nTotal) ? 0 : 1;
    }

    std::vector<BenchMetrics> summary;
    int failures = 0;

    auto runModeA = [&](const std::string& name, const gsMultiPatch<T>& inMp, const std::string& pvdPrefix,
                        gsMultiPatch<T>* outMpSave = nullptr) {
        T perim = computeMultiPatchBoundaryPerimeter(inMp);
        gsPatchSplitter<T> splitter(inMp);
        gsMultiPatch<T> res = splitter.compute();
        BenchMetrics bm;
        bool ok = verifySplitterResult(name, splitter, res, perim, nullptr, bm);
        if (!ok) ++failures;
        summary.push_back(bm);
        if (plot && res.nPatches() > 0)
            gsWriteParaview(res, pvdPrefix, 400, true);
        if (outMpSave) *outMpSave = res;
    };

    auto runModeB = [&](const std::string& name, const gsPlanarDomain<T>& dom, const std::string& pvdPrefix,
                        int strategyOpt = 0, gsMultiPatch<T>* outMpSave = nullptr) {
        T perim = computeDomainPerimeter(dom);
        gsPatchSplitter<T> splitter(dom);
        splitter.options().setInt("Strategy", strategyOpt);
        gsMultiPatch<T> res = splitter.compute();
        BenchMetrics bm;
        bool ok = verifySplitterResult(name, splitter, res, perim, &dom, bm);
        if (!ok) ++failures;
        summary.push_back(bm);
        if (plot && res.nPatches() > 0)
            gsWriteParaview(res, pvdPrefix, 400, true);
        if (outMpSave) *outMpSave = res;
    };

    // 1. Mode A: 1 x 5 stretched rectangle -> ~5 near-unit squares
    gsMultiPatch<T> rectSplitMp, diskModeAMp;
    runModeA("ModeA_Rect_1x5", makeStretchedRectMultiPatch(5.0, 1.0), "split_argyris_rect1x5", &rectSplitMp);
    if (rectSplitMp.nPatches() < 4 || rectSplitMp.nPatches() > 6) {
        std::cerr << "  [FAIL] ModeA_Rect_1x5 expected 4..6 patches, got " << rectSplitMp.nPatches() << "\n";
        ++failures;
    }

    // 2. Mode A: Single-patch unit disk (180-deg corners) -> 5-patch butterfly O-grid
    runModeA("ModeA_Disk_1Patch", makeSinglePatchUnitDisk(), "split_argyris_disk_modeA", &diskModeAMp);
    if (diskModeAMp.nPatches() != 5) {
        std::cerr << "  [FAIL] ModeA_Disk_1Patch expected 5-patch butterfly, got " << diskModeAMp.nPatches() << "\n";
        ++failures;
    }

    // 3. Mode B: L-shaped domain -> 3 near-square patches
    gsMultiPatch<T> lshapeMp, diskModeBMp, holeMp;
    {
        std::unique_ptr<gsPlanarDomain<T>> dom(makeLShapeDomain());
        runModeB("ModeB_LShape", *dom, "split_argyris_lshape", 0, &lshapeMp);
    }

    // 4. Mode B: T-shaped domain -> minimal near-square blocks
    {
        std::unique_ptr<gsPlanarDomain<T>> dom(makeTShapeDomain());
        runModeB("ModeB_TShape", *dom, "split_argyris_tshape");
    }

    // 5. Mode B: Smooth unit disk (Cross-Field / O-grid) -> 5-patch butterfly
    {
        std::unique_ptr<gsPlanarDomain<T>> dom(gsPlanarDomainGenerator<T>::preset("circle", 1.0, 3));
        runModeB("ModeB_UnitDisk_CF", *dom, "split_argyris_disk_modeB", 1, &diskModeBMp);
    }

    // 6. Mode B: Plate with circular hole -> conforming 4-patch annulus
    {
        std::unique_ptr<gsPlanarDomain<T>> dom(gsPlanarDomainGenerator<T>::preset("plate_with_hole", 1.0, 3));
        runModeB("ModeB_PlateWithHole", *dom, "split_argyris_plate_hole", 0, &holeMp);
    }

    // 7. Mode B: Optional file-based stress tests (previously stalled / stretched channel)
    const std::vector<std::pair<std::string, std::string>> fileCases = {
        {"ModeB_SerpentineSDuct", "domain2d/noHolesCounterExamples/02_serpentine_s_duct.xml"},
        {"ModeB_StarPentagon",    "domain2d/non_bilinear_no_holes_1to7corners/star_pentagon_5corner_nonconvex.xml"},
        {"ModeB_ConcaveHeptagon", "domain2d/non_bilinear_no_holes_1to7corners/concave_heptagon_7corner_nonconvex.xml"},
        {"ModeB_Heptagon7",       "domain2d/noHolesCounterExamples/04_heptagon_7corner_convex.xml"}
    };
    for (const auto& fc : fileCases) {
        std::string found = gsFileManager::find(fc.second);
        if (found.empty()) {
            std::string rel = "filedata/" + fc.second;
            if (gsFileManager::fileExists(rel)) found = rel;
        }
        if (!found.empty()) {
            gsReadFile<T> rf(found);
            memory::unique_ptr<gsPlanarDomain<T>> pDom = rf;
            if (pDom)
                runModeB(fc.first, *pDom, "split_argyris_" + fc.first);
        }
    }

    // Print Summary Quality Metrics Table
    std::cout << "\n========================================================================================================\n";
    std::cout << "  gsPatchSplitter Benchmark Quality Summary\n";
    std::cout << "========================================================================================================\n";
    std::cout << std::left << std::setw(24) << "Benchmark"
              << std::setw(6)  << "P"
              << std::setw(6)  << "IF"
              << std::setw(5)  << "EV"
              << std::setw(18) << "AspectRatio"
              << std::setw(20) << "AngleRange(deg)"
              << std::setw(12) << "MinDetJ"
              << "Status\n";
    std::cout << std::string(104, '-') << "\n";
    for (const auto& m : summary) {
        std::ostringstream arSS, angSS;
        arSS << std::fixed << std::setprecision(2) << "[" << m.minAR << ", " << m.maxAR << "]";
        angSS << std::fixed << std::setprecision(1) << "[" << m.minAng << ", " << m.maxAng << "]";
        std::cout << std::left << std::setw(24) << m.name
                  << std::setw(6)  << m.nPatches
                  << std::setw(6)  << m.nInterfaces
                  << std::setw(5)  << m.nEV
                  << std::setw(18) << arSS.str()
                  << std::setw(20) << angSS.str()
                  << std::scientific << std::setprecision(2) << std::setw(12) << m.minDetJ
                  << (m.passed ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "========================================================================================================\n";

    // End-to-End ArgyrisSplines Biharmonic PDE Convergence Tests
    T l2_1 = 0, h2_1 = 0, rL2_1 = 0, rH2_1 = 0;
    if (!runArgyrisBiharmonicTest("ModeA_Rect_1x5", rectSplitMp, l2_1, h2_1, rL2_1, rH2_1)) ++failures;

    T l2_2 = 0, h2_2 = 0, rL2_2 = 0, rH2_2 = 0;
    if (!runArgyrisBiharmonicTest("ModeB_LShape", lshapeMp, l2_2, h2_2, rL2_2, rH2_2)) ++failures;

    T l2_3 = 0, h2_3 = 0, rL2_3 = 0, rH2_3 = 0;
    if (!runArgyrisBiharmonicTest("ModeB_UnitDisk", diskModeBMp, l2_3, h2_3, rL2_3, rH2_3)) ++failures;

    std::cout << "\nOverall patch_splitting_argyris_example status: "
              << (failures == 0 ? "ALL PASSED" : "FAILURES DETECTED") << "\n";
    return (failures == 0) ? 0 : 1;
}
