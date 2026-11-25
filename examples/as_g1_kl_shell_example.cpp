/** @file as_g1_kl_shell_example.cpp

    @brief Multi-patch and geometric multigrid Kirchhoff-Love shell driver on
           open AS-G1 curved surfaces and non-isoparametric NURBS shells,
           grounded in F. Hasanova, S. Takacs, T. Takacs (2026).

    Capabilities:
      1. Manufactured-solution h-refinement convergence study (exact L^2, H^1,
         covariant H^2, shell energy errors + EOCs) compared directly against
         the paper's Theorem 6 patch-local quasi-interpolator \Pi_{AS}.
      2. (p, k)-robustness sweep across p = 3..6 and 1 <= k <= p - 2.
      3. Multi-level geometric multigrid (nodal_gs, ring1, line_ring1, scms,
         scms_hybrid) and single-level Patch-Schwarz DD-PCG solver comparison.
      4. Physical Scordelis-Lo cylindrical roof benchmark on exact NURBS.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsModeling/gsAsG1ShellSolver.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace gismo;

namespace {

void mkdirRecursive(const std::string& path)
{
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i)
    {
        cur.push_back(path[i]);
        if (path[i] == '/' && cur.size() > 1)
            ::mkdir(cur.c_str(), 0755);
    }
    if (!cur.empty())
        ::mkdir(cur.c_str(), 0755);
}

/// Normalize all patch knot supports to [0, 1]^2 and embed 2D planar domains in z = 0.
gsMultiPatch<real_t> loadOpenSurface3D(const std::string& filename)
{
    gsMultiPatch<real_t> mp;
    gsFileData<real_t> fd(filename);
    if (fd.has<gsMultiPatch<real_t>>())
    {
        mp = *fd.getFirst<gsMultiPatch<real_t>>();
    }
    else if (fd.has<gsGeometry<real_t>>())
    {
        auto geoms = fd.getAll<gsGeometry<real_t>>();
        for (size_t i = 0; i < geoms.size(); ++i)
            mp.addPatch(give(*geoms[i]));
    }
    else
    {
        GISMO_ERROR("No geometry found in file: " << filename);
    }

    // If input is a 2D planar domain, embed it in R^3 as (x, y, 0) for plate tests
    if (mp.targetDim() == 2)
    {
        gsMultiPatch<real_t> mp3;
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSpline<2, real_t>& bsp =
                dynamic_cast<const gsTensorBSpline<2, real_t>&>(mp.patch(p));
            gsMatrix<real_t> c3 = gsMatrix<real_t>::Zero(bsp.coefs().rows(), 3);
            c3.leftCols(2) = bsp.coefs();
            mp3.addPatch(typename gsGeometry<real_t>::uPtr(new gsTensorBSpline<2, real_t>(bsp.basis(), c3)));
        }
        mp = give(mp3);
    }

    // Normalize supports to [0, 1]^2
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        if (auto* bsp = dynamic_cast<gsTensorBSpline<2, real_t>*>(&mp.patch(p)))
        {
            bsp->basis().knots(0).affineTransformTo(0.0, 1.0);
            bsp->basis().knots(1).affineTransformTo(0.0, 1.0);
        }
        else if (auto* nurbs = dynamic_cast<gsTensorNurbs<2, real_t>*>(&mp.patch(p)))
        {
            nurbs->basis().source().knots(0).affineTransformTo(0.0, 1.0);
            nurbs->basis().source().knots(1).affineTransformTo(0.0, 1.0);
        }
    }
    mp.computeTopology();
    return mp;
}

/// Smooth 3D analytic manufactured displacement u_ex(x,y,z) with exact gradient and Hessian.
/// Designed so both tangential membrane and normal bending strains are active on any open 3D shell.
struct SmoothManufacturedShellField
{
    bool plateOnlyZ = false;
    real_t freq = 1.0;

    void operator()(index_t comp, const gsVector<real_t>& xP,
                    real_t& val, gsVector<real_t>& g, gsMatrix<real_t>& H) const
    {
        const real_t x = freq * xP(0);
        const real_t y = freq * xP(1);
        const real_t z = freq * xP(2);
        const real_t f = freq;
        const real_t f2 = freq * freq;
        g.setZero(3);
        H.setZero(3, 3);

        if (plateOnlyZ)
        {
            if (comp != 2)
            {
                val = 0.0;
                return;
            }
            // Transverse bending wave w(x,y) = sin(f*x) * sin(f*y)
            const real_t sx = std::sin(x), cx = std::cos(x);
            const real_t sy = std::sin(y), cy = std::cos(y);
            val = sx * sy;
            g(0) = f * cx * sy;
            g(1) = f * sx * cy;
            H(0, 0) = -f2 * sx * sy;
            H(1, 1) = -f2 * sx * sy;
            H(0, 1) =  f2 * cx * cy;
            H(1, 0) = H(0, 1);
            return;
        }

        if (comp == 0)
        {
            // u_x = sin(x)*cos(y) + 0.3*cos(z)
            const real_t sx = std::sin(x), cx = std::cos(x);
            const real_t sy = std::sin(y), cy = std::cos(y);
            const real_t sz = std::sin(z), cz = std::cos(z);
            val = sx * cy + 0.3 * cz;
            g(0) =  f * cx * cy;
            g(1) = -f * sx * sy;
            g(2) = -0.3 * f * sz;
            H(0, 0) = -f2 * sx * cy;
            H(1, 1) = -f2 * sx * cy;
            H(0, 1) = -f2 * cx * sy;
            H(1, 0) = H(0, 1);
            H(2, 2) = -0.3 * f2 * cz;
        }
        else if (comp == 1)
        {
            // u_y = cos(x)*sin(y) - 0.25*sin(z)
            const real_t sx = std::sin(x), cx = std::cos(x);
            const real_t sy = std::sin(y), cy = std::cos(y);
            const real_t sz = std::sin(z), cz = std::cos(z);
            val = cx * sy - 0.25 * sz;
            g(0) = -f * sx * sy;
            g(1) =  f * cx * cy;
            g(2) = -0.25 * f * cz;
            H(0, 0) = -f2 * cx * sy;
            H(1, 1) = -f2 * cx * sy;
            H(0, 1) = -f2 * sx * cy;
            H(1, 0) = H(0, 1);
            H(2, 2) =  0.25 * f2 * sz;
        }
        else
        {
            // u_z = sin(x + 0.7*y)*cos(0.8*z) + 0.2*cos(x - y)
            const real_t a = x + 0.7 * y;
            const real_t b = 0.8 * z;
            const real_t d = x - y;
            const real_t sa = std::sin(a), ca = std::cos(a);
            const real_t sb = std::sin(b), cb = std::cos(b);
            const real_t sd = std::sin(d), cd = std::cos(d);

            val = sa * cb + 0.2 * cd;
            g(0) = f * (ca * cb - 0.2 * sd);
            g(1) = f * (0.7 * ca * cb + 0.2 * sd);
            g(2) = f * (-0.8 * sa * sb);

            H(0, 0) = f2 * (-sa * cb - 0.2 * cd);
            H(1, 1) = f2 * (-0.49 * sa * cb - 0.2 * cd);
            H(2, 2) = f2 * (-0.64 * sa * cb);
            H(0, 1) = f2 * (-0.7 * sa * cb + 0.2 * cd);
            H(1, 0) = H(0, 1);
            H(0, 2) = f2 * (-0.8 * ca * sb);
            H(2, 0) = H(0, 2);
            H(1, 2) = f2 * (-0.56 * ca * sb);
            H(2, 1) = H(1, 2);
        }
    }
};

class SmoothExactGsFunction : public gsFunction<real_t>
{
    SmoothManufacturedShellField m_field;
public:
    explicit SmoothExactGsFunction(const SmoothManufacturedShellField& f) : m_field(f) {}
    GISMO_CLONE_FUNCTION(SmoothExactGsFunction)

    short_t domainDim() const override { return 3; }
    short_t targetDim() const override { return 3; }
    const gsFunction<real_t>& piece(const index_t) const override { return *this; }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        const index_t n = u.cols();
        result.resize(3, n);
        gsVector<real_t> g(3);
        gsMatrix<real_t> H(3, 3);
        for (index_t k = 0; k < n; ++k)
        {
            const gsVector<real_t> xP = u.col(k);
            for (index_t c = 0; c < 3; ++c)
            {
                real_t v = 0;
                m_field(c, xP, v, g, H);
                result(c, k) = v;
            }
        }
    }
};

void exportExactAndSolutionPVD(const gsMultiPatch<real_t>& mpConst,
                               const gsMultiBasis<real_t>& dbasis,
                               const gsMultiPatch<real_t>& dispMPConst,
                               const SmoothManufacturedShellField& exactField,
                               const std::string& exportPrefix,
                               index_t nPts = 1000)
{
    if (exportPrefix.empty())
        return;

    const size_t slash = exportPrefix.find_last_of("/\\");
    const std::string parentDir = (slash != std::string::npos) ? exportPrefix.substr(0, slash + 1) : "./";
    mkdirRecursive(parentDir + "_vtk");

    gsMultiPatch<real_t> mp = mpConst;
    gsMultiPatch<real_t> dispMP = dispMPConst;
    SmoothExactGsFunction exactFunc(exactField);

    // 1. Export Exact Function PVD
    {
        gsExprEvaluator<real_t> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G = ev.getMap(mp);
        auto uEx = ev.getVariable(exactFunc, G);

        gsParaviewCollection col(exportPrefix + "_exact.pvd", &ev);
        col.options().setInt("numPoints", nPts);
        col.options().setInt("precision", 6);
        col.options().setSwitch("plotElements", true);
        col.options().setInt("plotElements.resolution", -1);
        col.options().setSwitch("plotControlNet", false);
        col.options().setSwitch("makeSubfolder", true);
        col.options().setString("subfolder", "_vtk");

        col.newTimeStep(&mp);
        col.addField(uEx, "ExactDisplacement");
        col.addField(uEx.norm(), "ExactMagnitude");
        col.saveTimeStep();
        col.save();
        gsInfo << "Exported ParaView exact field    -> " << exportPrefix << "_exact.pvd\n";
    }

    // 2. Export Numerical Solution (+ Error) PVD
    {
        gsExprEvaluator<real_t> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G = ev.getMap(mp);
        auto uEx = ev.getVariable(exactFunc, G);
        auto uSol = ev.getVariable(dispMP);

        gsParaviewCollection col(exportPrefix + "_solution.pvd", &ev);
        col.options().setInt("numPoints", nPts);
        col.options().setInt("precision", 6);
        col.options().setSwitch("plotElements", true);
        col.options().setInt("plotElements.resolution", -1);
        col.options().setSwitch("plotControlNet", false);
        col.options().setSwitch("makeSubfolder", true);
        col.options().setString("subfolder", "_vtk");

        col.newTimeStep(&mp);
        col.addField(uSol, "Displacement");
        col.addField(uSol.norm(), "DisplacementMagnitude");
        col.addField(uEx - uSol, "ErrorVector");
        col.addField((uEx - uSol).norm(), "PointwiseError");
        col.saveTimeStep();
        col.save();
        gsInfo << "Exported ParaView solution field -> " << exportPrefix << "_solution.pvd\n";
    }
}

gsBoundaryConditions<real_t> buildAllBoundaryBC(const gsMultiPatch<real_t>& mp,
                                                const std::string& bcMode)
{
    gsBoundaryConditions<real_t> bc;
    gsConstantFunction<real_t> zero(0.0, 3);
    std::string label = "Clamped";
    if (bcMode == "simply_supported")
        label = "Values";
    else if (bcMode == "strong_c1" || bcMode == "values_and_derivatives")
        label = "ValuesAndDerivatives";
    for (size_t i = 0; i < mp.boundaries().size(); ++i)
    {
        const patchSide& ps = mp.boundaries()[i];
        for (short_t comp = 0; comp < 3; ++comp)
            bc.add(ps.patch, ps.side(), label, zero, comp);
    }
    return bc;
}

void runManufacturedStudy(const gsMultiPatch<real_t>& mp,
                          const gsKLShellMaterial<real_t>& mat,
                          short_t degreeP,
                          short_t regularityK,
                          index_t minRef,
                          index_t maxRef,
                          const std::string& bcMode,
                          const std::string& solverType,
                          const std::string& smootherType,
                          real_t freq,
                          bool plateOnlyZ,
                          bool computeQI,
                          const std::string& exportPrefix)
{
    const short_t effK = (regularityK < 0) ? static_cast<short_t>(degreeP - 2) : regularityK;
    gsBoundaryConditions<real_t> bc = buildAllBoundaryBC(mp, bcMode);

    gsInfo << "\n========================================================================================================\n"
           << "AS-G1 Kirchhoff-Love Shell Manufactured Convergence Study (Hasanova, S. Takacs, T. Takacs 2026)\n"
           << "  Patches=" << mp.nPatches() << "  Interfaces=" << mp.nInterfaces()
           << "  Boundaries=" << mp.boundaries().size() << "  (open=" << (!mp.boundaries().empty() ? "YES" : "NO") << ")\n"
           << "  Space: S_{p=" << degreeP << ", k=" << effK << "} (knot mult=" << (degreeP - effK)
           << "), BC=" << bcMode << ", t=" << mat.thickness << ", E=" << mat.E << ", nu=" << mat.nu << "\n"
           << "  Solver=" << solverType << "  Smoother=" << smootherType << "\n"
           << "========================================================================================================\n";

    SmoothManufacturedShellField exactField;
    exactField.plateOnlyZ = plateOnlyZ;
    exactField.freq = freq;

    gsInfo << std::left
           << std::setw(4)  << "Ref"
           << std::setw(9)  << "h"
           << std::setw(8)  << "nFree"
           << std::setw(11) << "L2_rel"
           << std::setw(6)  << "EOC"
           << std::setw(11) << "H1_rel"
           << std::setw(6)  << "EOC"
           << std::setw(11) << "H2_rel"
           << std::setw(6)  << "EOC"
           << std::setw(11) << "En_rel"
           << std::setw(6)  << "EOC";
    if (computeQI)
        gsInfo << std::setw(11) << "QI_H2_rel" << std::setw(6) << "EOC";
    gsInfo << std::setw(6)  << "Iter"
           << std::setw(10) << "RelRes"
           << std::setw(8)  << "Time(s)"
           << "\n" << std::string(computeQI ? 130 : 113, '-') << "\n";

    real_t prevH = 0, prevL2 = 0, prevH1 = 0, prevH2 = 0, prevEn = 0, prevQIH2 = 0;

    for (index_t r = minRef; r <= maxRef; ++r)
    {
        gsAsG1ShellHierarchy<real_t> hier(mp, bc, degreeP, effK, minRef, r, /*normalizeGluing=*/true);
        const auto& fin = hier.finest();

        gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
        gsSparseMatrix<real_t> K_discont, K_free;
        gsVector<real_t> F_discont, F_free;
        asmbl.assembleStiffness(K_discont);
        asmbl.assembleManufacturedRHS(exactField, F_discont);

        gsVector<real_t> u_bnd, u_qi_free;
        projectVectorAsG1Exact3D(
            mp, fin.dbasis, hier.gluingData(), hier.cornerFrames(),
            fin.unscaledArgBasis, fin.colScales, fin.mapper, 3,
            exactField, u_bnd, computeQI ? &u_qi_free : nullptr);

        asmbl.projectSystem(K_discont, F_discont, fin.T_free, fin.T_bnd, u_bnd, K_free, F_free);

        gsVector<real_t> u_free;
        const auto st = hier.solve(K_free, F_free, mat, u_free, solverType, smootherType,
                                   /*numCycles=*/1, /*pre=*/1, /*post=*/1, /*maxIter=*/250, 1e-11);

        gsVector<real_t> c_sol = fin.T_free * u_free;
        if (fin.T_bnd.cols() > 0)
            c_sol += fin.T_bnd * u_bnd;
        const gsKLShellErrors<real_t> err = asmbl.computeErrors(c_sol, exactField);

        gsKLShellErrors<real_t> qiErr;
        if (computeQI)
        {
            gsVector<real_t> c_qi = fin.T_free * u_qi_free;
            if (fin.T_bnd.cols() > 0)
                c_qi += fin.T_bnd * u_bnd;
            qiErr = asmbl.computeErrors(c_qi, exactField);
        }

        auto eoc = [&](real_t eCur, real_t ePrev) -> real_t {
            if (r == minRef || eCur <= 0 || ePrev <= 0 || prevH <= 0)
                return 0.0;
            return std::log(ePrev / eCur) / std::log(prevH / fin.hLevel);
        };

        gsInfo << std::left
               << std::setw(4)  << r
               << std::fixed << std::setprecision(4) << std::setw(9) << fin.hLevel
               << std::setw(8)  << st.nFree
               << std::scientific << std::setprecision(2) << std::setw(11) << err.l2Rel()
               << std::fixed << std::setprecision(2) << std::setw(6) << eoc(err.l2Rel(), prevL2)
               << std::scientific << std::setprecision(2) << std::setw(11) << err.h1Rel()
               << std::fixed << std::setprecision(2) << std::setw(6) << eoc(err.h1Rel(), prevH1)
               << std::scientific << std::setprecision(2) << std::setw(11) << err.h2Rel()
               << std::fixed << std::setprecision(2) << std::setw(6) << eoc(err.h2Rel(), prevH2)
               << std::scientific << std::setprecision(2) << std::setw(11) << err.energyRel()
               << std::fixed << std::setprecision(2) << std::setw(6) << eoc(err.energyRel(), prevEn);
        if (computeQI)
        {
            gsInfo << std::scientific << std::setprecision(2) << std::setw(11) << qiErr.h2Rel()
                   << std::fixed << std::setprecision(2) << std::setw(6) << eoc(qiErr.h2Rel(), prevQIH2);
        }
        gsInfo << std::setw(6)  << st.iterations
               << std::scientific << std::setprecision(1) << std::setw(10) << st.relResidual
               << std::fixed << std::setprecision(3) << std::setw(8) << (st.setupTimeSec + st.solveTimeSec)
               << "\n";

        prevH = fin.hLevel;
        prevL2 = err.l2Rel();
        prevH1 = err.h1Rel();
        prevH2 = err.h2Rel();
        prevEn = err.energyRel();
        prevQIH2 = qiErr.h2Rel();

        if (r == maxRef && !exportPrefix.empty())
        {
            gsMultiPatch<real_t> dispMP = asmbl.buildDisplacementMultiPatch(c_sol);
            exportExactAndSolutionPVD(mp, fin.dbasis, dispMP, exactField, exportPrefix);
        }
    }
}

void runSolverComparison(const gsMultiPatch<real_t>& mp,
                         const gsKLShellMaterial<real_t>& mat,
                         short_t degreeP,
                         short_t regularityK,
                         index_t minRef,
                         index_t maxRef,
                         const std::string& bcMode,
                         real_t freq)
{
    const short_t effK = (regularityK < 0) ? static_cast<short_t>(degreeP - 2) : regularityK;
    gsBoundaryConditions<real_t> bc = buildAllBoundaryBC(mp, bcMode);
    gsAsG1ShellHierarchy<real_t> hier(mp, bc, degreeP, effK, minRef, maxRef, true);
    const auto& fin = hier.finest();

    SmoothManufacturedShellField exactField;
    exactField.freq = freq;

    gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
    gsSparseMatrix<real_t> K_discont, K_free;
    gsVector<real_t> F_discont, F_free, u_bnd;
    asmbl.assembleStiffness(K_discont);
    asmbl.assembleManufacturedRHS(exactField, F_discont);
    projectVectorAsG1Exact3D(
        mp, fin.dbasis, hier.gluingData(), hier.cornerFrames(),
        fin.unscaledArgBasis, fin.colScales, fin.mapper, 3, exactField, u_bnd,
        static_cast<gsVector<real_t>*>(nullptr));
    asmbl.projectSystem(K_discont, F_discont, fin.T_free, fin.T_bnd, u_bnd, K_free, F_free);

    gsInfo << "\n==========================================================================================\n"
           << "Solver & Smoother Benchmark on Open 3D AS-G1 Shell (p=" << degreeP << ", k=" << effK
           << ", levels=" << minRef << ".." << maxRef << ", nFree=" << fin.mapper.freeSize()
           << ", t=" << mat.thickness << ")\n"
           << "==========================================================================================\n";
    gsInfo << std::left
           << std::setw(14) << "Solver"
           << std::setw(14) << "Smoother"
           << std::setw(8)  << "Iters"
           << std::setw(12) << "RelRes"
           << std::setw(11) << "Setup(s)"
           << std::setw(11) << "Solve(s)"
           << std::setw(12) << "H2_rel"
           << "\n" << std::string(82, '-') << "\n";

    const std::vector<std::pair<std::string, std::string>> configs = {
        {"direct",     "-"},
        {"schwarz_cg", "-"},
        {"cg",         "nodal_gs"},
        {"cg",         "ring1"},
        {"cg",         "line_ring1"},
        {"cg",         "scms"},
        {"cg",         "scms_hybrid"},
        {"mg",         "line_ring1"}
    };

    for (const auto& cfg : configs)
    {
        gsVector<real_t> u_free;
        const auto st = hier.solve(K_free, F_free, mat, u_free, cfg.first, cfg.second,
                                   /*numCycles=*/1, /*pre=*/1, /*post=*/1, /*maxIter=*/250, 1e-10);
        gsVector<real_t> c_sol = fin.T_free * u_free;
        if (fin.T_bnd.cols() > 0) c_sol += fin.T_bnd * u_bnd;
        const auto err = asmbl.computeErrors(c_sol, exactField);

        gsInfo << std::left
               << std::setw(14) << cfg.first
               << std::setw(14) << cfg.second
               << std::setw(8)  << st.iterations
               << std::scientific << std::setprecision(2) << std::setw(12) << st.relResidual
               << std::fixed << std::setprecision(4) << std::setw(11) << st.setupTimeSec
               << std::fixed << std::setprecision(4) << std::setw(11) << st.solveTimeSec
               << std::scientific << std::setprecision(2) << std::setw(12) << err.h2Rel()
               << "\n";
    }
}

void runScordelisLoBenchmark(const std::string& file,
                             short_t degreeP,
                             short_t regularityK,
                             index_t minRef,
                             index_t maxRef,
                             const std::string& solverType,
                             const std::string& smootherType)
{
    gsMultiPatch<real_t> mp = loadOpenSurface3D(file);
    gsKLShellMaterial<real_t> mat(4.32e8, 0.0, 0.25);

    // Compute overall geometry bounding box for scale-invariant boundary classification
    gsVector<real_t, 3> bbMin = gsVector<real_t, 3>::Constant(1e30);
    gsVector<real_t, 3> bbMax = gsVector<real_t, 3>::Constant(-1e30);
    {
        gsMatrix<real_t> grid(2, 25);
        for (int j = 0; j < 5; ++j)
            for (int i = 0; i < 5; ++i)
            {
                grid(0, i + 5 * j) = i * 0.25;
                grid(1, i + 5 * j) = j * 0.25;
            }
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsMatrix<real_t> X = mp.patch(p).eval(grid);
            for (int c = 0; c < 3; ++c)
            {
                bbMin(c) = std::min(bbMin(c), X.row(c).minCoeff());
                bbMax(c) = std::max(bbMax(c), X.row(c).maxCoeff());
            }
        }
    }
    const gsVector<real_t, 3> bbExt = bbMax - bbMin;
    const real_t scale = std::max(bbExt.maxCoeff(), real_t(1e-12));

    // Identify curved diaphragm boundaries (constant axial coordinate, varying Z) vs straight free edges.
    // On the exact MacNeal-Harder / Kiendl Scordelis-Lo benchmark, transverse displacements are zero
    // on both curved ends, while axial warping is free (only the 1-DOF rigid axial shift is removed).
    gsBoundaryConditions<real_t> bc;
    gsConstantFunction<real_t> zero(0.0, 3);
    std::vector<patchSide> freeSides;
    std::vector<patchSide> diaphragmSides;
    index_t axialComp = 0;
    for (size_t i = 0; i < mp.boundaries().size(); ++i)
    {
        const patchSide& ps = mp.boundaries()[i];
        const gsGeometry<real_t>& g = mp.patch(ps.patch);
        gsMatrix<real_t> pts(2, 3);
        const short_t nDir = ps.direction(), tDir = 1 - nDir;
        pts.row(nDir).setConstant(ps.parameter() ? 1.0 : 0.0);
        pts(tDir, 0) = 0.0; pts(tDir, 1) = 0.5; pts(tDir, 2) = 1.0;
        gsMatrix<real_t> X = g.eval(pts);
        const real_t dX = X.row(0).maxCoeff() - X.row(0).minCoeff();
        const real_t dY = X.row(1).maxCoeff() - X.row(1).minCoeff();
        const real_t dZ = X.row(2).maxCoeff() - X.row(2).minCoeff();
        if (dZ > 0.05 * bbExt(2) && (dY < 1e-3 * scale || dX < 1e-3 * scale))
        {
            diaphragmSides.push_back(ps);
            if (dY < 1e-3 * scale)
            {
                axialComp = 1;
                bc.add(ps.patch, ps.side(), "Values", zero, 0); // u_x = 0
                bc.add(ps.patch, ps.side(), "Values", zero, 2); // u_z = 0
            }
            else
            {
                axialComp = 0;
                bc.add(ps.patch, ps.side(), "Values", zero, 1); // u_y = 0
                bc.add(ps.patch, ps.side(), "Values", zero, 2); // u_z = 0
            }
        }
        else
        {
            freeSides.push_back(ps);
        }
    }
    GISMO_ENSURE(!diaphragmSides.empty(),
                 "runScordelisLoBenchmark: failed to identify curved diaphragm boundaries.");

    gsInfo << "\n================================================================================\n"
           << "Scordelis-Lo Roof Benchmark (Non-Isoparametric Polynomial S_{p=" << degreeP
           << "} on Exact NURBS)\n"
           << "  Reference free-edge midpoint vertical displacement u_z ~ -0.3006 (KL) .. -0.3024\n"
           << "================================================================================\n";

    auto gravityLoad = [](const gsVector<real_t>&, gsVector<real_t, 3>& fBody, real_t& pNorm) {
        fBody << 0.0, 0.0, -90.0;
        pNorm = 0.0;
    };

    const real_t midAx = 0.5 * (bbMin(axialComp) + bbMax(axialComp));

    for (index_t r = minRef; r <= maxRef; ++r)
    {
        gsAsG1ShellHierarchy<real_t> hier(mp, bc, degreeP, regularityK, minRef, r, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
        gsSparseMatrix<real_t> K_discont, K_free;
        gsVector<real_t> F_discont, F_free, u_bnd, u_free;
        u_bnd.setZero(fin.mapper.boundarySize());
        asmbl.assembleStiffness(K_discont);
        asmbl.assemblePhysicalRHS(gravityLoad, F_discont);
        asmbl.projectSystem(K_discont, F_discont, fin.T_free, fin.T_bnd, u_bnd, K_free, F_free);

        // Pin the single axial rigid-body translation mode on a guaranteed free corner value DOF
        bool axialPinned = false;
        for (size_t p = 0; p < mp.nPatches() && !axialPinned; ++p)
        {
            for (index_t c = 0; c < 4 && !axialPinned; ++c)
            {
                const index_t cornerOff = sumUntil(fin.argBasis[p].sizes, 9 + c);
                if (fin.mapper.is_free(cornerOff, static_cast<index_t>(p), axialComp))
                {
                    const index_t pinIdx = fin.mapper.index(cornerOff, static_cast<index_t>(p), axialComp);
                    K_free.coeffRef(pinIdx, pinIdx) += mat.E;
                    axialPinned = true;
                }
            }
        }
        GISMO_ENSURE(axialPinned,
                     "runScordelisLoBenchmark: could not find a free axial corner DOF to pin rigid translation.");

        const auto st = hier.solve(K_free, F_free, mat, u_free, solverType, smootherType,
                                   1, 1, 1, 250, 1e-10);
        gsVector<real_t> c_sol = fin.T_free * u_free;

        // Evaluate vertical displacement u_z at the physical axial midpoint along straight free edge(s)
        // (works identically for both single-patch and multi-patch split roofs).
        real_t minUz = 0.0;
        bool foundMid = false;
        for (const patchSide& ps : freeSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            const real_t samples[3] = {0.0, 0.5, 1.0};
            for (int s = 0; s < 3; ++s)
            {
                gsMatrix<real_t> pt(2, 1);
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = samples[s];
                const gsMatrix<real_t> X = mp.patch(ps.patch).eval(pt);
                if (std::abs(X(axialComp, 0) - midAx) <= 1e-3 * scale)
                {
                    const auto uVal = asmbl.evalDisplacement(c_sol, ps.patch, pt(0, 0), pt(1, 0));
                    if (!foundMid || uVal(2) < minUz)
                    {
                        minUz = uVal(2);
                        foundMid = true;
                    }
                }
            }
        }
        GISMO_ENSURE(foundMid,
                     "runScordelisLoBenchmark: failed to locate free-edge axial midpoint.");

        // Check diaphragm Dirichlet constraint residual across all diaphragm edges including corners
        real_t maxDiaphLeak = 0.0;
        const index_t transComp = 1 - axialComp;
        for (const patchSide& ps : diaphragmSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            for (int s = 0; s <= 20; ++s)
            {
                gsMatrix<real_t> pt(2, 1);
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = s / 20.0;
                const auto uVal = asmbl.evalDisplacement(c_sol, ps.patch, pt(0, 0), pt(1, 0));
                maxDiaphLeak = std::max(maxDiaphLeak, std::max(std::abs(uVal(transComp)), std::abs(uVal(2))));
            }
        }

        gsInfo << "  ref=" << r << "  nFree=" << std::setw(6) << st.nFree
               << "  u_z(free_mid)=" << std::fixed << std::setprecision(6) << minUz
               << "  max|u_diaph|=" << std::scientific << std::setprecision(2) << maxDiaphLeak
               << "  iters=" << std::setw(4) << st.iterations
               << "  relRes=" << std::scientific << std::setprecision(1) << st.relResidual << "\n";
    }
}

void runPhysicsVerification(const std::string& scordelisFile)
{
    gsInfo << "\n================================================================================\n"
           << "AS-G1 Kirchhoff-Love Shell Intrinsic & Physical Regression Suite\n"
           << "================================================================================\n";

    auto embed2DTo3D = [](const gsMultiPatch<real_t>& mp2) -> gsMultiPatch<real_t> {
        gsMultiPatch<real_t> mp3;
        for (size_t p = 0; p < mp2.nPatches(); ++p)
        {
            const auto& bsp = dynamic_cast<const gsTensorBSpline<2, real_t>&>(mp2.patch(p));
            gsMatrix<real_t> c3 = gsMatrix<real_t>::Zero(bsp.coefs().rows(), 3);
            c3.leftCols(2) = bsp.coefs();
            mp3.addPatch(typename gsGeometry<real_t>::uPtr(new gsTensorBSpline<2, real_t>(bsp.basis(), c3)));
        }
        mp3.computeTopology();
        return mp3;
    };

    // 1. Intrinsic checks (symmetry, 6 rigid-body modes, manufactured RHS consistency) on curved 2x2 paraboloid
    {
        gsMultiPatch<real_t> mp = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            auto& bsp = dynamic_cast<gsTensorBSpline<2, real_t>&>(mp.patch(p));
            bsp.degreeElevate(1); // quadratic in u, v so z = 0.25*(x^2 + y^2) is exact
            const gsMatrix<real_t> anc = bsp.basis().anchors();
            const gsMatrix<real_t> X = bsp.eval(anc);
            gsMatrix<real_t> vals(3, anc.cols());
            for (index_t j = 0; j < anc.cols(); ++j)
            {
                vals(0, j) = X(0, j);
                vals(1, j) = X(1, j);
                vals(2, j) = 0.25 * (X(0, j) * X(0, j) + X(1, j) * X(1, j));
            }
            bsp.coefs() = bsp.basis().interpolateAtAnchors(vals)->coefs();
        }
        mp.computeTopology();

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.1);
        gsMultiBasis<real_t> db = buildAsG1PolynomialMultiBasis(mp, 3, 1, 0);
        gsAsG1KLShellAssembler<real_t> asmbl(mp, db, mat);
        gsSparseMatrix<real_t> K;
        asmbl.assembleStiffness(K);
        const real_t kNorm = K.norm();
        const real_t symErr = (gsSparseMatrix<real_t>(K - gsSparseMatrix<real_t>(K.transpose()))).norm() / kNorm;
        GISMO_ENSURE(symErr < 1e-14, "Stiffness matrix asymmetry too large: " << symErr);

        const index_t N = asmbl.scalarDiscontSize();
        real_t maxRbErr = 0.0;
        for (int m = 0; m < 6; ++m)
        {
            gsVector<real_t> c(3 * N);
            for (size_t p = 0; p < mp.nPatches(); ++p)
            {
                const gsMatrix<real_t> anc = db[p].anchors();
                const gsMatrix<real_t> X = mp.patch(p).eval(anc);
                gsMatrix<real_t> vals(3, anc.cols());
                for (index_t j = 0; j < anc.cols(); ++j)
                {
                    gsVector<real_t, 3> u = gsVector<real_t, 3>::Zero();
                    if (m < 3)
                        u(m) = 1.0;
                    else
                    {
                        gsVector<real_t, 3> ax = gsVector<real_t, 3>::Zero();
                        ax(m - 3) = 1.0;
                        gsVector<real_t, 3> x3;
                        x3 << X(0, j), X(1, j), X(2, j);
                        u = ax.cross(x3);
                    }
                    vals.col(j) = u;
                }
                auto g = db[p].interpolateAtAnchors(vals);
                for (index_t comp = 0; comp < 3; ++comp)
                    c.segment(comp * N + asmbl.patchOffsets()[p], db[p].size()) = g->coefs().col(comp);
            }
            const real_t rbErr = (K * c).norm() / (kNorm * c.norm());
            maxRbErr = std::max(maxRbErr, rbErr);
        }
        GISMO_ENSURE(maxRbErr < 1e-13, "Rigid-body mode kernel residual too large: " << maxRbErr);
        gsInfo << "  [PASS] Intrinsic curved 2x2 paraboloid: sym=" << std::scientific << std::setprecision(2)
               << symErr << ", maxRigidBody=" << maxRbErr << "\n";
    }

    // 2. Clamped and Simply-Supported 2x2 Kirchhoff plates (verifies B2 smooth-boundary clamping fix)
    {
        gsMultiPatch<real_t> grid = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.01);
        const real_t D = mat.E * std::pow(mat.thickness, 3) / (12.0 * (1.0 - mat.nu * mat.nu));

        const std::pair<std::string, real_t> cases[3] = {
            {"ValuesAndDerivatives", 0.00126532 / D},
            {"Clamped",              0.00126532 / D},
            {"Values",               0.00406235 / D}
        };
        for (const auto& cs : cases)
        {
            gsBoundaryConditions<real_t> bc;
            gsConstantFunction<real_t> zero(0.0, 3);
            for (const auto& ps : grid.boundaries())
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), cs.first, zero, c);

            gsAsG1ShellHierarchy<real_t> hier(grid, bc, 3, 1, 0, 1, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(grid, fin.dbasis, mat, bc);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-12);
            gsVector<real_t> c_sol = fin.T_free * uf;

            real_t wc = 0.0;
            for (size_t p = 0; p < grid.nPatches(); ++p)
                for (int k = 0; k < 4; ++k)
                {
                    gsMatrix<real_t> u(2, 1);
                    u << (k % 2), (k / 2);
                    const gsMatrix<real_t> X = grid.patch(p).eval(u);
                    if (std::abs(X(0, 0) - 0.5) < 1e-10 && std::abs(X(1, 0) - 0.5) < 1e-10)
                        wc = asmbl.evalDisplacement(c_sol, p, u(0, 0), u(1, 0))(2);
                }
            const real_t ratio = wc / cs.second;
            GISMO_ENSURE(st.converged && std::abs(ratio - 1.0) < 1e-4,
                         "2x2 plate (" << cs.first << ") failed: w_c/w_ref = " << ratio);
            gsInfo << "  [PASS] 2x2 plate (" << std::left << std::setw(20) << cs.first
                   << "): w_c/w_ref=" << std::fixed << std::setprecision(6) << ratio
                   << ", iters=" << st.iterations << "\n";
        }
    }

    // 3. Multi-patch 2x2 Cantilever plate in both edge orientations (nu = 0 -> exact quartic beam w_tip = qL^4/(8D))
    {
        gsMultiPatch<real_t> grid = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        gsKLShellMaterial<real_t> mat(1e6, 0.0, 0.01);
        const real_t D = mat.E * std::pow(mat.thickness, 3) / 12.0;
        const real_t wTipRef = 1.0 / (8.0 * D);

        for (int clampAxis = 0; clampAxis < 2; ++clampAxis)
        {
            gsBoundaryConditions<real_t> bc;
            gsConstantFunction<real_t> zero(0.0, 3);
            for (const auto& ps : grid.boundaries())
            {
                gsMatrix<real_t> pt(2, 1);
                const short_t nDir = ps.direction(), tDir = 1 - nDir;
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = 0.5;
                const gsMatrix<real_t> X = grid.patch(ps.patch).eval(pt);
                if (std::abs(X(clampAxis, 0)) < 1e-10)
                    for (short_t c = 0; c < 3; ++c)
                        bc.add(ps.patch, ps.side(), "ValuesAndDerivatives", zero, c);
            }

            gsAsG1ShellHierarchy<real_t> hier(grid, bc, 4, 2, 0, 1, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(grid, fin.dbasis, mat);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "direct");
            gsVector<real_t> c_sol = fin.T_free * uf;

            real_t maxTipErr = 0.0;
            real_t maxClampLeak = 0.0;
            for (size_t p = 0; p < grid.nPatches(); ++p)
            {
                for (int j = 0; j <= 4; ++j)
                    for (int i = 0; i <= 4; ++i)
                    {
                        gsMatrix<real_t> pt(2, 1);
                        pt << i * 0.25, j * 0.25;
                        const gsMatrix<real_t> X = grid.patch(p).eval(pt);
                        const auto uVal = asmbl.evalDisplacement(c_sol, p, pt(0, 0), pt(1, 0));
                        if (std::abs(X(clampAxis, 0) - 1.0) < 1e-10)
                            maxTipErr = std::max(maxTipErr, std::abs(uVal(2) - wTipRef) / wTipRef);
                        if (std::abs(X(clampAxis, 0)) < 1e-10)
                            maxClampLeak = std::max(maxClampLeak, uVal.norm());
                    }
            }
            GISMO_ENSURE(st.converged && maxTipErr < 1e-8 && maxClampLeak < 1e-12,
                         "2x2 cantilever (axis=" << clampAxis << ") failed: tipRelErr="
                         << maxTipErr << ", clampLeak=" << maxClampLeak);
            gsInfo << "  [PASS] 2x2 cantilever (clamp axis=" << clampAxis
                   << "): tipRelErr=" << std::scientific << std::setprecision(2) << maxTipErr
                   << ", clampLeak=" << maxClampLeak << "\n";
        }
    }

    // 4. Curved sharp-corner one-sided Dirichlet test on quarter-annulus (r in [1,2], theta in [0,pi/2])
    {
        gsMultiPatch<real_t> ann;
        ann.addPatch(gsNurbsCreator<real_t>::NurbsQuarterAnnulus(1.0, 2.0));
        // Embed 2D NURBS in 3D (z = 0)
        const auto& r2 = dynamic_cast<const gsTensorNurbs<2, real_t>&>(ann.patch(0));
        gsMatrix<real_t> c3 = gsMatrix<real_t>::Zero(r2.coefs().rows(), 3);
        c3.leftCols(2) = r2.coefs();
        gsMultiPatch<real_t> ann3;
        ann3.addPatch(typename gsGeometry<real_t>::uPtr(new gsTensorNurbs<2, real_t>(r2.basis(), c3)));
        ann3.computeTopology();

        // Constrain Values ONLY on inner curved arc (r = 1) and outer curved arc (r = 2), leaving radial edges free
        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        std::vector<patchSide> curvedSides;
        for (const auto& ps : ann3.boundaries())
        {
            gsMatrix<real_t> pts(2, 3);
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            pts.row(nDir).setConstant(ps.parameter() ? 1.0 : 0.0);
            pts(tDir, 0) = 0.0; pts(tDir, 1) = 0.5; pts(tDir, 2) = 1.0;
            const gsMatrix<real_t> X = ann3.patch(0).eval(pts);
            const real_t r0 = std::hypot(X(0, 0), X(1, 0));
            const real_t r1 = std::hypot(X(0, 1), X(1, 1));
            const real_t r2v = std::hypot(X(0, 2), X(1, 2));
            if (std::abs(r0 - r1) < 1e-4 && std::abs(r1 - r2v) < 1e-4)
            {
                curvedSides.push_back(ps);
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), "Values", zero, c);
            }
        }
        GISMO_ENSURE(curvedSides.size() == 2, "Expected 2 curved arcs on quarter annulus.");

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
        gsAsG1ShellHierarchy<real_t> hier(ann3, bc, 4, 2, 0, 1, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(ann3, fin.dbasis, mat);
        gsSparseMatrix<real_t> K, Kf;
        gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
        asmbl.assembleStiffness(K);
        asmbl.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
        asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
        const auto st = hier.solve(Kf, Ff, mat, uf, "direct");
        gsVector<real_t> c_sol = fin.T_free * uf;

        real_t maxArcLeak = 0.0;
        for (const patchSide& ps : curvedSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            for (int s = 0; s <= 20; ++s)
            {
                gsMatrix<real_t> pt(2, 1);
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = s / 20.0;
                maxArcLeak = std::max(maxArcLeak, asmbl.evalDisplacement(c_sol, 0, pt(0, 0), pt(1, 0)).norm());
            }
        }
        GISMO_ENSURE(st.converged && maxArcLeak < 1e-11,
                     "Curved sharp-corner quarter-annulus Dirichlet leak too large: " << maxArcLeak);
        gsInfo << "  [PASS] Quarter-annulus curved sharp-corner BCs: maxArcLeak="
               << std::scientific << std::setprecision(2) << maxArcLeak << "\n";
    }

    // 5. Scordelis-Lo cylindrical NURBS roof (verifies B1 sharp-corner one-sided diaphragm fix)
    {
        gsMultiPatch<real_t> mp = loadOpenSurface3D(scordelisFile);
        gsKLShellMaterial<real_t> mat(4.32e8, 0.0, 0.25);
        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        std::vector<patchSide> diaphragmSides, freeSides;
        index_t axialComp = 0;
        for (const auto& ps : mp.boundaries())
        {
            gsMatrix<real_t> pts(2, 3);
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            pts.row(nDir).setConstant(ps.parameter() ? 1.0 : 0.0);
            pts(tDir, 0) = 0.0; pts(tDir, 1) = 0.5; pts(tDir, 2) = 1.0;
            const gsMatrix<real_t> X = mp.patch(ps.patch).eval(pts);
            const real_t dX = X.row(0).maxCoeff() - X.row(0).minCoeff();
            const real_t dY = X.row(1).maxCoeff() - X.row(1).minCoeff();
            const real_t dZ = X.row(2).maxCoeff() - X.row(2).minCoeff();
            if (dZ > 1.0 && (dY < 1e-3 || dX < 1e-3))
            {
                diaphragmSides.push_back(ps);
                const index_t tComp = (dY < 1e-3) ? 0 : 1;
                axialComp = 1 - tComp;
                bc.add(ps.patch, ps.side(), "Values", zero, tComp);
                bc.add(ps.patch, ps.side(), "Values", zero, 2);
            }
            else
            {
                freeSides.push_back(ps);
            }
        }
        gsAsG1ShellHierarchy<real_t> hier(mp, bc, 4, 2, 0, 2, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
        gsSparseMatrix<real_t> K, Kf;
        gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
        asmbl.assembleStiffness(K);
        asmbl.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f << 0.0, 0.0, -90.0; pn = 0.0; }, F);
        asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
        const index_t c0Off = sumUntil(fin.argBasis[0].sizes, 9);
        GISMO_ENSURE(fin.mapper.is_free(c0Off, 0, axialComp), "Scordelis axial corner DOF unexpectedly constrained.");
        const index_t pinIdx = fin.mapper.index(c0Off, 0, axialComp);
        Kf.coeffRef(pinIdx, pinIdx) += mat.E;

        const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 250, 1e-11);
        gsVector<real_t> c_sol = fin.T_free * uf;

        real_t minUz = 0.0;
        for (const patchSide& ps : freeSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            gsMatrix<real_t> pt(2, 1);
            pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
            pt(tDir, 0) = 0.5;
            minUz = std::min(minUz, asmbl.evalDisplacement(c_sol, ps.patch, pt(0, 0), pt(1, 0))(2));
        }
        real_t maxDiaphLeak = 0.0;
        const index_t transComp = 1 - axialComp;
        for (const patchSide& ps : diaphragmSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            for (int s = 0; s <= 20; ++s)
            {
                gsMatrix<real_t> pt(2, 1);
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = s / 20.0;
                const auto uVal = asmbl.evalDisplacement(c_sol, ps.patch, pt(0, 0), pt(1, 0));
                maxDiaphLeak = std::max(maxDiaphLeak, std::max(std::abs(uVal(transComp)), std::abs(uVal(2))));
            }
        }
        GISMO_ENSURE(st.converged && std::abs(minUz + 0.3006) < 1e-3 && maxDiaphLeak < 1e-11,
                     "Scordelis-Lo verification failed: u_z=" << minUz << ", diaphLeak=" << maxDiaphLeak);
        gsInfo << "  [PASS] Scordelis-Lo NURBS roof (p=4, ref=2): u_z=" << std::fixed << std::setprecision(5)
               << minUz << " (ref -0.3006), max|u_diaph|=" << std::scientific << std::setprecision(2)
               << maxDiaphLeak << ", iters=" << st.iterations << "\n";
    }

    // 6. Curved 2x2 paraboloid physically Clamped shell (Symmetric Nitsche theta_nu = 0, nonzero membrane strain eps_nunu != 0)
    {
        gsMultiPatch<real_t> mp = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        const real_t alpha = 0.25;
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            auto& bsp = dynamic_cast<gsTensorBSpline<2, real_t>&>(mp.patch(p));
            bsp.degreeElevate(1);
            const gsMatrix<real_t> anc = bsp.basis().anchors();
            const gsMatrix<real_t> X = bsp.eval(anc);
            gsMatrix<real_t> vals(3, anc.cols());
            for (index_t j = 0; j < anc.cols(); ++j)
            {
                vals(0, j) = X(0, j);
                vals(1, j) = X(1, j);
                vals(2, j) = alpha * (X(0, j) * X(0, j) + X(1, j) * X(1, j));
            }
            bsp.coefs() = bsp.basis().interpolateAtAnchors(vals)->coefs();
        }
        mp.computeTopology();

        // Exact physical clamped field on z = alpha*(x^2 + y^2):
        //   u = (phi, phi, 2*alpha*(x+y)*phi + psi), phi = sin(pi x)sin(pi y), psi = sin^2(pi x)sin^2(pi y).
        // Satisfies u|_Gamma = 0 AND A_3 . d_nu u|_Gamma = 0 identically, while eps_nunu|_Gamma != 0!
        auto exactClampedCurved = [&](index_t comp, const gsVector<real_t>& X,
                                      real_t& val, gsVector<real_t>& g, gsMatrix<real_t>& H)
        {
            const real_t pi = EIGEN_PI;
            const real_t x = X(0), y = X(1);
            const real_t sx = std::sin(pi * x), cx = std::cos(pi * x);
            const real_t sy = std::sin(pi * y), cy = std::cos(pi * y);
            const real_t phi = sx * sy;
            const real_t px = pi * cx * sy, py = pi * sx * cy;
            const real_t pxx = -pi * pi * phi, pyy = -pi * pi * phi, pxy = pi * pi * cx * cy;

            g.setZero(3);
            H.setZero(3, 3);
            if (comp == 0 || comp == 1)
            {
                val = phi;
                g(0) = px; g(1) = py;
                H(0, 0) = pxx; H(1, 1) = pyy; H(0, 1) = H(1, 0) = pxy;
            }
            else
            {
                const real_t c2 = 2.0 * alpha;
                const real_t gx = sx * sx, gy = sy * sy;
                const real_t dgx = pi * std::sin(2.0 * pi * x), dgy = pi * std::sin(2.0 * pi * y);
                const real_t ddgx = 2.0 * pi * pi * std::cos(2.0 * pi * x);
                const real_t ddgy = 2.0 * pi * pi * std::cos(2.0 * pi * y);
                const real_t psi = gx * gy;
                const real_t wx = dgx * gy, wy = gx * dgy;
                const real_t wxx = ddgx * gy, wyy = gx * ddgy, wxy = dgx * dgy;

                val = c2 * (x + y) * phi + psi;
                g(0) = c2 * phi + c2 * (x + y) * px + wx;
                g(1) = c2 * phi + c2 * (x + y) * py + wy;
                H(0, 0) = 2.0 * c2 * px + c2 * (x + y) * pxx + wxx;
                H(1, 1) = 2.0 * c2 * py + c2 * (x + y) * pyy + wyy;
                H(0, 1) = H(1, 0) = c2 * (px + py) + c2 * (x + y) * pxy + wxy;
            }
        };

        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        for (const auto& ps : mp.boundaries())
            for (short_t c = 0; c < 3; ++c)
                bc.add(ps.patch, ps.side(), "Clamped", zero, c);

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
        real_t prevL2 = 0.0, prevH2 = 0.0, prevRotErr = 0.0;
        real_t lastL2 = 0.0, lastH2 = 0.0, lastRotErr = 0.0;
        index_t lastIters = 0;
        for (index_t r = 0; r <= 1; ++r)
        {
            gsAsG1ShellHierarchy<real_t> hier(mp, bc, 3, 1, 0, r, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub, uf;
            asmbl.assembleStiffness(K);
            asmbl.assembleManufacturedRHS(exactClampedCurved, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-12);
            GISMO_ENSURE(st.converged, "Curved Clamped Nitsche MG-PCG failed to converge at ref=" << r);
            gsVector<real_t> c_sol = fin.T_free * uf;
            const auto err = asmbl.computeErrors(c_sol, exactClampedCurved);
            if (r == 0)
            {
                prevL2 = err.l2Rel();
                prevH2 = err.h2Rel();
                prevRotErr = asmbl.computeRotationBCErrorL2(c_sol);
            }
            else
            {
                lastL2 = err.l2Rel();
                lastH2 = err.h2Rel();
                lastRotErr = asmbl.computeRotationBCErrorL2(c_sol);
                lastIters = st.iterations;
            }
        }
        const real_t eocL2 = std::log2(prevL2 / lastL2);
        const real_t eocH2 = std::log2(prevH2 / lastH2);
        const real_t eocRot = std::log2(prevRotErr / lastRotErr);
        GISMO_ENSURE(eocL2 > 3.5 && eocH2 > 1.7 && eocRot > 2.5 && lastH2 < 5e-2 && lastRotErr < 1e-3,
                     "Curved Clamped Nitsche check failed: eocL2=" << eocL2 << ", eocH2=" << eocH2
                     << ", eocRot=" << eocRot << ", H2=" << lastH2 << ", rotErr=" << lastRotErr);
        gsInfo << "  [PASS] Curved 2x2 Clamped Nitsche (eps_nunu!=0): L2=" << std::scientific << std::setprecision(2)
               << lastL2 << " (EOC=" << std::fixed << std::setprecision(2) << eocL2 << "), H2="
               << std::scientific << lastH2 << " (EOC=" << std::fixed << eocH2
               << "), rotL2=" << std::scientific << lastRotErr << " (EOC=" << std::fixed << eocRot
               << "), iters=" << lastIters << "\n";
    }

    // 7. Inhomogeneous NormalRotation BC patch test (exact tilted plane u_z(x,y) = slope * x under zero load)
    {
        gsMultiPatch<real_t> sq;
        sq.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0));
        sq = embed2DTo3D(sq);

        const real_t slope = 0.05;
        // On West side (x = 0), outward conormal nu = (-1, 0, 0), A3 = (0, 0, 1) -> theta_nu = -du_z/dx = -slope
        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        gsConstantFunction<real_t> rotVal(-slope, 3);
        bc.add(0, boundary::west, "Values", zero, 0);
        bc.add(0, boundary::west, "Values", zero, 1);
        bc.add(0, boundary::west, "Values", zero, 2);
        bc.add(0, boundary::west, "NormalRotation", rotVal, 0);

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
        gsAsG1ShellHierarchy<real_t> hier(sq, bc, 3, 1, 0, 0, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(sq, fin.dbasis, mat, bc);
        gsSparseMatrix<real_t> K, Kf;
        gsVector<real_t> F, Ff, ub, uf;
        asmbl.assembleStiffness(K);
        asmbl.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 0.0; }, F);
        asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
        const auto st = hier.solve(Kf, Ff, mat, uf, "direct");
        gsVector<real_t> c_sol = fin.T_free * uf;

        real_t maxPlaneErr = 0.0;
        for (int j = 0; j <= 4; ++j)
            for (int i = 0; i <= 4; ++i)
            {
                const real_t x = i * 0.25, y = j * 0.25;
                const auto u = asmbl.evalDisplacement(c_sol, 0, x, y);
                gsVector<real_t, 3> uRef; uRef << 0.0, 0.0, slope * x;
                maxPlaneErr = std::max(maxPlaneErr, (u - uRef).norm());
            }
        const real_t rotErr = asmbl.computeRotationBCErrorL2(c_sol);
        GISMO_ENSURE(st.converged && maxPlaneErr < 1e-12 && rotErr < 1e-12,
                     "Inhomogeneous NormalRotation patch test failed: planeErr=" << maxPlaneErr
                     << ", rotErr=" << rotErr);
        gsInfo << "  [PASS] Inhomogeneous NormalRotation patch test: maxPlaneErr="
               << std::scientific << std::setprecision(2) << maxPlaneErr << ", rotL2=" << rotErr << "\n";
    }

    gsInfo << "All intrinsic & physical verification checks PASSED.\n";
}

} // anonymous namespace

int main(int argc, char* argv[])
{
    std::string geomFile = "asg1_reparametrized_domains/farahat_kl_shells/xml/surfHalfSphere_asg1.xml";
    std::string testMode = "manufactured";
    std::string bcMode = "clamped";
    std::string solverType = "cg";
    std::string smootherType = "line_ring1";
    std::string exportPrefix = "";
    int degreeP = 3;
    int regularityK = -1; // -1 => maximal AS-G1 smoothness k = p - 2
    int minRef = 0;
    int maxRef = 3;
    real_t E = 1e6;
    real_t nu = 0.3;
    real_t thickness = 0.1;
    real_t freq = 1.0;
    bool plateOnlyZ = false;
    bool skipQI = false;

    gsCmdLine cmd("AS-G1 Multi-Patch & Multigrid Kirchhoff-Love Shell Solver (Hasanova-Takacs-Takacs 2026).");
    cmd.addString("f", "file", "Input open 3D surface or 2D domain XML file.", geomFile);
    cmd.addString("t", "test", "Mode: manufactured | benchmark_solvers | p_sweep | scordelis_lo | verify_physics.", testMode);
    cmd.addString("b", "bc", "Boundary condition: clamped | simply_supported.", bcMode);
    cmd.addString("s", "solver", "Solver: direct | schwarz_cg | mg | cg (MG-PCG).", solverType);
    cmd.addString("m", "smoother", "MG smoother: nodal_gs | ring1 | line_ring1 | scms | scms_hybrid.", smootherType);
    cmd.addString("o", "export", "Optional ParaView output prefix.", exportPrefix);
    cmd.addInt("p", "degree", "Spline degree p >= 3.", degreeP);
    cmd.addInt("k", "regularity", "Interior regularity 1 <= k <= p-2 (-1 for k=p-2).", regularityK);
    cmd.addInt("r", "minRef", "Coarsest multigrid refinement level.", minRef);
    cmd.addInt("R", "maxRef", "Finest refinement level.", maxRef);
    cmd.addReal("", "young", "Young's modulus E.", E);
    cmd.addReal("", "nu", "Poisson ratio nu.", nu);
    cmd.addReal("", "thickness", "Shell thickness t.", thickness);
    cmd.addReal("w", "freq", "Manufactured solution wave frequency.", freq);
    cmd.addSwitch("plateZ", "Use purely transverse u_z manufactured wave (for flat plate checks).", plateOnlyZ);
    cmd.addSwitch("noQI", "Skip standalone Theorem 6 quasi-interpolant computation.", skipQI);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    const bool computeQI = !skipQI;

    gsKLShellMaterial<real_t> mat;
    mat.E = E;
    mat.nu = nu;
    mat.thickness = thickness;

    if (testMode == "verify_physics")
    {
        const std::string scordFile = (geomFile == "asg1_reparametrized_domains/farahat_kl_shells/xml/surfHalfSphere_asg1.xml")
                                    ? "filedata/surfaces/scordelis_lo_roof.xml"
                                    : geomFile;
        runPhysicsVerification(scordFile);
        return 0;
    }

    if (testMode == "scordelis_lo")
    {
        runScordelisLoBenchmark(geomFile, degreeP, regularityK, minRef, maxRef, solverType, smootherType);
        return 0;
    }

    gsMultiPatch<real_t> mp = loadOpenSurface3D(geomFile);

    if (testMode == "manufactured")
    {
        runManufacturedStudy(mp, mat, degreeP, regularityK, minRef, maxRef,
                             bcMode, solverType, smootherType, freq, plateOnlyZ,
                             computeQI, exportPrefix);
    }
    else if (testMode == "benchmark_solvers")
    {
        runSolverComparison(mp, mat, degreeP, regularityK, minRef, maxRef, bcMode, freq);
    }
    else if (testMode == "p_sweep")
    {
        for (short_t p = 3; p <= 6; ++p)
        {
            runManufacturedStudy(mp, mat, p, static_cast<short_t>(p - 2), minRef, maxRef,
                                 bcMode, solverType, smootherType, freq, plateOnlyZ,
                                 computeQI, "");
        }
    }
    else
    {
        GISMO_ERROR("Unknown --test mode: " << testMode);
    }

    return 0;
}
