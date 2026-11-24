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
    const std::string label = (bcMode == "simply_supported") ? "Values" : "ValuesAndDerivatives";
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

        gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat);
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

    gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat);
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
        {"schwarz_cg", "patch_skel"},
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
    gsKLShellMaterial<real_t> mat;
    mat.E = 4.32e8;
    mat.nu = 0.0;
    mat.thickness = 0.25;

    // Identify curved diaphragm boundaries (where X is constant at x=0 or x=50) vs free straight edges.
    // On the exact MacNeal-Harder / Kiendl Scordelis-Lo benchmark, u_y = u_z = 0 on both curved ends,
    // while axial warping u_x is FREE on both ends (only the 1-DOF rigid axial shift is removed by
    // pinning u_x = 0 at a single point via a tiny spring/corner elimination).
    gsBoundaryConditions<real_t> bc;
    gsConstantFunction<real_t> zero(0.0, 3);
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
        if (dZ > 1.0 && (dY < 1e-3 || dX < 1e-3))
        {
            if (dY < 1e-3)
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
    }

    gsInfo << "\n================================================================================\n"
           << "Scordelis-Lo Roof Benchmark (Non-Isoparametric Polynomial S_{p=" << degreeP
           << "} on Exact NURBS)\n"
           << "  Reference free-edge midpoint vertical displacement u_z ~ -0.3006 .. -0.3024\n"
           << "================================================================================\n";

    auto gravityLoad = [](const gsVector<real_t>&, gsVector<real_t, 3>& fBody, real_t& pNorm) {
        fBody << 0.0, 0.0, -90.0;
        pNorm = 0.0;
    };

    for (index_t r = minRef; r <= maxRef; ++r)
    {
        gsAsG1ShellHierarchy<real_t> hier(mp, bc, degreeP, regularityK, minRef, r, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat);
        gsSparseMatrix<real_t> K_discont, K_free;
        gsVector<real_t> F_discont, F_free, u_bnd, u_free;
        u_bnd.setZero(fin.mapper.boundarySize());
        asmbl.assembleStiffness(K_discont);
        asmbl.assemblePhysicalRHS(gravityLoad, F_discont);
        asmbl.projectSystem(K_discont, F_discont, fin.T_free, fin.T_bnd, u_bnd, K_free, F_free);

        // Pin the single axial rigid-body translation mode at corner 0 of patch 0
        const index_t corner0Off = sumUntil(fin.argBasis[0].sizes, 9);
        if (fin.mapper.is_free(corner0Off, 0, axialComp))
        {
            const index_t pinIdx = fin.mapper.index(corner0Off, 0, axialComp);
            K_free.coeffRef(pinIdx, pinIdx) += mat.E;
        }

        const auto st = hier.solve(K_free, F_free, mat, u_free, solverType, smootherType,
                                   1, 1, 1, 250, 1e-10);
        gsVector<real_t> c_sol = fin.T_free * u_free;

        // Evaluate vertical displacement u_z at midpoints of all 4 edges (free edge has max |u_z|)
        real_t minUz = 0.0;
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const real_t pts[4][2] = {{0.0, 0.5}, {1.0, 0.5}, {0.5, 0.0}, {0.5, 1.0}};
            for (int k = 0; k < 4; ++k)
            {
                const auto uVal = asmbl.evalDisplacement(c_sol, p, pts[k][0], pts[k][1]);
                if (uVal(2) < minUz) minUz = uVal(2);
            }
        }
        gsInfo << "  ref=" << r << "  nFree=" << std::setw(6) << st.nFree
               << "  u_z(free_mid)=" << std::fixed << std::setprecision(6) << minUz
               << "  iters=" << std::setw(4) << st.iterations
               << "  relRes=" << std::scientific << std::setprecision(1) << st.relResidual << "\n";
    }
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
    cmd.addString("t", "test", "Mode: manufactured | benchmark_solvers | p_sweep | scordelis_lo.", testMode);
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
