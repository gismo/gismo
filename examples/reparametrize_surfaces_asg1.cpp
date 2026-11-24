/** @file reparametrize_surfaces_asg1.cpp

    @brief Batch utility to reparametrize 2D or 3D multi-patch geometries
           (B-spline or NURBS) into certified polynomial AS-G1 multi-patch
           surfaces using gsASG1LocalFitter, exporting both XML and ParaView
           (.pvd) previews along with AS-G1 residual and geometric error metrics.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsModeling/gsASG1LocalFitter.h>

#include <algorithm>
#include <cmath>
#include <dirent.h>
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

std::string stripPathAndExt(const std::string& path)
{
    std::string base = path;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos)
        base = base.substr(slash + 1);
    const size_t dot = base.find_last_of('.');
    if (dot != std::string::npos)
        base = base.substr(0, dot);
    return base;
}

bool endsWith(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> listXmlFiles(const std::string& dirPath)
{
    std::vector<std::string> files;
    DIR* dir = opendir(dirPath.c_str());
    if (!dir)
        return files;
    while (struct dirent* ent = readdir(dir))
    {
        std::string name(ent->d_name);
        if (name == "." || name == "..")
            continue;
        if (endsWith(name, ".xml"))
        {
            std::string full = dirPath;
            if (!full.empty() && full.back() != '/')
                full += '/';
            full += name;
            files.push_back(full);
        }
    }
    closedir(dir);
    std::sort(files.begin(), files.end());
    return files;
}

void normalizeSupportsToUnitSquare(gsMultiPatch<>& mp)
{
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        if (auto* bsp = dynamic_cast<gsTensorBSpline<2, real_t>*>(&mp.patch(p)))
        {
            bsp->basis().knots(0).affineTransformTo(real_t(0), real_t(1));
            bsp->basis().knots(1).affineTransformTo(real_t(0), real_t(1));
        }
        else if (auto* nurbs = dynamic_cast<gsTensorNurbs<2, real_t>*>(&mp.patch(p)))
        {
            nurbs->basis().source().knots(0).affineTransformTo(real_t(0), real_t(1));
            nurbs->basis().source().knots(1).affineTransformTo(real_t(0), real_t(1));
        }
    }
    mp.computeTopology();
}

bool readAsMultiPatch(const std::string& filename, gsMultiPatch<>& mp)
{
    mp.clear();
    gsFileData<> fd(filename);
    if (fd.has<gsMultiPatch<>>())
    {
        mp = *fd.getFirst<gsMultiPatch<>>();
        normalizeSupportsToUnitSquare(mp);
        return mp.nPatches() > 0;
    }
    if (fd.has<gsGeometry<>>())
    {
        std::vector<gsGeometry<>::uPtr> geoms = fd.getAll<gsGeometry<>>();
        for (size_t i = 0; i < geoms.size(); ++i)
            mp.addPatch(give(*geoms[i]));
        normalizeSupportsToUnitSquare(mp);
        return mp.nPatches() > 0;
    }
    return false;
}

struct AsG1CheckStats
{
    bool valid = true;
    real_t maxAlphaRes = 0;
    real_t maxBetaRes = 0;
    real_t maxNormalSin = 0;
};

AsG1CheckStats evaluateAsG1Stats(const gsMultiPatch<>& mp, real_t eps)
{
    AsG1CheckStats st;
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        auto S = sampleInterface(mp, *it, 5);
        st.maxNormalSin = std::max(st.maxNormalSin, S.maxNormalAngle);
        if (mp.targetDim() == 3 && S.maxNormalAngle * S.maxNormalAngle > eps)
            st.valid = false;
        if (S.D1.minCoeff() * S.D1.maxCoeff() <= 0 ||
            S.D2.minCoeff() * S.D2.maxCoeff() <= 0)
        {
            st.valid = false;
            continue;
        }
        const bool sameSign = (S.D1.minCoeff() * S.D2.minCoeff() > 0);
        if (sameSign) S.D1 = -S.D1;

        const index_t N = S.size();
        gsMatrix<real_t> Phi(N, 4);
        gsVector<real_t> p0 = 1.0 - S.nodes.row(0).array();
        gsVector<real_t> p1 = S.nodes.row(0).array();
        Phi.col(0) = p0.cwiseProduct(S.D2);
        Phi.col(1) = p1.cwiseProduct(S.D2);
        Phi.col(2) = p0.cwiseProduct(S.D1);
        Phi.col(3) = p1.cwiseProduct(S.D1);
        gsMatrix<real_t> G = Phi.transpose() * S.weights.asDiagonal() * Phi;
        const real_t tikh = 1e-10 * G.diagonal().cwiseAbs().maxCoeff();
        G(0,0) += tikh; G(1,1) += tikh; G(0,1) -= tikh; G(1,0) -= tikh;
        G(2,2) += tikh; G(3,3) += tikh; G(2,3) -= tikh; G(3,2) -= tikh;
        gsVector<real_t> c(4); c.setConstant(0.5);
        gsMatrix<real_t> K(5,5); K.setZero();
        K.block(0,0,4,4) = 2.0 * G;
        K.block(0,4,4,1) = c;
        K.block(4,0,1,4) = c.transpose();
        gsVector<real_t> r(5); r.setZero(); r(4) = 1.0;
        gsVector<real_t> alpha = K.fullPivLu().solve(r).head(4);
        const real_t resA = S.integrate((Phi * alpha).array().square());

        gsMatrix<real_t> Psi = Phi;
        Psi.rightCols(2) *= -1.0;
        gsMatrix<real_t> H = Psi.transpose() * S.weights.asDiagonal() * Psi;
        gsVector<real_t> d = Psi.transpose() * S.weights.asDiagonal() * S.D3;
        const real_t tikhB = 1e-10 * H.diagonal().cwiseAbs().maxCoeff();
        H(0,0) += tikhB; H(2,2) += tikhB; H(0,2) -= tikhB; H(2,0) -= tikhB;
        H(1,1) += tikhB; H(3,3) += tikhB; H(1,3) -= tikhB; H(3,1) -= tikhB;
        gsVector<real_t> e(4);
        gsVector<real_t> alpha1 = alpha[0] * p0 + alpha[1] * p1;
        gsVector<real_t> alpha2 = alpha[2] * p0 + alpha[3] * p1;
        e(0) =  S.integrate(alpha1.cwiseProduct(p0));
        e(1) =  S.integrate(alpha1.cwiseProduct(p1));
        e(2) = -S.integrate(alpha2.cwiseProduct(p0));
        e(3) = -S.integrate(alpha2.cwiseProduct(p1));
        gsMatrix<real_t> Kb(5,5); Kb.setZero();
        Kb.block(0,0,4,4) = 2.0 * H;
        Kb.block(0,4,4,1) = e;
        Kb.block(4,0,1,4) = e.transpose();
        gsVector<real_t> rb(5); rb.head(4) = 2.0 * d; rb(4) = 0.0;
        gsVector<real_t> beta = Kb.fullPivLu().solve(rb).head(4);
        const real_t resB = S.integrate((Psi * beta - S.D3).array().square());

        st.maxAlphaRes = std::max(st.maxAlphaRes, resA);
        st.maxBetaRes  = std::max(st.maxBetaRes,  resB);
        if (resA > eps || resB > eps)
            st.valid = false;
    }
    return st;
}

void computeGeometricError(const gsMultiPatch<>& target,
                           const gsMultiPatch<>& fitted,
                           real_t& maxErr,
                           real_t& relErr)
{
    const int nS = 11;
    gsMatrix<real_t> pts(2, nS * nS);
    for (int j = 0; j < nS; ++j)
        for (int i = 0; i < nS; ++i)
        {
            pts(0, i + j * nS) = real_t(i) / real_t(nS - 1);
            pts(1, i + j * nS) = real_t(j) / real_t(nS - 1);
        }

    const short_t gd = target.targetDim();
    gsVector<real_t> bbMin = gsVector<real_t>::Constant(gd,  1e30);
    gsVector<real_t> bbMax = gsVector<real_t>::Constant(gd, -1e30);
    maxErr = 0.0;

    for (size_t p = 0; p < target.nPatches(); ++p)
    {
        gsMatrix<real_t> vT = target.patch(p).eval(pts);
        gsMatrix<real_t> vF = fitted.patch(p).eval(pts);
        for (index_t k = 0; k < pts.cols(); ++k)
        {
            const real_t err = (vF.col(k) - vT.col(k)).norm();
            if (err > maxErr) maxErr = err;
            for (short_t d = 0; d < gd; ++d)
            {
                if (vT(d, k) < bbMin(d)) bbMin(d) = vT(d, k);
                if (vT(d, k) > bbMax(d)) bbMax(d) = vT(d, k);
            }
        }
    }
    const real_t diag = (bbMax - bbMin).norm();
    relErr = (diag > 1e-14) ? (maxErr / diag) : maxErr;
}

void exportPvd(gsMultiPatch<>& mp,
               const std::string& pvdDir,
               const std::string& baseName,
               index_t numPoints)
{
    std::string dir = pvdDir;
    if (!dir.empty() && dir.back() != '/')
        dir += '/';
    mkdirRecursive(dir + "_vtk");

    const std::string pvdBase = dir + baseName;
    try
    {
        gsMultiBasis<> mb(mp);
        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

        gsParaviewCollection col(pvdBase + ".pvd", &ev);
        col.options().setInt("numPoints", numPoints);
        col.options().setInt("precision", 6);
        col.options().setSwitch("plotElements", true);
        col.options().setInt("plotElements.resolution", -1);
        col.options().setSwitch("plotControlNet", false);
        col.options().setSwitch("makeSubfolder", true);
        col.options().setString("subfolder", "_vtk");

        col.newTimeStep(&mp);
        col.addField(meas(G), "Measure");
        if (mp.targetDim() == 3)
            col.addField(usn(G), "SurfaceNormal");
        col.saveTimeStep();
        col.save();
    }
    catch (...)
    {
        gsWriteParaview(mp, pvdBase, static_cast<unsigned>(numPoints), true, false);
    }
}

} // anonymous namespace

int main(int argc, char* argv[])
{
    std::string inputFile("");
    std::string inputDir("");
    std::string outDir("asg1_reparametrized_domains");
    int degree = 3;
    int regularity = 1;
    int innerKnots = 2;
    real_t eps = 1e-6;
    index_t numPoints = 1000;
    std::string modeStr = "global";

    gsCmdLine cmd("Reparametrize 2D/3D multi-patch geometries into certified AS-G1 domains.");
    cmd.addString("f", "file", "Single input geometry XML file.", inputFile);
    cmd.addString("d", "dir", "Directory of input geometry XML files.", inputDir);
    cmd.addString("o", "outdir", "Output directory (creates xml/ and pvd/ subfolders).", outDir);
    cmd.addInt("p", "degree", "Target AS-G1 spline degree p >= 3.", degree);
    cmd.addInt("r", "regularity", "Target interior regularity 1 <= r <= p-2.", regularity);
    cmd.addInt("k", "knots", "Number of inner knot intervals k.", innerKnots);
    cmd.addReal("e", "eps", "AS-G1 residual tolerance.", eps);
    cmd.addInt("n", "points", "ParaView sampling points per patch.", numPoints);
    cmd.addString("m", "mode", "Fitting mode: 'global' (default) or 'local'.", modeStr);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    std::vector<std::string> targets;
    if (!inputFile.empty())
        targets.push_back(inputFile);
    if (!inputDir.empty())
    {
        auto list = listXmlFiles(inputDir);
        targets.insert(targets.end(), list.begin(), list.end());
    }
    if (targets.empty())
    {
        gsInfo << "Usage: reparametrize_surfaces_asg1 -d <xml_dir> -o <outdir> [-p 3 -r 1 -k 2]\n";
        return EXIT_FAILURE;
    }

    const std::string xmlOutDir = outDir + "/xml";
    const std::string pvdOutDir = outDir + "/pvd";
    mkdirRecursive(xmlOutDir);
    mkdirRecursive(pvdOutDir);

    std::ofstream csv(outDir + "/reparametrization_report.csv");
    csv << "domain,patches,interfaces,geoDim,was_asg1,in_res_a,in_res_b,in_norm_sin,"
        << "fit_mode,is_asg1,out_res_a,out_res_b,out_norm_sin,max_geo_err,rel_geo_err\n";

    int nSuccess = 0;
    int nFailed = 0;

    gsInfo << "Reparametrizing " << targets.size() << " domain(s) with p=" << degree
           << ", r=" << regularity << ", k=" << innerKnots << " -> " << outDir << "/\n";

    for (size_t idx = 0; idx < targets.size(); ++idx)
    {
        const std::string base = stripPathAndExt(targets[idx]);
        gsMultiPatch<> target;
        try
        {
            if (!readAsMultiPatch(targets[idx], target))
            {
                gsWarn << "  [SKIP] " << base << " (no geometry)\n";
                continue;
            }
        }
        catch (const std::exception& e)
        {
            gsWarn << "  [SKIP] " << base << " (read failed: " << e.what() << ")\n";
            continue;
        }

        // Unify 3D surface normals / 2D Jacobian signs
        gsASG1LocalFitter<real_t>::unifyPatchOrientations(target);

        const AsG1CheckStats inStats = evaluateAsG1Stats(target, eps);

        // Run fitter
        gsMultiPatch<> fitted;
        bool okAsG1 = false;
        real_t outResA = 0, outResB = 0, outResN = 0;
        std::string usedMode = modeStr;

        try
        {
            gsASG1LocalFitter<real_t> fitter(target, degree, regularity, innerKnots);
            fitter.setMode(modeStr == "local"
                           ? gsASG1LocalFitter<real_t>::FittingMode::Local
                           : gsASG1LocalFitter<real_t>::FittingMode::Global);
            fitter.fit();
            okAsG1 = fitter.checkASG1(eps);
            outResA = fitter.maxAlphaResidual();
            outResB = fitter.maxBetaResidual();
            outResN = fitter.maxNormalResidual();
            fitted = fitter.result();

            // Fallback to Local if Global failed on a tricky geometry
            if (!okAsG1 && modeStr == "global")
            {
                gsASG1LocalFitter<real_t> fitterLoc(target, degree, regularity, innerKnots);
                fitterLoc.setMode(gsASG1LocalFitter<real_t>::FittingMode::Local);
                fitterLoc.fit();
                bool okLoc = fitterLoc.checkASG1(eps);
                if (okLoc || (fitterLoc.maxAlphaResidual() + fitterLoc.maxBetaResidual() < outResA + outResB))
                {
                    okAsG1 = okLoc;
                    outResA = fitterLoc.maxAlphaResidual();
                    outResB = fitterLoc.maxBetaResidual();
                    outResN = fitterLoc.maxNormalResidual();
                    fitted = fitterLoc.result();
                    usedMode = "local_fallback";
                }
            }
        }
        catch (const std::exception& e)
        {
            gsWarn << "  [FAIL] " << base << ": " << e.what() << "\n";
            ++nFailed;
            continue;
        }

        real_t maxErr = 0, relErr = 0;
        computeGeometricError(target, fitted, maxErr, relErr);

        // Save XML
        const std::string outXml = xmlOutDir + "/" + base + "_asg1.xml";
        gsFileData<real_t> outFd;
        outFd << fitted;
        outFd.save(outXml);

        // Save ParaView PVD
        exportPvd(fitted, pvdOutDir, base + "_asg1", numPoints);

        if (okAsG1) ++nSuccess; else ++nFailed;

        csv << base << "," << fitted.nPatches() << "," << fitted.nInterfaces() << ","
            << fitted.targetDim() << "," << (inStats.valid ? 1 : 0) << ","
            << std::scientific << std::setprecision(3)
            << inStats.maxAlphaRes << "," << inStats.maxBetaRes << "," << inStats.maxNormalSin << ","
            << usedMode << "," << (okAsG1 ? 1 : 0) << ","
            << outResA << "," << outResB << "," << outResN << ","
            << maxErr << "," << relErr << "\n";

        gsInfo << "  [" << (okAsG1 ? "OK  " : "WARN") << "] "
               << std::left << std::setw(40) << base
               << " P=" << std::setw(2) << fitted.nPatches()
               << " I=" << std::setw(2) << fitted.nInterfaces()
               << " inAS=" << (inStats.valid ? "Y" : "N")
               << " resA=" << std::scientific << std::setprecision(1) << outResA
               << " resB=" << outResB
               << " normSin=" << outResN
               << " relErr=" << relErr
               << " (" << usedMode << ")\n";
    }

    csv.close();
    gsInfo << "Summary: " << nSuccess << " certified AS-G1, " << nFailed << " warnings/failures.\n"
           << "Saved XMLs -> " << xmlOutDir << "/, PVDs -> " << pvdOutDir << "/\n";
    return 0;
}
