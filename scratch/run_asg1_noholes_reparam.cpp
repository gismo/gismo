#include <gismo.h>
#include <gsModeling/gsASG1LocalFitter.h>
#include <dirent.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

static real_t evalMinDetJ(const gsMultiPatch<real_t>& mp, int nSample = 21)
{
    real_t minDet = 1e30;
    gsMatrix<real_t> pts(2, nSample * nSample);
    int idx = 0;
    for (int j = 0; j < nSample; ++j) {
        real_t v = static_cast<real_t>(j) / (nSample - 1);
        for (int i = 0; i < nSample; ++i) {
            real_t u = static_cast<real_t>(i) / (nSample - 1);
            pts(0, idx) = u;
            pts(1, idx) = v;
            ++idx;
        }
    }
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        gsMatrix<real_t> jac;
        mp.patch(p).jacobian_into(pts, jac);
        for (index_t k = 0; k < pts.cols(); ++k) {
            real_t detJ = jac(0, 2 * k) * jac(1, 2 * k + 1) - jac(0, 2 * k + 1) * jac(1, 2 * k);
            if (detJ < minDet) minDet = detJ;
        }
    }
    return minDet;
}

static std::vector<std::pair<real_t, real_t>> sampleOuterBoundary(const gsMultiPatch<real_t>& mp, int nPerSide = 41)
{
    std::vector<std::pair<real_t, real_t>> res;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it) {
        const gsGeometry<real_t>& patch = mp.patch(it->patch);
        boxSide s = it->side();
        gsMatrix<real_t> uv(2, nPerSide);
        for (int k = 0; k < nPerSide; ++k) {
            real_t t = static_cast<real_t>(k) / (nPerSide - 1);
            if (s.direction() == 0) {
                uv(0, k) = s.parameter() ? 1.0 : 0.0;
                uv(1, k) = t;
            } else {
                uv(0, k) = t;
                uv(1, k) = s.parameter() ? 1.0 : 0.0;
            }
        }
        gsMatrix<real_t> pts;
        patch.eval_into(uv, pts);
        for (int k = 0; k < nPerSide; ++k)
            res.emplace_back(pts(0, k), pts(1, k));
    }
    return res;
}

static real_t computeBoundaryDiff(const gsMultiPatch<real_t>& a, const gsMultiPatch<real_t>& b)
{
    auto pa = sampleOuterBoundary(a);
    auto pb = sampleOuterBoundary(b);
    auto maxOneWay = [](const std::vector<std::pair<real_t, real_t>>& src,
                        const std::vector<std::pair<real_t, real_t>>& dst) -> real_t {
        real_t worst = 0;
        for (const auto& p : src) {
            real_t best2 = 1e30;
            for (const auto& q : dst) {
                real_t dx = p.first - q.first;
                real_t dy = p.second - q.second;
                real_t d2 = dx * dx + dy * dy;
                if (d2 < best2) best2 = d2;
            }
            worst = std::max(worst, std::sqrt(best2));
        }
        return worst;
    };
    return std::max(maxOneWay(pa, pb), maxOneWay(pb, pa));
}

static real_t computeInteriorRmsDiff(const gsMultiPatch<real_t>& a, const gsMultiPatch<real_t>& b, int nSample = 13)
{
    gsMatrix<real_t> pts(2, nSample * nSample);
    int idx = 0;
    for (int j = 0; j < nSample; ++j) {
        real_t v = static_cast<real_t>(j) / (nSample - 1);
        for (int i = 0; i < nSample; ++i) {
            real_t u = static_cast<real_t>(i) / (nSample - 1);
            pts(0, idx) = u;
            pts(1, idx) = v;
            ++idx;
        }
    }
    real_t sum2 = 0;
    index_t count = 0;
    for (size_t p = 0; p < a.nPatches(); ++p) {
        gsMatrix<real_t> va, vb;
        a.patch(p).eval_into(pts, va);
        b.patch(p).eval_into(pts, vb);
        for (index_t k = 0; k < pts.cols(); ++k) {
            real_t dx = va(0, k) - vb(0, k);
            real_t dy = va(1, k) - vb(1, k);
            sum2 += dx * dx + dy * dy;
            ++count;
        }
    }
    return std::sqrt(sum2 / std::max<index_t>(1, count));
}

static std::vector<std::string> listXmls(const std::string& dir)
{
    std::vector<std::string> files;
    DIR* dp = opendir(dir.c_str());
    if (!dp) return files;
    while (struct dirent* ent = readdir(dp)) {
        std::string name(ent->d_name);
        if (name.size() > 4 && name.substr(name.size() - 4) == ".xml")
            files.push_back(name);
    }
    closedir(dp);
    std::sort(files.begin(), files.end());
    return files;
}

int main(int argc, char* argv[])
{
    int p = 3;
    int r = 1;
    int k = 2;
    real_t eps = 1e-4;
    std::string outRoot = "filedata/domain2d/as_g1_no_holes_reparametrized";
    if (argc > 1) k = std::stoi(argv[1]);
    if (argc > 2) outRoot = argv[2];

    gsFileManager::mkdir(outRoot);

    struct Suite {
        std::string label;
        std::string dir;
        std::vector<std::string> exclude;
    };
    std::vector<Suite> suites = {
        {"benchmark_6", "filedata/domain2d/benchmark_6_split", {"plate_with_hole.xml"}},
        {"comb_no_holes", "filedata/domain2d/comb_no_holes_split", {}},
        {"1to7corners", "filedata/domain2d/split_no_holes_1to7corners", {}},
        {"noHolesCounterExamples", "filedata/domain2d/split_noHolesCounterExamples", {}}
    };

    std::ofstream csv(outRoot + "/asg1_reparametrization_report.csv");
    csv << "suite,name,patches,interfaces,reoriented,asg1_ok,max_alpha_res,max_beta_res,min_detJ_in,min_detJ_asg1,bdr_diff,interior_rms,xml_path\n";

    std::cout << std::left
              << std::setw(24) << "Suite"
              << std::setw(36) << "Geometry"
              << std::setw(4)  << "P"
              << std::setw(4)  << "IF"
              << std::setw(8)  << "ASG1"
              << std::setw(11) << "MaxRes"
              << std::setw(12) << "MinDetJ"
              << std::setw(11) << "BdrDiff"
              << std::setw(11) << "FitRMS"
              << "\n";
    std::cout << std::string(121, '-') << "\n";

    int total = 0, asg1Pass = 0, detPass = 0;

    for (const auto& s : suites) {
        std::string subOut = outRoot + "/" + s.label;
        gsFileManager::mkdir(subOut);
        for (const auto& fn : listXmls(s.dir)) {
            if (std::find(s.exclude.begin(), s.exclude.end(), fn) != s.exclude.end())
                continue;
            ++total;
            std::string base = fn.substr(0, fn.size() - 4);
            std::string inPath = s.dir + "/" + fn;

            gsMultiPatch<real_t> target;
            gsFileData<real_t> fd;
            fd.read(inPath);
            fd.getFirst(target);
            target.computeTopology();

            bool reor = gsASG1LocalFitter<real_t>::unifyPatchOrientations(target);
            real_t minDetIn = evalMinDetJ(target);

            gsASG1LocalFitter<real_t> fitter(target, p, r, k);
            fitter.setMode(gsASG1LocalFitter<real_t>::FittingMode::Global);
            fitter.fit();

            bool okASG1 = fitter.checkASG1(eps);
            real_t resA = fitter.maxAlphaResidual();
            real_t resB = fitter.maxBetaResidual();
            real_t maxRes = std::max(resA, resB);

            const gsMultiPatch<real_t>& fitted = fitter.result();
            real_t minDetOut = evalMinDetJ(fitted);
            real_t bdrDiff = computeBoundaryDiff(target, fitted);
            real_t rmsDiff = computeInteriorRmsDiff(target, fitted);

            if (okASG1) ++asg1Pass;
            if (minDetOut > 0) ++detPass;

            std::string outXml = subOut + "/" + base + "_asg1.xml";
            gsFileData<real_t> outFd;
            outFd << fitted;
            outFd.save(outXml);
            gsWriteParaview(fitted, subOut + "/" + base + "_asg1", 900, true);

            std::cout << std::left
                      << std::setw(24) << s.label
                      << std::setw(36) << base
                      << std::setw(4)  << fitted.nPatches()
                      << std::setw(4)  << fitted.nInterfaces()
                      << std::setw(8)  << (okASG1 ? "PASS" : "FAIL")
                      << std::scientific << std::setprecision(2)
                      << std::setw(11) << maxRes
                      << std::setw(12) << minDetOut
                      << std::setw(11) << bdrDiff
                      << std::setw(11) << rmsDiff
                      << "\n";

            csv << s.label << "," << base << ","
                << fitted.nPatches() << "," << fitted.nInterfaces() << ","
                << (reor ? 1 : 0) << "," << (okASG1 ? 1 : 0) << ","
                << std::scientific << std::setprecision(6)
                << resA << "," << resB << ","
                << minDetIn << "," << minDetOut << ","
                << bdrDiff << "," << rmsDiff << ","
                << outXml << "\n";
        }
    }

    std::cout << std::string(121, '-') << "\n";
    std::cout << "Total: " << total
              << " | AS-G1 pass: " << asg1Pass << "/" << total
              << " | minDetJ > 0: " << detPass << "/" << total << "\n";
    return 0;
}
