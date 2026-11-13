#include <gismo.h>
#include <gsModeling/gsPatchSplitter.h>
#include <dirent.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace gismo;

struct Pt2 {
    real_t x, y;
};

static real_t cross2(const Pt2& a, const Pt2& b, const Pt2& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static bool isStrictlyConvexQuad(const Pt2 c[4], real_t& minTriArea2)
{
    // Each vertex must have the same strict orientation wrt the preceding edge;
    // Equivalently all 4 triangles (c[i], c[i+1], c[i+2]) have nonzero signed area of same sign,
    // which holds iff every corner is strictly outside the convex hull of the other 3.
    int sign = 0;
    minTriArea2 = 1e30;
    for (int i = 0; i < 4; ++i)
    {
        real_t cr = cross2(c[i], c[(i + 1) % 4], c[(i + 2) % 4]);
        real_t acr = std::abs(cr);
        if (acr < minTriArea2) minTriArea2 = acr;
        if (acr <= 1e-8) return false;
        int s = (cr > 0) ? 1 : -1;
        if (sign == 0) sign = s;
        else if (s != sign) return false;
    }
    return true;
}

static bool segProperIntersect(const Pt2& a, const Pt2& b, const Pt2& c, const Pt2& d)
{
    real_t d1 = cross2(a, b, c), d2 = cross2(a, b, d);
    real_t d3 = cross2(c, d, a), d4 = cross2(c, d, b);
    const real_t eps = 1e-9;
    return (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
            ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps)));
}

static bool pointInPolyStrict(const Pt2& p, const std::vector<Pt2>& poly)
{
    // Winding / ray casting; returns false if near boundary
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i)
    {
        const Pt2& a = poly[i];
        const Pt2& b = poly[(i + 1) % n];
        real_t dx = b.x - a.x, dy = b.y - a.y;
        real_t L2 = dx * dx + dy * dy;
        if (L2 > 1e-20)
        {
            real_t t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2;
            if (t < 0) t = 0; if (t > 1) t = 1;
            real_t qx = a.x + t * dx - p.x;
            real_t qy = a.y + t * dy - p.y;
            if (qx * qx + qy * qy < 1e-8) return false; // on boundary
        }
    }
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
            (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y + 1e-30) + poly[i].x))
        {
            inside = !inside;
        }
    }
    return inside;
}

static std::vector<Pt2> samplePatchBoundary(const gsGeometry<real_t>& patch, int nPerSide = 20)
{
    std::vector<Pt2> poly;
    gsMatrix<real_t> u(2, 1), val;
    auto pushPt = [&](real_t xi, real_t eta) {
        u(0, 0) = xi; u(1, 0) = eta;
        patch.eval_into(u, val);
        poly.push_back({val(0, 0), val(1, 0)});
    };
    for (int i = 0; i < nPerSide; ++i) pushPt(static_cast<real_t>(i) / nPerSide, 0.0);
    for (int i = 0; i < nPerSide; ++i) pushPt(1.0, static_cast<real_t>(i) / nPerSide);
    for (int i = 0; i < nPerSide; ++i) pushPt(1.0 - static_cast<real_t>(i) / nPerSide, 1.0);
    for (int i = 0; i < nPerSide; ++i) pushPt(0.0, 1.0 - static_cast<real_t>(i) / nPerSide);
    return poly;
}

int main(int argc, char** argv)
{
    std::string dirPath = (argc > 1) ? argv[1] : "filedata/domain2d/non_bilinear_no_holes_1to7corners";
    std::string outDir  = (argc > 2) ? argv[2] : "";

    std::vector<std::string> files;
    DIR* dir = opendir(dirPath.c_str());
    if (!dir)
    {
        std::cerr << "Cannot open directory: " << dirPath << "\n";
        return 1;
    }
    struct dirent* ent;
    while ((ent = readdir(dir)) != nullptr)
    {
        std::string name(ent->d_name);
        if (name.size() > 4 && name.substr(name.size() - 4) == ".xml")
            files.push_back(name);
    }
    closedir(dir);
    std::sort(files.begin(), files.end());

    int total = 0, passedAll = 0;
    std::cout << std::left
              << std::setw(36) << "Domain"
              << " | Split | P | BdrErr   | DetSignOK | QuadCvxOK | NoSelfInt | Status\n";
    std::cout << std::string(110, '-') << "\n";

    for (const auto& fname : files)
    {
        ++total;
        std::string fullPath = dirPath + "/" + fname;
        std::string base = fname.substr(0, fname.size() - 4);

        gsReadFile<real_t> rf(fullPath);
        memory::unique_ptr<gsPlanarDomain<real_t>> pDom = rf;
        if (!pDom)
        {
            std::cout << std::setw(36) << base << " | READ ERROR\n";
            continue;
        }
        const gsPlanarDomain<real_t>& dom = *pDom;

        gsPatchSplitter<real_t> splitter(dom);
        splitter.setStrategy(gsPatchSplitter<real_t>::Strategy::AutomaticBest);
        splitter.setDegree(3);
        bool splitOk = splitter.computeSplitting();
        const gsMultiPatch<real_t>& mp = splitter.getMultiPatch();

        int nP = static_cast<int>(mp.nPatches());
        real_t maxBdrDist = -1.0;
        bool bdrOk = (nP > 0) && gsPatchSplitter<real_t>::verifyBoundaryMatch(dom, mp, 1e-5, &maxBdrDist);

        // Per-patch det(J) non-zero & constant sign
        bool detSignOk = (nP > 0);
        real_t worstMinAbsDet = 1e30;
        bool anySignFlipWithinPatch = false;

        // Convex corners check
        bool quadCvxOk = (nP > 0);
        real_t worstTriArea = 1e30;

        // Sample boundary polygons for intersection/overlap check
        std::vector<std::vector<Pt2>> polys(nP);
        bool noSelfInt = (nP > 0);

        for (int p = 0; p < nP; ++p)
        {
            const auto& patch = mp.patch(p);
            real_t pMinDet = 1e30, pMaxDet = -1e30, pMinAbs = 1e30;
            const int N = 24;
            gsMatrix<real_t> u(2, 1), J;
            for (int j = 0; j <= N; ++j)
            {
                for (int i = 0; i <= N; ++i)
                {
                    u(0, 0) = static_cast<real_t>(i) / N;
                    u(1, 0) = static_cast<real_t>(j) / N;
                    J = patch.jacobian(u);
                    real_t dJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
                    pMinDet = std::min(pMinDet, dJ);
                    pMaxDet = std::max(pMaxDet, dJ);
                    pMinAbs = std::min(pMinAbs, std::abs(dJ));
                }
            }
            if (pMinAbs < worstMinAbsDet) worstMinAbsDet = pMinAbs;
            if (!((pMinDet > 1e-7) || (pMaxDet < -1e-7)))
            {
                detSignOk = false;
                if (pMinDet < -1e-9 && pMaxDet > 1e-9)
                    anySignFlipWithinPatch = true;
            }

            // 4 corners
            Pt2 c[4];
            gsMatrix<real_t> val;
            real_t uv[4][2] = {{0,0}, {1,0}, {1,1}, {0,1}};
            for (int k = 0; k < 4; ++k)
            {
                u(0, 0) = uv[k][0]; u(1, 0) = uv[k][1];
                patch.eval_into(u, val);
                c[k] = {val(0, 0), val(1, 0)};
            }
            real_t triA = 0;
            if (!isStrictlyConvexQuad(c, triA))
                quadCvxOk = false;
            if (triA < worstTriArea) worstTriArea = triA;

            // Self-intersection within patch boundary
            polys[p] = samplePatchBoundary(patch, 20);
            int m = static_cast<int>(polys[p].size());
            for (int i = 0; i < m && noSelfInt; ++i)
            {
                for (int j = i + 2; j < m; ++j)
                {
                    if (i == 0 && j == m - 1) continue;
                    if (segProperIntersect(polys[p][i], polys[p][(i + 1) % m],
                                           polys[p][j], polys[p][(j + 1) % m]))
                    {
                        noSelfInt = false;
                        break;
                    }
                }
            }
        }

        // Cross-patch boundary crossings & interior overlap check
        for (int p = 0; p < nP && noSelfInt; ++p)
        {
            for (int q = p + 1; q < nP && noSelfInt; ++q)
            {
                int mp1 = static_cast<int>(polys[p].size());
                int mq1 = static_cast<int>(polys[q].size());
                for (int i = 0; i < mp1 && noSelfInt; ++i)
                {
                    for (int j = 0; j < mq1; ++j)
                    {
                        if (segProperIntersect(polys[p][i], polys[p][(i + 1) % mp1],
                                               polys[q][j], polys[q][(j + 1) % mq1]))
                        {
                            noSelfInt = false;
                            break;
                        }
                    }
                }
                // Check interior sample points of patch p inside patch q
                gsMatrix<real_t> u(2, 1), val;
                real_t pts[3] = {0.25, 0.5, 0.75};
                for (int a = 0; a < 3 && noSelfInt; ++a)
                {
                    for (int b = 0; b < 3 && noSelfInt; ++b)
                    {
                        u(0, 0) = pts[a]; u(1, 0) = pts[b];
                        mp.patch(p).eval_into(u, val);
                        if (pointInPolyStrict({val(0, 0), val(1, 0)}, polys[q]))
                            noSelfInt = false;
                    }
                }
            }
        }

        // Check single-patch rule if domain is already a convex quad
        bool singleQuadRuleOk = true;
        if (dom.numHoles() == 0 && dom.outer().size() == 4)
        {
            Pt2 dc[4];
            for (int k = 0; k < 4; ++k)
            {
                gsMatrix<real_t> u0(1, 1);
                u0(0, 0) = dom.outer().curve(k).parameterRange()(0, 0);
                gsMatrix<real_t> v = dom.outer().curve(k).eval(u0);
                dc[k] = {v(0, 0), v(1, 0)};
            }
            real_t da = 0;
            if (isStrictlyConvexQuad(dc, da) && nP != 1)
                singleQuadRuleOk = false;
        }

        bool allOk = splitOk && (nP > 0) && bdrOk && detSignOk && quadCvxOk && noSelfInt && singleQuadRuleOk;
        if (allOk) ++passedAll;

        if (!outDir.empty() && nP > 0)
        {
            splitter.exportParaview(outDir + "/" + base);
            gsFileData<real_t> fd;
            fd << mp;
            fd.save(outDir + "/" + base + ".xml");
        }

        std::cout << std::setw(36) << base
                  << " | " << (splitOk ? " YES " : " NO  ")
                  << " | " << std::setw(2) << nP
                  << " | " << std::scientific << std::setprecision(1) << std::setw(8) << maxBdrDist
                  << " | " << (detSignOk ? "  PASS   " : (anySignFlipWithinPatch ? " FLIP!   " : " ZERO/NA "))
                  << " | " << (quadCvxOk ? "  PASS   " : "  FAIL   ")
                  << " | " << (noSelfInt ? "  PASS   " : " CROSS!  ")
                  << " | " << (allOk ? "PASS" : "FAIL")
                  << (!singleQuadRuleOk ? " (split convex quad!)" : "")
                  << "\n";
    }

    std::cout << std::string(110, '-') << "\n";
    std::cout << "Summary: " << passedAll << " / " << total << " domains passed all checks.\n";
    return 0;
}
