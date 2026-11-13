#include <gismo.h>
#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace gismo;

int main()
{
    const std::string initDir = "filedata/domain2d/comb_no_holes_initial";
    const std::string splitDir = "filedata/domain2d/comb_no_holes_split";
    gsFileManager::mkdir(initDir);
    gsFileManager::mkdir(splitDir);

    const std::vector<std::string> presets = { "comb", "curved_comb" };
    const int degree = 3;
    int passed = 0;

    for (const auto& name : presets)
    {
        std::unique_ptr<gsPlanarDomain<real_t>> dom(
            gsPlanarDomainGenerator<real_t>::preset(name, 1.0, degree));
        if (!dom)
        {
            std::cerr << "[ERROR] Failed to generate preset: " << name << "\n";
            continue;
        }

        std::cout << "Preset " << name
                  << ": outer_curves=" << dom->outer().size()
                  << ", holes=" << dom->numHoles() << "\n";

        // Save initial non-split domain
        gsPlanarDomainGenerator<real_t>::saveParaView(*dom, initDir + "/" + name, 40);
        {
            gsFileData<real_t> fd;
            fd << *dom;
            fd.save(initDir + "/" + name + ".xml");
        }

        // Split domain
        gsPatchSplitter<real_t> splitter(*dom);
        splitter.setStrategy(gsPatchSplitter<real_t>::Strategy::AutomaticBest);
        splitter.setDegree(degree);

        bool splitOk = splitter.computeSplitting();
        const auto& mp = splitter.getMultiPatch();

        real_t bdrErr = 0;
        bool bdrOk = gsPatchSplitter<real_t>::verifyBoundaryMatch(*dom, mp, 1e-4, &bdrErr);
        bool confOk = splitter.isConforming();
        bool cvxOk = splitter.checkPatchCornersConvex();
        bool noIntOk = splitter.checkNoSelfIntersections();

        // Check per-patch nonzero single-sign detJ
        bool detSignOk = (mp.nPatches() > 0);
        real_t minAbsDet = 1e30;
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            real_t pMin = 1e30, pMax = -1e30, pMinAbs = 1e30;
            for (int j = 0; j <= 20; ++j)
            {
                for (int i = 0; i <= 20; ++i)
                {
                    gsMatrix<real_t> u(2, 1);
                    u(0, 0) = i / 20.0;
                    u(1, 0) = j / 20.0;
                    gsMatrix<real_t> J = mp.patch(p).jacobian(u);
                    real_t dJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
                    pMin = std::min(pMin, dJ);
                    pMax = std::max(pMax, dJ);
                    pMinAbs = std::min(pMinAbs, std::abs(dJ));
                }
            }
            minAbsDet = std::min(minAbsDet, pMinAbs);
            if (!((pMin > 1e-7) || (pMax < -1e-7)))
                detSignOk = false;
        }

        splitter.exportParaview(splitDir + "/" + name);
        {
            gsFileData<real_t> fd;
            fd << mp;
            fd.save(splitDir + "/" + name + ".xml");
        }

        bool allOk = splitOk && (mp.nPatches() > 0) && bdrOk && confOk && cvxOk && noIntOk && detSignOk;
        if (allOk) ++passed;

        std::cout << "  method=" << splitter.bestCandidate().methodUsed
                  << " | patches=" << mp.nPatches()
                  << " | bdrErr=" << bdrErr
                  << " | min|detJ|=" << minAbsDet
                  << " | conf=" << confOk
                  << " | detSign=" << detSignOk
                  << " | cvx=" << cvxOk
                  << " | noSelfInt=" << noIntOk
                  << " | status=" << (allOk ? "PASS" : "FAIL") << "\n";
    }

    std::cout << "Summary: " << passed << " / " << presets.size() << " comb domains passed.\n";
    return (passed == static_cast<int>(presets.size())) ? 0 : 1;
}
