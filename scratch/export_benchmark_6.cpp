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
    const std::string initDir = "filedata/domain2d/benchmark_6_initial";
    const std::string splitDir = "filedata/domain2d/benchmark_6_split";
    gsFileManager::mkdir(initDir);
    gsFileManager::mkdir(splitDir);

    const std::vector<std::string> benchPresets = {
        "plate_with_hole",
        "lshape",
        "curved_lshape",
        "curved_quad",
        "circle",
        "quad"
    };

    const int degree = 3;
    int passed = 0;

    for (const auto& name : benchPresets)
    {
        std::unique_ptr<gsPlanarDomain<real_t>> dom(
            gsPlanarDomainGenerator<real_t>::preset(name, 1.0, degree));
        if (!dom)
        {
            std::cerr << "[ERROR] Failed to generate preset: " << name << "\n";
            continue;
        }

        // 1. Export initial non-split domain (ParaView + XML)
        const std::string initPrefix = initDir + "/" + name;
        gsPlanarDomainGenerator<real_t>::saveParaView(*dom, initPrefix, 40);
        {
            gsFileData<real_t> fd;
            fd << *dom;
            fd.save(initPrefix + ".xml");
        }

        // 2. Run splitter and export split domain (ParaView + XML)
        gsPatchSplitter<real_t> splitter(*dom);
        splitter.setStrategy(gsPatchSplitter<real_t>::Strategy::AutomaticBest);
        splitter.setDegree(degree);

        bool ok = splitter.computeSplitting();
        real_t minDet = 0;
        bool regular = splitter.checkJacobians(&minDet);
        bool conforming = splitter.isConforming();
        real_t bdrErr = 0;
        bool bdrOk = gsPatchSplitter<real_t>::verifyBoundaryMatch(*dom, splitter.getMultiPatch(), 1e-4, &bdrErr);

        const std::string splitPrefix = splitDir + "/" + name;
        splitter.exportParaview(splitPrefix);
        {
            gsFileData<real_t> fd;
            fd << splitter.getMultiPatch();
            fd.save(splitPrefix + ".xml");
        }

        if (ok && regular && conforming && bdrOk)
            ++passed;

        std::cout << name
                  << " | method=" << splitter.bestCandidate().methodUsed
                  << " | patches=" << splitter.getMultiPatch().nPatches()
                  << " | minDet=" << minDet
                  << " | bdrErr=" << bdrErr
                  << " | ok=" << (ok && regular && conforming && bdrOk ? "PASS" : "FAIL")
                  << "\n";
    }

    std::cout << "Exported " << passed << " / " << benchPresets.size() << " benchmark domains.\n";
    return (passed == static_cast<int>(benchPresets.size())) ? 0 : 1;
}
