/** @file as_g1_local_fitting_example.cpp
    @brief Example demonstrating the local AS-G1 multi-patch fitting algorithm.
    Reference: Farahat, Kapl, Kosmac, Vitrih (2024), Computers and Mathematics with Applications 168, 46-57.
*/

#include <gismo.h>
#include <gsModeling/gsASG1LocalFitter.h>
#include <iostream>
#include <iomanip>

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string inputFile = "filedata/domain2d/yeti_mp2.xml";
    std::string outputFile = "filedata/domain2d/as_g1_generator/yeti_asg1.xml";
    int degree = 3;
    int regularity = 1;
    int knots = 2;
    real_t eps = 1e-4;

    gsCmdLine cmd("Local Analysis-Suitable G1 (AS-G1) Multi-Patch Spline Surface Fitter.");
    cmd.addString("f", "file", "Input multi-patch geometry XML file.", inputFile);
    cmd.addString("o", "output", "Output fitted AS-G1 multi-patch XML file.", outputFile);
    cmd.addInt("p", "degree", "Target spline degree p >= 3.", degree);
    cmd.addInt("r", "regularity", "Target interface regularity 1 <= r <= p-2.", regularity);
    cmd.addInt("k", "knots", "Number of inner knot intervals k.", knots);
    cmd.addReal("e", "eps", "Tolerance for AS-G1 gluing residual check.", eps);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    gsInfo << "======================================================================\n";
    gsInfo << "Local Analysis-Suitable G1 (AS-G1) Multi-Patch Fitter\n";
    gsInfo << "Input file:  " << inputFile << "\n";
    gsInfo << "Output file: " << outputFile << "\n";
    gsInfo << "Parameters:  p = " << degree << ", r = " << regularity << ", k = " << knots << "\n";
    gsInfo << "======================================================================\n";

    gsMultiPatch<real_t> target;
    gsFileData<real_t> fd;
    fd.read(inputFile);
    if (!fd.has<gsMultiPatch<real_t>>()) {
        gsWarn << "File " << inputFile << " does not contain a multi-patch!\n";
        return 1;
    }
    fd.getFirst(target);
    target.computeTopology();

    gsInfo << "Target surface has " << target.nPatches() << " patches and "
           << target.nInterfaces() << " interfaces.\n";

    gsASG1LocalFitter<real_t> fitter(target, degree, regularity, knots);

    gsInfo << "Executing local 3-step Lagrange multiplier fitting...\n";
    fitter.fit();

    const gsMultiPatch<real_t>& fitted = fitter.result();
    gsInfo << "Fitting complete. Result has " << fitted.nPatches() << " patches and "
           << fitted.nInterfaces() << " interfaces.\n";

    bool is_asg1 = fitter.checkASG1(eps);
    gsInfo << "Max alpha residual: " << std::scientific << fitter.maxAlphaResidual() << "\n";
    gsInfo << "Max beta residual:  " << std::scientific << fitter.maxBetaResidual() << "\n";

    if (is_asg1) {
        gsInfo << ">>> SUCCESS: All interfaces satisfy the AS-G1 condition within tolerance "
               << eps << "! <<<\n";
    } else {
        gsInfo << ">>> WARNING: Some interfaces exceed AS-G1 tolerance " << eps << ". <<<\n";
    }

    gsFileData<real_t> outFd;
    outFd << fitted;
    outFd.save(outputFile);
    gsInfo << "Saved fitted AS-G1 geometry to: " << outputFile << "\n";

    return 0;
}

