/** @file as_g1_geometry_generation_example.cpp
    @brief Example demonstrating the AS-G1 Multi-Patch Geometry Generator in G+Smo.

    Implements the "Bilinear Skeleton + Linearized Constrained Optimization" pipeline
    around an extraordinary vertex (e.g. valence V = 5) with curved outer boundaries.

    Author(s): F. Hasanova, Antigravity
*/

#include <gismo.h>
#include <gsModeling/gsAsG1GeometryGenerator.h>
#include <iostream>
#include <iomanip>

using namespace gismo;

int main(int argc, char* argv[])
{
    int valence = 5;
    int degree = 3;
    int numRefine = 0;
    real_t radius = 1.0;
    real_t midRadius = 1.2;
    std::string shape = "curved";
    std::string outPrefix = "as_g1_patch5";
    int numSamples = 100;
    bool writeFiles = true;
    bool randomize = false;
    int seed = 42;
    real_t anglePerturb = 0.15;
    real_t radiusPerturb = 0.2;
    std::string boundaryFile = "";

    gsCmdLine cmd("AS-G1 Multi-Patch Geometry Generator Example.");
    cmd.addInt("v", "valence", "Valence of central extraordinary vertex", valence);
    cmd.addInt("p", "degree", "B-spline degree (p >= 3)", degree);
    cmd.addInt("r", "refine", "Number of uniform knot refinements", numRefine);
    cmd.addReal("R", "radius", "Outer corner radius", radius);
    cmd.addReal("m", "midRadius", "Outer mid-point radius", midRadius);
    cmd.addString("s", "shape", "Boundary shape: curved, polygon, star", shape);
    cmd.addString("o", "output", "Output filename prefix", outPrefix);
    cmd.addString("b", "boundaryFile", "Path to boundary curves definition file (from Mathematica designer)", boundaryFile);
    cmd.addInt("q", "samples", "Number of test samples per interface", numSamples);
    cmd.addSwitch("w", "write", "Export Paraview and XML files", writeFiles);
    cmd.addSwitch("x", "random", "Generate randomized (asymmetric) geometry with perturbed control points", randomize);
    cmd.addInt("", "seed", "Random seed for perturbations (default: 42)", seed);
    cmd.addReal("", "angPerturb", "Magnitude of angle perturbation (fraction of sector angle, e.g. 0.15)", anglePerturb);
    cmd.addReal("", "radPerturb", "Magnitude of radius perturbation (fraction of radius, e.g. 0.2)", radiusPerturb);

    try {
        cmd.getValues(argc, argv);
    } catch (int) {
        return 0;
    }

    gsAsG1GeometryGenerator<real_t> generator(valence, degree, numRefine);
    if (!boundaryFile.empty()) {
        if (!generator.loadBoundaryDefinition(boundaryFile)) {
            gsWarn << "Failed to load boundary definition from: " << boundaryFile << "\n";
            return -1;
        }
        valence = generator.valence();
        degree = generator.degree();
    } else {
        generator.setRadius(radius);
        generator.setMidRadius(midRadius);
        generator.setRandomize(randomize, anglePerturb, radiusPerturb, seed);
    }

    gsInfo << "==================================================================\n";
    gsInfo << "  Analysis-Suitable G1 (AS-G1) Multi-Patch Geometry Generator    \n";
    gsInfo << "==================================================================\n";
    gsInfo << "Configuration:\n";
    gsInfo << "  Valence (patches): " << valence << "\n";
    gsInfo << "  B-spline degree:   " << degree << "\n";
    gsInfo << "  Knot refinements:  " << numRefine << " (" << (1 << numRefine) << " span(s) per direction)\n";
    if (!boundaryFile.empty()) {
        gsInfo << "  Boundary source:   External file (" << boundaryFile << ")\n";
    } else {
        gsInfo << "  Base radius R:     " << radius << "\n";
        gsInfo << "  Mid radius R_mid:  " << midRadius << "\n";
        gsInfo << "  Boundary shape:    " << shape << "\n";
        gsInfo << "  Symmetry mode:     " << (randomize ? "Randomized / Asymmetric" : "Symmetric") << "\n";
        if (randomize) {
            gsInfo << "    Seed:            " << seed << "\n";
            gsInfo << "    Angle perturb:   " << anglePerturb << "\n";
            gsInfo << "    Radius perturb:  " << radiusPerturb << "\n";
        }
    }
    gsInfo << "  Sample points:     " << numSamples << " per interface\n\n";

    // Setup custom boundary curves if requested
    std::vector<gsAsG1GeometryGenerator<real_t>::CurveFunction> eastCurves(valence, nullptr);
    std::vector<gsAsG1GeometryGenerator<real_t>::CurveFunction> northCurves(valence, nullptr);

    if (boundaryFile.empty()) {
        if (shape == "polygon") {
            // Straight chords
            for (int i = 0; i < valence; ++i) {
                eastCurves[i] = [&generator, i](real_t t) -> gsVector<real_t> {
                    // Linear segment
                    real_t th = 2.0 * M_PI * i / generator.valence();
                    real_t th_mid = 2.0 * M_PI * (i + 0.5) / generator.valence();
                    gsVector<real_t> p0(2), p1(2);
                    p0 << std::cos(th), std::sin(th);
                    p1 << 1.2 * std::cos(th_mid), 1.2 * std::sin(th_mid);
                    return (1.0 - t) * p0 + t * p1;
                };
            }
        } else if (shape == "star") {
            generator.setMidRadius(1.4 * radius);
        }
    }

    gsStopwatch timer;
    timer.restart();

    gsInfo << "[Step 1] Constructing bilinear skeleton...\n";
    generator.setupBilinearSkeleton();

    gsInfo << "[Step 2] Extracting exact linear gluing data (alpha, beta)...\n";
    generator.extractGluingData();
    const auto& gd = generator.gluingData();
    for (int i = 0; i < valence; ++i) {
        gsInfo << "  Interface " << i << ": D1=[ " << gd[i].D1_0 << " .. " << gd[i].D1_1
               << " ], D2=[ " << gd[i].D2_0 << " .. " << gd[i].D2_1
               << " ], b_glue=[ " << gd[i].b_glue(0.0) << " .. " << gd[i].b_glue(1.0) << " ]\n";
    }

    gsInfo << "[Step 3] Elevating degree to p=" << degree << " and refining basis...\n";
    generator.refineGeometry();

    gsInfo << "[Step 4] Assembling deformation stiffness matrix and linear constraints...\n";
    generator.assembleStiffnessMatrix();
    generator.assembleConstraints(eastCurves, northCurves);

    gsInfo << "[Step 5] Solving linear KKT system...\n";
    bool success = generator.solveKKT();
    double elapsed = timer.stop();

    if (!success) {
        gsInfo << "ERROR: KKT solve failed!\n";
        return 1;
    }
    gsInfo << "  Solve completed successfully in " << std::fixed << std::setprecision(4) << elapsed << " s.\n\n";

    gsInfo << "[Step 6] Verification Test Harness:\n";
    std::vector<real_t> res_list;
    real_t max_res = generator.verifyAsG1(numSamples, &res_list);

    gsInfo << "--- AS-G1 Interface Residuals across " << numSamples << " sample points ---\n";
    for (int i = 0; i < valence; ++i) {
        gsInfo << "  Interface " << i << " max residual: " << std::scientific << std::setprecision(4) << res_list[i] << "\n";
    }
    gsInfo << "  GLOBAL MAX RESIDUAL = " << std::scientific << std::setprecision(4) << max_res << "\n";

    const real_t tol = 1e-10;
    if (max_res < tol) {
        gsInfo << "  [SUCCESS] AS-G1 residual satisfies tolerance < " << std::scientific << tol
               << " (machine precision achieved!)\n\n";
    } else {
        gsInfo << "  [WARNING] AS-G1 residual exceeds tolerance " << std::scientific << tol << "\n\n";
    }

    gsInfo << "--- Patch Regularity Verification (det Jac > 0) ---\n";
    real_t minDet = 0, maxDet = 0;
    bool regular = generator.checkJacobianRegularity(25, &minDet, &maxDet);
    gsInfo << "  Grid sampled (25x25 per patch):\n";
    gsInfo << "  Minimum det(J) = " << std::fixed << std::setprecision(6) << minDet << "\n";
    gsInfo << "  Maximum det(J) = " << std::fixed << std::setprecision(6) << maxDet << "\n";
    if (regular) {
        gsInfo << "  [SUCCESS] All patches are strictly regular (det(J) > 0 everywhere)!\n\n";
    } else {
        gsInfo << "  [WARNING] Non-positive Jacobian determinant detected!\n\n";
    }

    if (writeFiles) {
        gsInfo << "[Step 7] Exporting geometry files...\n";
        std::string xmlName = outPrefix + ".xml";
        generator.exportXml(xmlName);
        gsInfo << "  XML file saved:      " << xmlName << "\n";

        generator.exportParaview(outPrefix, 1000);
        gsInfo << "  Paraview file saved: " << outPrefix << ".vts\n";
    }

    gsInfo << "==================================================================\n";
    gsInfo << "  AS-G1 Generation Complete!                                      \n";
    gsInfo << "==================================================================\n";

    return 0;
}
