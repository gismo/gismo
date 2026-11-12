/** @file planar_domain_generator_example.cpp
    @brief Demonstration and CLI for the 2D Planar Boundary Domain Generator in G+Smo.

    Generates planar B-spline boundary curves, domains with holes, and non-bilinear
    multi-patch geometries for isogeometric analysis (IGA) and patch splitting.

    Exports:
    - G+Smo native <PlanarDomain> XML (.xml)
    - ParaView collection (.pvd) referencing:
      * _surface.vtp (triangulated 2D interior domain via robust Earcut)
      * _curves.vtp (B-spline boundary polylines)
      * _corners.vtp (C^0 corners)
      * _cnet.vtp (control polygons)
    - ParaView UnstructuredGrid (.vtu)
    - Mathematica visualization script (.m)

    Usage:
      ./planar_domain_generator_example -p <preset> -o <output_prefix>
      ./planar_domain_generator_example --benchmark
      ./planar_domain_generator_example --multipatch

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#include <gismo.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <iostream>
#include <string>
#include <vector>

using namespace gismo;

void printUsage(const char* prog) {
    std::cout << "========================================================================\n"
              << "  G+Smo 2D Planar Boundary Domain & Non-Bilinear Geometry Generator     \n"
              << "========================================================================\n"
              << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  -p, --preset <name>       Preset domain name (default: curved_lshape)\n"
              << "  -k, --corners <int>       Number of corners (C^0 reduction points)\n"
              << "  -H, --holes <int>         Number of interior holes (default: 0)\n"
              << "  --hole-shape <shape>      Hole shape: circle, square, triangle, diamond,\n"
              << "                            hexagon, star, cross, mixed (default: mixed)\n"
              << "  --convex <0|1>            Convexity flag (default: 1)\n"
              << "  --star <0|1>              Star-shaped flag (default: 1)\n"
              << "  --curved <0|1>            Curved non-bilinear boundaries (default: 0)\n"
              << "  --curvature <float>       Boundary curvature factor (default: 0.3)\n"
              << "  --scale <float>           Domain scale factor (default: 1.0)\n"
              << "  --degree <int>            B-spline degree (default: 3)\n"
              << "  -o, --output <prefix>     Output file prefix (default: domain_output)\n"
              << "  -b, --benchmark           Run full benchmark suite of non-bilinear domains\n"
              << "  -m, --multipatch          Generate non-bilinear multi-patch geometries\n"
              << "  -h, --help                Show this help message\n\n"
              << "Available Presets:\n"
              << "  Smooth (0 corners):   circle, ellipse, bean, horseshoe_smooth, wavy_ring\n"
              << "  1 Corner:             teardrop, heart\n"
              << "  2 Corners:            lens, crescent, crescent_3holes\n"
              << "  3 Corners:            triangle, curved_triangle\n"
              << "  4 Corners:            quad, curved_quad, bulged_quad, indented_quad,\n"
              << "                        plate_with_hole, plate_with_2holes, plate_5holes, plate_8holes\n"
              << "  5 Corners:            pentagon, curved_pentagon\n"
              << "  6 Corners (L-shape):  lshape, curved_lshape, lshape_with_hole, lshape_3holes\n"
              << "  8 Corners (U-channel):u_channel, curved_uchannel, uchannel_2holes, uchannel_4holes\n"
              << "  10 Corners:           star5\n"
              << "  12 Corners:           swiss_cross, bulged_cross, starburst, cross_5holes, starburst_3holes\n"
              << "  16 Corners:           comb, curved_comb, comb_4holes\n"
              << "========================================================================\n";
}

int runBenchmarkSuite(const std::string& out_dir) {
    std::cout << "\n========================================================================\n"
              << "  RUNNING FULL BENCHMARK SUITE (NON-BILINEAR & MULTI-HOLE GEOMETRIES)   \n"
              << "========================================================================\n";

    gsFileManager::mkdir(out_dir);

    struct BenchmarkCase {
        std::string name;
        std::string preset;
        std::string desc;
    };

    std::vector<BenchmarkCase> cases = {
        {"uchannel_8corner_2holes", "uchannel_2holes", "8-corner U-channel with 2 circular holes"},
        {"quad_4corner_2holes", "plate_with_2holes", "4-corner quad plate with 2 circular holes"},
        {"curved_lshape_6corner_1hole", "curved_lshape_1hole", "6-corner curved L-shape with square hole"},
        {"plate_4corner_8holes_weird", "plate_8holes", "4-corner plate with 8 diverse non-bilinear holes"},
        {"cross_12corner_5holes_star", "cross_5holes", "12-corner Swiss cross with 5 interior holes"},
        {"comb_16corner_4holes_curved", "comb_4holes", "16-corner 4-finger comb with 4 interior holes"},
        {"bulged_quad_non_bilinear", "bulged_quad", "4-corner bulged quad with curved cubic boundaries"},
        {"curved_lshape_non_bilinear", "curved_lshape", "6-corner curved non-bilinear L-shape"},
        {"wavy_ring_annulus", "wavy_ring", "0-corner wavy annular ring with circular hole"},
        {"crescent_2corner_3holes", "crescent_3holes", "2-corner crescent with 3 interior holes"},
        {"starburst_12corner_3holes", "starburst_3holes", "12-corner starburst with 3 interior holes"}
    };

    int success_count = 0;

    for (const auto& bc : cases) {
        std::cout << "\n--> Generating: " << bc.name << " (" << bc.desc << ")...\n";
        try {
            auto domain = gsPlanarDomainGenerator<real_t>::preset(bc.preset, 1.0, 3);

            std::string prefix = out_dir + "/" + bc.name;
            std::string xml_file = prefix + ".xml";
            std::string pvd_file = prefix + ".pvd";
            std::string m_file = prefix + ".m";

            // Save XML, ParaView suite, and Mathematica script
            gsPlanarDomainGenerator<real_t>::saveXml(*domain, xml_file);
            gsPlanarDomainGenerator<real_t>::saveParaView(*domain, prefix, 30);
            gsPlanarDomainGenerator<real_t>::saveMathematica(*domain, m_file);

            // Compute triangulation for mesh validation
            auto mesh = gsPlanarDomainGenerator<real_t>::triangulate(*domain, 30);
            auto corners = gsPlanarDomainGenerator<real_t>::extractCorners(*domain);

            std::cout << "    [OK] Loops: " << domain->numLoops()
                      << " (1 outer CCW, " << domain->numHoles() << " holes CW)"
                      << " | Corners: " << corners.size()
                      << " | Mesh: " << mesh.vertices.size() << " vertices, "
                      << mesh.elements.size() << " triangles (Area = " << mesh.area() << ")\n"
                      << "    [OK] Files: " << xml_file << ", " << pvd_file << ", " << m_file << "\n";

            delete domain;
            success_count++;
        } catch (const std::exception& e) {
            std::cerr << "    [ERROR] Failed to generate " << bc.name << ": " << e.what() << "\n";
        }
    }

    std::cout << "\n========================================================================\n"
              << "Benchmark completed: " << success_count << "/" << cases.size() << " passed successfully!\n"
              << "Output files located in: " << out_dir << "\n"
              << "========================================================================\n";
    return (success_count == static_cast<int>(cases.size())) ? 0 : 1;
}

int runMultiPatchSuite(const std::string& out_dir) {
    std::cout << "\n========================================================================\n"
              << "  GENERATING NON-BILINEAR MULTI-PATCH GEOMETRIES (gsMultiPatch)         \n"
              << "========================================================================\n";

    gsFileManager::mkdir(out_dir);

    std::vector<std::string> mp_presets = {
        "curved_2patch",
        "curved_lshape_3patch",
        "curved_ring_4patch"
    };

    for (const auto& name : mp_presets) {
        std::cout << "\n--> Generating MultiPatch: " << name << "...\n";
        try {
            auto mp = gsPlanarDomainGenerator<real_t>::generateNonBilinearMultiPatch(name, 3);

            std::string prefix = out_dir + "/" + name;
            std::string xml_file = prefix + ".xml";

            gsPlanarDomainGenerator<real_t>::saveMultiPatchXml(*mp, xml_file);
            gsWriteParaview(*mp, prefix, 25);

            std::cout << "    [OK] Patches: " << mp->nPatches()
                      << " | Interfaces: " << mp->nInterfaces()
                      << " | Boundaries: " << mp->nBoundary() << "\n"
                      << "    [OK] Saved XML: " << xml_file << "\n"
                      << "    [OK] Saved ParaView: " << prefix << ".pvd\n";

            delete mp;
        } catch (const std::exception& e) {
            std::cerr << "    [ERROR] Failed to generate " << name << ": " << e.what() << "\n";
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "-h" || arg1 == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }

    std::string preset_name = "curved_lshape";
    std::string output_prefix = "domain_output";
    int num_corners = -1;
    int num_holes = 0;
    std::string hole_shape = "mixed";
    bool convex = true;
    bool star_shaped = true;
    bool curved = false;
    real_t curvature = 0.3;
    real_t scale = 1.0;
    int degree = 3;
    bool run_benchmark = false;
    bool run_multipatch = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-p" || arg == "--preset") {
            if (i + 1 < argc) preset_name = argv[++i];
        } else if (arg == "-k" || arg == "--corners") {
            if (i + 1 < argc) num_corners = std::stoi(argv[++i]);
        } else if (arg == "-H" || arg == "--holes") {
            if (i + 1 < argc) num_holes = std::stoi(argv[++i]);
        } else if (arg == "--hole-shape") {
            if (i + 1 < argc) hole_shape = argv[++i];
        } else if (arg == "--convex") {
            if (i + 1 < argc) convex = (std::stoi(argv[++i]) != 0);
        } else if (arg == "--star") {
            if (i + 1 < argc) star_shaped = (std::stoi(argv[++i]) != 0);
        } else if (arg == "--curved") {
            if (i + 1 < argc) curved = (std::stoi(argv[++i]) != 0);
        } else if (arg == "--curvature") {
            if (i + 1 < argc) curvature = std::stod(argv[++i]);
        } else if (arg == "--scale") {
            if (i + 1 < argc) scale = std::stod(argv[++i]);
        } else if (arg == "--degree") {
            if (i + 1 < argc) degree = std::stoi(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) output_prefix = argv[++i];
        } else if (arg == "-b" || arg == "--benchmark") {
            run_benchmark = true;
        } else if (arg == "-m" || arg == "--multipatch") {
            run_multipatch = true;
        }
    }

    if (run_benchmark) {
        return runBenchmarkSuite("filedata/domain2d/generated_domains");
    }

    if (run_multipatch) {
        return runMultiPatchSuite("filedata/domain2d/generated_domains");
    }

    // Single domain generation
    std::cout << "========================================================================\n"
              << "  G+Smo Planar Boundary Domain Generator                                \n"
              << "========================================================================\n";

    gsPlanarDomain<real_t>* domain = nullptr;

    if (num_corners >= 0) {
        // Explicit parameter mode
        gsPlanarDomainGenerator<real_t>::Parameters p;
        p.num_corners = num_corners;
        p.convex = convex;
        p.star_shaped = star_shaped;
        p.num_holes = num_holes;
        p.hole_shape = hole_shape;
        p.curved = curved;
        p.curvature = curvature;
        p.scale = scale;
        p.degree = degree;
        domain = gsPlanarDomainGenerator<real_t>::generate(p);
        std::cout << "Generating custom domain: " << num_corners << " corners, "
                  << num_holes << " holes (" << hole_shape << "), curved=" << curved << "...\n";
    } else {
        // Preset mode
        std::cout << "Generating preset domain: " << preset_name << " (scale=" << scale << ", degree=" << degree << ")...\n";
        domain = gsPlanarDomainGenerator<real_t>::preset(preset_name, scale, degree);
    }

    std::string xml_file = output_prefix + ".xml";
    std::string pvd_file = output_prefix + ".pvd";
    std::string m_file = output_prefix + ".m";

    gsPlanarDomainGenerator<real_t>::saveXml(*domain, xml_file);
    gsPlanarDomainGenerator<real_t>::saveParaView(*domain, output_prefix, 30);
    gsPlanarDomainGenerator<real_t>::saveMathematica(*domain, m_file);

    auto mesh = gsPlanarDomainGenerator<real_t>::triangulate(*domain, 30);
    auto corners = gsPlanarDomainGenerator<real_t>::extractCorners(*domain);

    std::cout << "\nDomain Generated Successfully!\n"
              << "------------------------------------------------------------------------\n"
              << "Outer loop curves:    " << domain->outer().numCurves() << " (CCW)\n"
              << "Holes:                " << domain->numHoles() << " (CW)\n"
              << "C^0 Corners:          " << corners.size() << "\n"
              << "Triangulation Mesh:   " << mesh.vertices.size() << " vertices, "
              << mesh.elements.size() << " triangles (Area = " << mesh.area() << ")\n"
              << "G+Smo XML Output:     " << xml_file << "\n"
              << "ParaView PVD Output:  " << pvd_file << "\n"
              << "Mathematica Output:   " << m_file << "\n"
              << "========================================================================\n";

    delete domain;
    return 0;
}

