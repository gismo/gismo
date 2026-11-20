/** @file planar_domain_generator_example.cpp
    @brief CLI for the 2D Planar Boundary Domain & Global-Distortion Generator in G+Smo.

    Usage examples:
      ./planar_domain_generator_example -p lshape_3holes -d weird -s 0.45 -o warped_lshape
      ./planar_domain_generator_example -p quad --expr-x "x+0.25*sin(2*y)" --expr-y "y+0.2*x^2" -o custom_warp
      ./planar_domain_generator_example --benchmark -d weird -s 0.4 -o filedata/domain2d/my_domains
      ./planar_domain_generator_example --multipatch -d swirl -s 0.4

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#include <gismo.h>
#include <gsModeling/gsPlanarDomainGenerator.h>

using namespace gismo;
using Gen = gsPlanarDomainGenerator<real_t>;

static bool exportDomain(const gsPlanarDomain<real_t>& domain, const std::string& prefix) {
    Gen::saveXml(domain, prefix + ".xml");
    auto mesh    = Gen::saveParaView(domain, prefix, 30);
    auto corners = Gen::extractCorners(domain);
    Gen::saveMathematica(domain, prefix + ".m");

    gsInfo << "    [OK] Loops: " << domain.numLoops()
           << " (1 outer CCW, " << domain.numHoles() << " holes CW)"
           << " | Sharp C^0 corners: " << corners.size()
           << " | Mesh: " << mesh.vertices.size() << " verts, "
           << mesh.elements.size() << " tris (Area = " << mesh.area() << ")\n"
           << "    [OK] Files: " << prefix << ".{xml,pvd,vtu,m}\n";
    return !mesh.elements.empty();
}

int main(int argc, char* argv[]) {
    std::string preset_name   = "curved_lshape";
    std::string output_prefix = "";
    std::string hole_shape    = "mixed";
    std::string distortion    = "none";
    std::string expr_x, expr_y;

    index_t num_corners = -1;
    index_t num_holes   = 0;
    index_t degree      = 3;
    index_t edge_waves  = 0;
    index_t warp_refine = 3;
    index_t seed        = 42;

    real_t curvature  = 0.3;
    real_t scale      = 1.0;
    real_t randomness = 0.0;
    real_t strength   = 0.75;
    real_t frequency  = 2.6;

    bool non_convex     = false;
    bool non_star       = false;
    bool curved         = false;
    bool run_benchmark  = false;
    bool run_multipatch = false;

    gsCmdLine cmd("G+Smo 2D Planar Boundary Domain & Global-Distortion Generator");
    cmd.addString("p", "preset",      "Preset name (circle, bean, amoeba, quad, lshape, plate_8holes, comb_4holes, ...)", preset_name);
    cmd.addInt   ("k", "corners",     "Custom corner count (>=0 overrides preset)", num_corners);
    cmd.addInt   ("H", "holes",       "Number of interior holes (custom mode)", num_holes);
    cmd.addString("",  "hole-shape",  "Hole shape: circle,ellipse,square,triangle,diamond,hexagon,star,cross,mixed", hole_shape);
    cmd.addSwitch("",  "non-convex",  "Generate non-convex custom domain", non_convex);
    cmd.addSwitch("",  "non-star",    "Generate non-star-shaped custom domain", non_star);
    cmd.addSwitch("",  "curved",      "Use curved edges in custom mode", curved);
    cmd.addReal  ("",  "curvature",   "Edge curvature factor (custom mode)", curvature);
    cmd.addInt   ("w", "waves",       "Edge wave / inflection count per chord (custom mode)", edge_waves);
    cmd.addReal  ("",  "randomness",  "Corner & hole placement randomness in [0, 1]", randomness);
    cmd.addReal  ("",  "scale",       "Domain scale factor", scale);
    cmd.addInt   ("",  "degree",      "B-spline degree", degree);
    cmd.addString("d", "distortion",  "Global warp: none,wave,swirl,bend,pinch,shear,weird,ffd", distortion);
    cmd.addReal  ("s", "strength",    "Global distortion strength", strength);
    cmd.addReal  ("",  "frequency",   "Spatial frequency for wave/swirl/weird distortion", frequency);
    cmd.addInt   ("r", "warp-refine", "Uniform curve refinements before Greville warp interpolation", warp_refine);
    cmd.addInt   ("",  "seed",        "Random seed for randomness & weird/ffd warps", seed);
    cmd.addString("",  "expr-x",      "Custom analytic warp x'(x,y) via gsFunctionExpr (overrides -d)", expr_x);
    cmd.addString("",  "expr-y",      "Custom analytic warp y'(x,y) via gsFunctionExpr (overrides -d)", expr_y);
    cmd.addString("o", "output",      "Output file prefix (single mode) or directory (benchmark/multipatch)", output_prefix);
    cmd.addSwitch("b", "benchmark",   "Run benchmark suite of distorted & multi-hole domains", run_benchmark);
    cmd.addSwitch("m", "multipatch",  "Generate non-bilinear multi-patch suite", run_multipatch);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    Gen::DistortionParams dp;
    dp.type      = distortion;
    dp.strength  = strength;
    dp.frequency = frequency;
    dp.refine    = static_cast<int>(warp_refine);
    dp.seed      = static_cast<unsigned>(seed);

    const bool hasExpr = (!expr_x.empty() && !expr_y.empty());
    if (hasExpr) dp.type = "none";

    if (run_benchmark) {
        const std::string out_dir = output_prefix.empty() ? "filedata/domain2d/generated_domains" : output_prefix;
        gsFileManager::mkdir(out_dir);

        struct BenchSpec { std::string name, preset, warp; real_t str; unsigned sd; };
        const std::vector<BenchSpec> cases = {
            {"uchannel_8corner_2holes",       "uchannel_2holes",      dp.active() ? dp.type : "weird", dp.active() ? dp.strength : 0.85, 101},
            {"quad_4corner_2holes_swirl",     "plate_with_2holes",    dp.active() ? dp.type : "swirl", dp.active() ? dp.strength : 0.95, 102},
            {"curved_lshape_6corner_1hole",   "curved_lshape_1hole",  dp.active() ? dp.type : "bend",  dp.active() ? dp.strength : 0.85, 103},
            {"plate_4corner_8holes_weird",    "plate_8holes",         dp.active() ? dp.type : "weird", dp.active() ? dp.strength : 0.80, 104},
            {"cross_12corner_5holes_ffd",     "cross_5holes",         dp.active() ? dp.type : "ffd",   dp.active() ? dp.strength : 0.90, 105},
            {"comb_16corner_4holes_wave",     "curved_comb_4holes",   dp.active() ? dp.type : "wave",  dp.active() ? dp.strength : 0.75, 106},
            {"bulged_quad_non_bilinear",      "bulged_quad",          dp.active() ? dp.type : "pinch", dp.active() ? dp.strength : 0.95, 107},
            {"warped_lshape_3holes",          "lshape_3holes",        dp.active() ? dp.type : "weird", dp.active() ? dp.strength : 0.90, 108},
            {"pinched_annulus_non_bilinear",  "pinched_annulus",      dp.active() ? dp.type : "pinch", dp.active() ? dp.strength : 0.90, 109},
            {"crescent_2corner_3holes_bend",  "crescent_3holes",      dp.active() ? dp.type : "bend",  dp.active() ? dp.strength : 0.85, 110},
            {"starburst_12corner_3holes_swirl","starburst_3holes",    dp.active() ? dp.type : "swirl", dp.active() ? dp.strength : 0.85, 111},
            {"amoeba_smooth_irregular",       "amoeba",               dp.active() ? dp.type : "weird", dp.active() ? dp.strength : 0.90, 112}
        };

        int passed = 0;
        for (const auto& bc : cases) {
            Gen::DistortionParams cdp = dp;
            cdp.type     = hasExpr ? "none" : bc.warp;
            cdp.strength = bc.str;
            cdp.seed     = bc.sd;

            gsInfo << "\n--> Generating: " << bc.name << " (preset=" << bc.preset
                   << ", warp=" << (hasExpr ? "expr" : cdp.type) << ", s=" << cdp.strength << ")\n";
            std::unique_ptr<gsPlanarDomain<real_t>> dom(
                Gen::preset(bc.preset, scale, static_cast<int>(degree), cdp));
            if (hasExpr) {
                gsFunctionExpr<real_t> phi(expr_x, expr_y, 2);
                dom.reset(Gen::distortDomain(*dom, phi, dp.refine, static_cast<int>(degree)));
            }
            if (exportDomain(*dom, out_dir + "/" + bc.name)) ++passed;
        }
        gsInfo << "\nBenchmark completed: " << passed << "/" << cases.size()
               << " domains saved in " << out_dir << "\n";
        return (passed == static_cast<int>(cases.size())) ? 0 : 1;
    }

    if (run_multipatch) {
        const std::string out_dir = output_prefix.empty() ? "filedata/domain2d/generated_domains" : output_prefix;
        gsFileManager::mkdir(out_dir);
        for (const std::string& name : {"curved_2patch", "curved_lshape_3patch", "curved_ring_4patch", "grid_4patch"}) {
            gsInfo << "\n--> MultiPatch: " << name << "...\n";
            std::unique_ptr<gsMultiPatch<real_t>> mp(
                Gen::generateNonBilinearMultiPatch(name, static_cast<int>(degree), dp));
            if (hasExpr) {
                gsFunctionExpr<real_t> phi(expr_x, expr_y, 2);
                mp.reset(Gen::distortMultiPatch(*mp, phi, 1, static_cast<int>(degree)));
            }
            const std::string prefix = out_dir + "/" + name;
            Gen::saveMultiPatchXml(*mp, prefix + ".xml");
            gsWriteParaview(*mp, prefix, 25);
            gsInfo << "    [OK] Patches: " << mp->nPatches() << " | Interfaces: " << mp->nInterfaces()
                   << " | Boundaries: " << mp->nBoundary() << " -> " << prefix << ".{xml,pvd}\n";
        }
        return 0;
    }

    // Single domain generation
    const std::string prefix = output_prefix.empty() ? "domain_output" : output_prefix;
    std::unique_ptr<gsPlanarDomain<real_t>> domain;
    if (num_corners >= 0) {
        Gen::Parameters p;
        p.num_corners = static_cast<int>(num_corners);
        p.convex      = !non_convex;
        p.star_shaped = !non_star;
        p.num_holes   = static_cast<int>(num_holes);
        p.hole_shape  = hole_shape;
        p.curved      = curved;
        p.curvature   = curvature;
        p.edge_waves  = static_cast<int>(edge_waves);
        p.randomness  = randomness;
        p.scale       = scale;
        p.degree      = static_cast<int>(degree);
        p.seed        = static_cast<unsigned>(seed);
        p.distortion  = dp;
        domain.reset(Gen::generate(p));
    } else {
        domain.reset(Gen::preset(preset_name, scale, static_cast<int>(degree), dp));
    }

    if (hasExpr) {
        gsFunctionExpr<real_t> phi(expr_x, expr_y, 2);
        domain.reset(Gen::distortDomain(*domain, phi, dp.refine, static_cast<int>(degree)));
    }

    gsInfo << "Generated domain (" << (num_corners >= 0 ? "custom" : preset_name)
           << ", warp=" << (hasExpr ? "expr" : dp.type) << "):\n";
    return exportDomain(*domain, prefix) ? 0 : 1;
}
