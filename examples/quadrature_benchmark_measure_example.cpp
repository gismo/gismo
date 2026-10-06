/** @file quadrature_benchmark_measure_example.cpp

    @brief Sweep all cut-cell quadrature methods across geometries and report
           area (2D) / volume (3D) accuracy and cost.

    One method (|sum w - exact|) x (geometry) table per refinement level, plus
    a final experimental-order-of-convergence table (one row per method).

    Methods:
      CutCell           gsQuadrature::CutCellRule (=11)
      Octree_L<k>       gsQuadrature::OctreeRule (=13), octLevels=k
      Algoim            gsAlgoimGenericRule
      Uniform_d<k>      gsAlgoimAdaptiveRule, indicator="uniform", maxDepth=k
      Divi_d<k>_t<tol>  gsAlgoimAdaptiveRule, indicator="integralChange"
      Moment_qA<k>      gsMomentRule over gsAlgoimAdaptiveRule, quA=k

    Usage:
      ./bin/quadrature_benchmark_measure_example --geo circle -r 4
      ./bin/quadrature_benchmark_measure_example --geo spain -f spain_mesh.msh -r 4
      ./bin/quadrature_benchmark_measure_example --geo cow -f obj/spot.obj -r 3

    This file is part of the G+Smo library.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  SpainLevelSet: signed-distance to a closed coastline polygon.
//  (Copied verbatim from immersed_fcm_octree_adaptive_surface_example2.cpp)
// =============================================================================

struct Seg { real_t x0, y0, x1, y1; };

static std::vector<Seg> polygonSegments(const std::vector<std::pair<real_t,real_t>>& pts)
{
    std::vector<Seg> segs;
    const size_t n = pts.size();
    segs.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % n];
        segs.push_back(Seg{a.first, a.second, b.first, b.second});
    }
    return segs;
}

static inline void closestOnSeg(real_t px, real_t py, const Seg& s, real_t& cx, real_t& cy)
{
    const real_t dx = s.x1 - s.x0, dy = s.y1 - s.y0;
    const real_t len2 = dx * dx + dy * dy;
    real_t t = (len2 > 0) ? ((px - s.x0) * dx + (py - s.y0) * dy) / len2 : (real_t)0;
    t = std::max<real_t>(0, std::min<real_t>(1, t));
    cx = s.x0 + t * dx;
    cy = s.y0 + t * dy;
}

static inline bool pointInside(real_t px, real_t py, const std::vector<Seg>& segs)
{
    int cn = 0;
    for (size_t i = 0; i < segs.size(); ++i)
    {
        const Seg& s = segs[i];
        const bool crosses = (s.y0 <= py) != (s.y1 <= py);
        if (crosses)
        {
            const real_t xint = s.x0 + (py - s.y0) / (s.y1 - s.y0) * (s.x1 - s.x0);
            if (px < xint) ++cn;
        }
    }
    return (cn & 1) != 0;
}

template<class T>
class SpainLevelSet : public gsFunction<T>
{
public:
    GISMO_CLONE_FUNCTION(SpainLevelSet)

    SpainLevelSet(std::vector<Seg> segs, const gsMatrix<T>& bbox)
        : m_segs(give(segs)), m_bbox(bbox) {}

    short_t     domainDim() const override { return 2; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            const T x = u(0, k), y = u(1, k);
            T best = std::numeric_limits<T>::max();
            for (size_t i = 0; i < m_segs.size(); ++i)
            {
                T cx, cy;
                closestOnSeg(x, y, m_segs[i], cx, cy);
                const T dx = x - cx, dy = y - cy;
                best = std::min(best, dx * dx + dy * dy);
            }
            const T dist = math::sqrt(best);
            result(0, k) = pointInside(x, y, m_segs) ? -dist : dist;
        }
    }

private:
    std::vector<Seg> m_segs;
    gsMatrix<T>      m_bbox;
};

static std::vector<std::pair<real_t,real_t>> loadPolygonTxt(const std::string& file)
{
    std::vector<std::pair<real_t,real_t>> pts;
    std::ifstream in(file.c_str());
    GISMO_ENSURE(in.good(), "Cannot open polygon text file: " << file);
    real_t x, y;
    while (in >> x >> y) pts.emplace_back(x, y);
    return pts;
}

static std::vector<std::pair<real_t,real_t>> loadPolygonGmsh(const std::string& file)
{
    std::ifstream in(file.c_str());
    GISMO_ENSURE(in.good(), "Cannot open Gmsh mesh file: " << file);

    std::string line;
    real_t version = 0;
    while (std::getline(in, line))
        if (line.rfind("$MeshFormat", 0) == 0)
        {
            in >> version;
            std::getline(in, line);
            break;
        }
    GISMO_ENSURE(version >= 4.0,
        "loadPolygonGmsh expects MSH 4.x; got version " << version);

    std::unordered_map<index_t, std::pair<real_t,real_t>> coord;
    while (std::getline(in, line))
        if (line.rfind("$Nodes", 0) == 0) break;

    index_t numBlocks = 0, numNodes = 0, minTag = 0, maxTag = 0;
    {
        std::getline(in, line);
        std::istringstream hs(line);
        hs >> numBlocks >> numNodes >> minTag >> maxTag;
    }
    coord.reserve(static_cast<size_t>(numNodes) * 2);
    for (index_t b = 0; b < numBlocks; ++b)
    {
        index_t edim, etag, parametric, nInBlock;
        std::getline(in, line);
        std::istringstream bh(line);
        bh >> edim >> etag >> parametric >> nInBlock;

        std::vector<index_t> tags(nInBlock);
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            std::istringstream ts(line);
            ts >> tags[i];
        }
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            std::istringstream cs(line);
            real_t x, y, z;
            cs >> x >> y >> z;
            coord[tags[i]] = std::make_pair(x, y);
        }
    }

    while (std::getline(in, line))
        if (line.rfind("$Elements", 0) == 0) break;

    index_t eNumBlocks = 0, eNumElems = 0, eMin = 0, eMax = 0;
    {
        std::getline(in, line);
        std::istringstream hs(line);
        hs >> eNumBlocks >> eNumElems >> eMin >> eMax;
    }

    std::vector<std::pair<index_t,index_t>> edges;
    for (index_t b = 0; b < eNumBlocks; ++b)
    {
        index_t edim, etag, etype, nInBlock;
        std::getline(in, line);
        std::istringstream bh(line);
        bh >> edim >> etag >> etype >> nInBlock;

        const bool boundaryLines = (edim == 1 && etype == 1);
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            if (!boundaryLines) continue;
            std::istringstream es(line);
            index_t etagId, n1, n2;
            es >> etagId >> n1 >> n2;
            edges.emplace_back(n1, n2);
        }
    }
    GISMO_ENSURE(!edges.empty(), "No boundary line elements (dim=1) found in " << file);

    std::unordered_map<index_t, std::vector<index_t>> adj;
    adj.reserve(edges.size() * 2);
    for (size_t i = 0; i < edges.size(); ++i)
    {
        adj[edges[i].first ].push_back(edges[i].second);
        adj[edges[i].second].push_back(edges[i].first );
    }

    const index_t start = edges.front().first;
    std::vector<std::pair<real_t,real_t>> pts;
    pts.reserve(edges.size() + 1);

    index_t prev = -1, cur = start;
    for (size_t guard = 0; guard <= edges.size() + 1; ++guard)
    {
        auto it = coord.find(cur);
        GISMO_ENSURE(it != coord.end(), "Boundary node " << cur << " missing from $Nodes.");
        pts.push_back(it->second);

        const auto& nb = adj[cur];
        index_t next = -1;
        for (size_t j = 0; j < nb.size(); ++j)
            if (nb[j] != prev) { next = nb[j]; break; }

        if (next == -1 || next == start) break;
        prev = cur;
        cur  = next;
    }
    return pts;
}

// =============================================================================
//  Method table
// =============================================================================

enum MethodType { MT_CUTCELL, MT_OCTREE, MT_ALGOIM, MT_UNIFORM, MT_DIVI, MT_MOMENT };

struct MethodConfig
{
    const char*  name;
    MethodType   type;
    int          octLevels     = 0;
    int          maxDepth      = 2;
    real_t       indicatorTol  = 1e-2;
    real_t       quA           = 1;
};

static const std::vector<MethodConfig> g_methods = {
    {"CutCell",         MT_CUTCELL},
    {"Octree_L1",        MT_OCTREE,  1},
    {"Octree_L2",        MT_OCTREE,  2},
    {"Octree_L3",        MT_OCTREE,  3},
    {"Algoim",          MT_ALGOIM},
    {"Uniform_d2",      MT_UNIFORM, 0, 2},
    {"Uniform_d3",      MT_UNIFORM, 0, 3},
    {"Divi_d2_t1e-2",   MT_DIVI,    0, 2, 1e-2},
    {"Divi_d2_t1e-3",   MT_DIVI,    0, 2, 1e-3},
    {"Divi_d3_t1e-2",   MT_DIVI,    0, 3, 1e-2},
    {"Divi_d3_t1e-3",   MT_DIVI,    0, 3, 1e-3},
    {"Moment_qA1",      MT_MOMENT,  0, 3, 1e-2, 1},
    {"Moment_qA2",      MT_MOMENT,  0, 3, 1e-2, 2},
};

// =============================================================================
//  Geometry descriptors
// =============================================================================

struct GeoInfo
{
    const char* key;
    int         dim;
    real_t      exactValue;
    bool        needsFile;
};

static const std::map<std::string, GeoInfo> g_geos = {
    {"circle", {"circle", 2, EIGEN_PI * 0.4 * 0.4, false}},
    {"sphere", {"sphere", 3, 4.0/3.0 * EIGEN_PI * 0.4 * 0.4 * 0.4, false}},
    {"spain",  {"spain",  2, -1, true}},
    {"cow",    {"cow",    3, 0.7182587880998567, true}},
};

// =============================================================================
//  Quadrature sweep helpers
// =============================================================================

// Integration closure: given a phi function and background cells, integrate 1.
template<class T, int D>
static T sweepMeasure(const gsFunction<T>& phi,
                      gsTensorBSplineBasis<D,T>& bkgBasis,
                      const MethodConfig& mc,
                      index_t deg,
                      index_t* nQuadPts = nullptr)
{
    gsImplicitTrimmedDomain<D,T> tr_domain(phi, bkgBasis);

    switch (mc.type)
    {
    case MT_CUTCELL:
    {
        gsOptionList opts;
        opts.addInt  ("quRule", "Quadrature rule id", gsQuadrature::CutCellRule);
        opts.addInt  ("quB",    "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal ("quA",    "quA: nodes = quA*deg + quB", 1.0);
        auto rule = gsQuadrature::getPtr<T>(tr_domain, opts);
        T area = 0;
        for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
        {
            gsMatrix<T> pts; gsVector<T> wts;
            rule->mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (nQuadPts) *nQuadPts += wts.size();
            for (index_t k = 0; k < wts.size(); ++k) area += wts[k];
        }
        return area;
    }
    case MT_OCTREE:
    {
        gsOptionList opts;
        opts.addInt  ("quRule",    "Quadrature rule id", gsQuadrature::OctreeRule);
        opts.addInt  ("octLevels", "Number of octree subdivision levels", mc.octLevels);
        opts.addInt  ("quB",       "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal ("quA",       "quA: nodes = quA*deg + quB", 1.0);
        auto rule = gsQuadrature::getPtr<T>(tr_domain, opts);
        T area = 0;
        for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
        {
            gsMatrix<T> pts; gsVector<T> wts;
            rule->mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (nQuadPts) *nQuadPts += wts.size();
            for (index_t k = 0; k < wts.size(); ++k) area += wts[k];
        }
        return area;
    }
    case MT_ALGOIM:
    {
        gsOptionList opts;
        opts.addInt ("quB", "quB: nodes = quA*deg + quB", deg + 1);
        opts.addReal("quA", "quA: nodes = quA*deg + quB", 1.0);
        gsAlgoimGenericRule<T> rule(phi, bkgBasis, opts);
        T area = 0;
        for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
        {
            gsMatrix<T> pts; gsVector<T> wts;
            rule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (nQuadPts) *nQuadPts += wts.size();
            for (index_t k = 0; k < wts.size(); ++k) area += wts[k];
        }
        return area;
    }
    case MT_UNIFORM:
    case MT_DIVI:
    {
        gsOptionList opts = gsAlgoimAdaptiveRule<T>::defaultOptions();
        opts.setInt  ("maxDepth",  mc.maxDepth);
        opts.setInt  ("nFallback", deg + 1);
        opts.setString("indicator", mc.type == MT_UNIFORM ? "uniform" : "integralChange");
        if (mc.type == MT_DIVI)
            opts.setReal("indicatorTol", mc.indicatorTol);
        gsAlgoimAdaptiveRule<T> rule(phi, bkgBasis, opts);
        T area = 0;
        for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
        {
            gsMatrix<T> pts; gsVector<T> wts;
            rule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (nQuadPts) *nQuadPts += wts.size();
            for (index_t k = 0; k < wts.size(); ++k) area += wts[k];
        }
        return area;
    }
    case MT_MOMENT:
    {
        gsOptionList opts = gsAlgoimAdaptiveRule<T>::defaultOptions();
        opts.setInt  ("maxDepth",  mc.maxDepth);
        opts.setInt  ("nFallback", deg + 1);
        opts.setString("indicator", "integralChange");
        opts.setReal ("indicatorTol", mc.indicatorTol);
        opts.setReal ("quA",     mc.quA);

        // Same output order as the deleted gsAlgoimMomentFittingRule::_numNodes1d():
        // n = round(quA*maxDegree) + quB, clamped to >= 1. NOT exactnessOrder(),
        // which returns 2*deg+1 and would redefine the Moment_qA1/qA2 columns.
        const long nRaw = std::lround(opts.askReal("quA", 1.0)
                                      * static_cast<double>(deg))
                        + static_cast<long>(opts.askInt("quB", 1));
        const index_t mfOrder1d = nRaw > 0 ? static_cast<index_t>(nRaw) : 1;

        // The compressor OWNS the adaptive rule: the gsQuadRule hierarchy has no
        // clone(), so ownership transfer is the only non-slicing option.
        typename gsMomentRule<T>::uPtr rule = gsMomentRule<T>::make(
            typename gsQuadRule<T>::uPtr(
                new gsAlgoimAdaptiveRule<T>(phi, bkgBasis, opts)),
            mfOrder1d);

        T area = 0;
        for (auto it = tr_domain.beginAll(); it != tr_domain.endAll(); ++it)
        {
            gsMatrix<T> pts; gsVector<T> wts;
            rule->mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            for (index_t k = 0; k < wts.size(); ++k) area += wts[k];
        }
        if (nQuadPts)
        {
            const auto& st = rule->stats();
            *nQuadPts += st.nOutputQPs + st.nPassThroughQPs;
        }
        return area;
    }
    }
    GISMO_ERROR("Unknown method type");
    return 0;
}

// Helper: filter methods by csv list (empty = all)
static std::vector<size_t> filterMethods(const std::string& csv)
{
    std::vector<size_t> idx;
    if (csv.empty())
    {
        for (size_t i = 0; i < g_methods.size(); ++i) idx.push_back(i);
        return idx;
    }

    std::istringstream iss(csv);
    std::string tok;
    while (std::getline(iss, tok, ','))
    {
        bool found = false;
        for (size_t i = 0; i < g_methods.size(); ++i)
        {
            if (g_methods[i].name == tok)
            {
                idx.push_back(i);
                found = true;
                break;
            }
        }
        if (!found) gsWarn << "Unknown method '" << tok << "' in --methods, ignored.\n";
    }
    return idx;
}

// =============================================================================
//  main
// =============================================================================

int main(int argc, char* argv[])
{
    std::string geoKey    = "circle";
    std::string filename  = "";
    std::string outDir    = "output_measure_benchmark";
    std::string methodsCsv= "";
    index_t     numRefine = 4;
    index_t     degree    = 2;
    real_t      fill      = 0.9;
    bool        plot      = false;

    gsCmdLine cmd("Quadrature benchmark: area/volume accuracy and cost sweep "
                  "across all cut-cell quadrature methods and geometries.");
    cmd.addString("",  "geo",     "Geometry: circle, sphere, spain, cow", geoKey);
    cmd.addString("f", "file",    "Input mesh file (required for spain, cow)", filename);
    cmd.addInt   ("r", "refine",  "Number of uniform refinement steps",    numRefine);
    cmd.addInt   ("e", "degree",  "B-spline degree of the background mesh",  degree);
    cmd.addReal  ("",  "fill",    "Fill fraction of [0,1]^d (0..1)",        fill);
    cmd.addString("o", "output",  "Output folder for results/plots",        outDir);
    cmd.addString("",  "methods", "CSV of method names to include (empty=all)", methodsCsv);
    cmd.addSwitch("plot", "Write ParaView output", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    auto gio = g_geos.find(geoKey);
    GISMO_ENSURE(gio != g_geos.end(), "Unknown geometry '" << geoKey
        << "'. Use circle, sphere, spain, or cow.");

    GeoInfo geo = gio->second;
    if (geo.needsFile)
        GISMO_ENSURE(!filename.empty(), "Geometry '" << geoKey << "' requires -f/--file.");

    const std::vector<size_t> methodIdx = filterMethods(methodsCsv);
    GISMO_ENSURE(!methodIdx.empty(), "No methods selected.");

    // -------------------------------------------------------------------------
    //  Setup geometry: phi function, exact value, invScale factors
    // -------------------------------------------------------------------------
    real_t exactVal = geo.exactValue;
    real_t scale = 1.0, invScaleVal = 1.0;  // invScale for value: v_phys = v_param * invScale

    std::unique_ptr<gsFunction<real_t>> phiPtr;
    std::function<void(index_t)> createBasis;

    if (geoKey == "circle")
    {
        real_t R = 0.4;
        exactVal = EIGEN_PI * R * R;
        scale = fill / (2 * R);
        invScaleVal = 1.0 / (scale * scale); // area
        std::string phiStr = "(x-0.5)^2 + (y-0.5)^2 - " + std::to_string(R*R);
        phiPtr = memory::make_unique(new gsFunctionExpr<real_t>(phiStr, 2));

        gsInfo << "Geometry: circle (R=" << R << ", fill=" << fill << ")\n";
        gsInfo << "Exact area (param) = " << std::setprecision(12) << exactVal
               << ", exact area (phys) = " << exactVal * invScaleVal << "\n\n";
    }
    else if (geoKey == "sphere")
    {
        real_t R = 0.4;
        exactVal = 4.0/3.0 * EIGEN_PI * R * R * R;
        scale = fill / (2 * R);
        invScaleVal = 1.0 / (scale * scale * scale); // volume
        std::string phiStr = "(x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2 - " + std::to_string(R*R);
        phiPtr = memory::make_unique(new gsFunctionExpr<real_t>(phiStr, 3));

        gsInfo << "Geometry: sphere (R=" << R << ", fill=" << fill << ")\n";
        gsInfo << "Exact volume (param) = " << std::setprecision(12) << exactVal
               << ", exact volume (phys) = " << exactVal * invScaleVal << "\n\n";
    }
    else if (geoKey == "spain")
    {
        const bool isMsh = filename.size() >= 4 &&
                           filename.compare(filename.size() - 4, 4, ".msh") == 0;
        auto pts = isMsh ? loadPolygonGmsh(filename) : loadPolygonTxt(filename);
        GISMO_ENSURE(pts.size() >= 3, "Coastline needs at least 3 vertices.");

        real_t xmin =  std::numeric_limits<real_t>::max(), ymin = xmin;
        real_t xmax = -std::numeric_limits<real_t>::max(), ymax = xmax;
        for (size_t i = 0; i < pts.size(); ++i)
        {
            xmin = std::min(xmin, pts[i].first);  xmax = std::max(xmax, pts[i].first);
            ymin = std::min(ymin, pts[i].second); ymax = std::max(ymax, pts[i].second);
        }
        const real_t cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax);
        const real_t extent = std::max(xmax - xmin, ymax - ymin);
        scale = fill / extent;

        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i].first  = (pts[i].first  - cx) * scale + 0.5;
            pts[i].second = (pts[i].second - cy) * scale + 0.5;
        }

        auto segs = polygonSegments(pts);
        exactVal = 0;
        for (size_t i = 0; i < segs.size(); ++i)
            exactVal += segs[i].x0 * segs[i].y1 - segs[i].x1 * segs[i].y0;
        exactVal = 0.5 * std::abs(exactVal);

        invScaleVal = 1.0 / (scale * scale);
        gsMatrix<real_t> bbox(2, 2);
        bbox.col(0) << 0, 0;
        bbox.col(1) << 1, 1;
        auto* spain = new SpainLevelSet<real_t>(segs, bbox);
        phiPtr = memory::make_unique(spain);

        gsInfo << "Geometry: spain (" << pts.size() << " vertices, fill=" << fill << ")\n";
        gsInfo << "Exact area (param) = " << std::setprecision(12) << exactVal
               << ", exact area (phys) = " << exactVal * invScaleVal << " (shoelace)\n\n";
    }
    else if (geoKey == "cow")
    {
        // Scoped: gsMeshSignedDist's BVH keeps its own copy of the triangles,
        // so the mesh need not outlive this branch.
        gsSurfMesh mesh;
        GISMO_ENSURE(gsReadSurfMesh(filename, mesh),
                     "Failed to read triangle mesh: " << filename);
        const std::size_t nVert = mesh.n_vertices();

        const real_t scale = gsNormalizeToUnitBox(mesh, fill);
        invScaleVal = 1.0 / (scale * scale * scale);

        exactVal = 0.7182587880998567;
        phiPtr = memory::make_unique(new gsMeshSignedDist<real_t>(mesh, gsUnitBox3()));

        gsInfo << "Geometry: cow (" << gsFileManager::getBasename(filename) << ", "
               << nVert << " verts, fill=" << fill << ")\n";
        gsInfo << "Reference volume (param) = " << std::setprecision(12) << exactVal
               << ", volume (phys) = " << exactVal * invScaleVal << "\n\n";
    }

    // -------------------------------------------------------------------------
    //  Header
    // -------------------------------------------------------------------------
    gsInfo << std::setw(14) << "method"
           << std::setw(16) << "value(param)"
           << std::setw(16) << "value(phys)"
           << std::setw(13) << "absErr"
           << std::setw(13) << "relErr"
           << std::setw(12) << "quadPts"
           << std::setw(12) << "time_ms\n";

    gsFileManager::mkdir(outDir);
    gsFileManager::mkdir(outDir + "/results");
    std::ofstream fout((outDir + "/results/measure.txt").c_str());
    fout << "# Quadrature benchmark: area/volume\n# geo=" << geoKey
         << " file=" << filename << " fill=" << fill
         << " degree=" << degree << " exact(param)=" << exactVal << "\n\n";

    // -------------------------------------------------------------------------
    //  Error storage for EoC
    // -------------------------------------------------------------------------
    gsMatrix<real_t> errMat(methodIdx.size(), numRefine + 1);

    // -------------------------------------------------------------------------
    //  Refinement loop
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        gsInfo << "\n--- Refinement " << r << " ---\n";

        if (geo.dim == 2)
        {
            const int nElem = 1 << r;
            gsKnotVector<> kvx(0, 1, nElem - 1, static_cast<short_t>(degree + 1), 1, static_cast<short_t>(degree));
            gsKnotVector<> kvy(0, 1, nElem - 1, static_cast<short_t>(degree + 1), 1, static_cast<short_t>(degree));
            gsTensorBSplineBasis<2,real_t> bkgBasis(kvx, kvy);
            const index_t deg = bkgBasis.maxDegree();

            for (size_t mi = 0; mi < methodIdx.size(); ++mi)
            {
                const auto& mc = g_methods[methodIdx[mi]];
                real_t val = 0;
                index_t nQP = 0;
                double tMs = 0;
                bool failed = false;

                try
                {
                    gsStopwatch clock;
                    val = sweepMeasure<real_t, 2>(*phiPtr, bkgBasis, mc, deg, &nQP);
                    tMs = clock.stop() * 1000.0;
                }
                catch (std::exception& e)
                {
                    gsInfo << "  " << std::setw(12) << mc.name << " FAILED: " << e.what() << "\n";
                    errMat(mi, r) = std::numeric_limits<real_t>::quiet_NaN();
                    continue;
                }

                const bool isNan = !std::isfinite(val);
                if (isNan)
                {
                    gsInfo << "  " << std::setw(12) << mc.name
                           << "  nan\n";
                    errMat(mi, r) = std::numeric_limits<real_t>::quiet_NaN();
                    continue;
                }

                const real_t absErr = std::abs(val - exactVal);
                const real_t relErr = (exactVal > 0) ? absErr / exactVal : 0;
                errMat(mi, r) = absErr;

                gsInfo << "  " << std::setw(12) << mc.name
                       << std::setw(16) << std::fixed << std::setprecision(8) << val
                       << std::setw(16) << std::fixed << std::setprecision(8) << val * invScaleVal
                       << std::setw(13) << std::scientific << std::setprecision(3) << absErr
                       << std::setw(13) << std::scientific << std::setprecision(3) << relErr
                       << std::setw(12) << nQP
                       << std::setw(12) << std::fixed << std::setprecision(2) << tMs << "\n";
            }
        }
        else // 3D
        {
            const int nElem = 1 << r;
            gsKnotVector<> kvx(0, 1, nElem - 1, static_cast<short_t>(degree + 1), 1, static_cast<short_t>(degree));
            gsKnotVector<> kvy(0, 1, nElem - 1, static_cast<short_t>(degree + 1), 1, static_cast<short_t>(degree));
            gsKnotVector<> kvz(0, 1, nElem - 1, static_cast<short_t>(degree + 1), 1, static_cast<short_t>(degree));
            gsTensorBSplineBasis<3,real_t> bkgBasis(kvx, kvy, kvz);
            const index_t deg = bkgBasis.maxDegree();

            for (size_t mi = 0; mi < methodIdx.size(); ++mi)
            {
                const auto& mc = g_methods[methodIdx[mi]];
                real_t val = 0;
                index_t nQP = 0;
                double tMs = 0;

                try
                {
                    gsStopwatch clock;
                    val = sweepMeasure<real_t, 3>(*phiPtr, bkgBasis, mc, deg, &nQP);
                    tMs = clock.stop() * 1000.0;
                }
                catch (std::exception& e)
                {
                    gsInfo << "  " << std::setw(12) << mc.name << " FAILED: " << e.what() << "\n";
                    errMat(mi, r) = std::numeric_limits<real_t>::quiet_NaN();
                    continue;
                }

                if (!std::isfinite(val))
                {
                    gsInfo << "  " << std::setw(12) << mc.name << "  nan\n";
                    errMat(mi, r) = std::numeric_limits<real_t>::quiet_NaN();
                    continue;
                }

                const real_t absErr = std::abs(val - exactVal);
                const real_t relErr = (exactVal > 0) ? absErr / exactVal : 0;
                errMat(mi, r) = absErr;

                gsInfo << "  " << std::setw(12) << mc.name
                       << std::setw(16) << std::fixed << std::setprecision(8) << val
                       << std::setw(16) << std::fixed << std::setprecision(8) << val * invScaleVal
                       << std::setw(13) << std::scientific << std::setprecision(3) << absErr
                       << std::setw(13) << std::scientific << std::setprecision(3) << relErr
                       << std::setw(12) << nQP
                       << std::setw(12) << std::fixed << std::setprecision(2) << tMs << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    //  EoC table
    // -------------------------------------------------------------------------
    if (numRefine > 0)
    {
        gsInfo << "\n\nExperimental order of convergence (|error| in area/volume):\n"
               << std::setw(14) << "method";
        for (int r = 1; r <= numRefine; ++r)
        {
            std::ostringstream oss;
            oss << "r" << (r-1) << "->" << r;
            gsInfo << std::setw(10) << oss.str();
        }
        gsInfo << "\n";

        for (size_t mi = 0; mi < methodIdx.size(); ++mi)
        {
            gsInfo << "  " << std::setw(12) << g_methods[methodIdx[mi]].name;
            for (int r = 1; r <= numRefine; ++r)
            {
                const real_t e0 = errMat(mi, r-1), e1 = errMat(mi, r);
                if (e0 > 0 && e1 > 0 && std::isfinite(e0) && std::isfinite(e1))
                {
                    const real_t rate = std::log(e0 / e1) / std::log(2.0);
                    gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << rate;
                }
                else
                    gsInfo << std::setw(10) << "  -";
            }
            gsInfo << "\n";
        }
    }

    fout.close();
    gsInfo << "\nResults written to: " << outDir << "/results/measure.txt\n";

    return EXIT_SUCCESS;
}
