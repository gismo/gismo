/** @file fischerkolmogrov_prostate_example_bodyfitted.cpp

    @brief Fisher-Kolmogorov tumour growth on a body-fitted, PATIENT-SPECIFIC
           prostate volume (single TENSOR B-spline patch), obtained by mapping
           a solid torus onto a triangulated prostate-with-urethra surface.

    Geometry (Lorenzo et al., PNAS 2016). A solid torus with a POLAR
    cross-section is the parametric domain; degree p (default 2),
    nCirc x nCirc x nRad elements (default 32 x 32 x 8):
        u : toroidal angle, around the urethra       (seam u = 0 = 1)
        v : poloidal angle, around the cross-section (seam v = 0 = 1)
        w : radial; w = 0 is the core line (collapsed polar axis), w = 1 the
            surface (prostate capsule AND urethra wall: the hole is the urethra)

    Construction (see fischerkolmogrov_prostate_common.h):
      1. Urethra centreline, straightened coordinates (s, r, phi).
      2. Half-plane cross-sections {phi = 2 pi u}; polar rays theta = 2 pi v
         from an interior centre c(phi) give the boundary b(phi, theta); the
         interior is x = Centreline(c + w (b - c), phi). With --fit ls, v is
         the ARC LENGTH of the section contour (--curvWeight 0, default): the
         polar angle (--curvWeight -1) puts very few points at the far ends
         of the elongated sections, i.e. at the apex/base rims where the
         urethra leaves the prostate.
      3. --fit interp : interpolation at the Greville points (legacy, 34 x 34
                        boundary samples, raw centres, C0 seams).
         --fit ls     : (default) least-squares fit of densely sampled
                        sections (--fitSamples per element), centre c(phi)
                        smoothed by a Fourier series, C1 across the u- and
                        v-seams, core line exactly collapsed.
      4. --fitIts k  : k closest-point iterations that pull the w = 1 face onto
                        the triangulation (input vertices at their foot points
                        + spline samples at their closest points), regularized
                        by the parametric Dirichlet energy, which also extends
                        the correction harmonically into the volume. Updates
                        that would fold the volume are halved.
      5. --refine r  : r uniform knot insertions (32x32x8 -> 256x256x64 for r = 3).
    The scaled Jacobian is checked at the quadrature points of the coarse and
    refined mesh; the surface deviation (Hausdorff, mean, RMS) to the
    triangulation is printed and written to ParaView.

    Solver: generalized-alpha + Newton (fischerkolmogrov_torus_example_bodyfitted.cpp),
    seams glued C0, the DoFs of the collapsed face w = 0 merged per u.
    UNITS: mm, day (the OBJ is in mm).

    Geometry check, legacy vs improved fit:
      ./build/bin/fischerkolmogrov_prostate_example_bodyfitted -r 1 --geometryOnly --fit interp --fitIts 0 -o 00_prostate_interp
      ./build/bin/fischerkolmogrov_prostate_example_bodyfitted -r 1 --geometryOnly -o 00_prostate_ls

    Best fit: ./build/bin/fischerkolmogrov_prostate_example_bodyfitted -r 1 --geometryOnly -o 00_prostate_ls

    Full run (256 x 256 x 64, ~4.4M DoFs: run on the server):
      ./build/bin/fischerkolmogrov_prostate_example_bodyfitted \
          -D 1.26e-3 -p 1.86e-3 -t 1 -N 730 -R 2.5 -s 1 -a 0.5 -r 3 \
          -o 00_output_FK_prostate_bf --plot

    See also fischerkolmogrov_prostate_example_thb.cpp (adaptive THB fit).

    Author(s): L. Venta-Viñuela
*/
// =============================================================================
// Self-contained: the geometry fit, quality checks and the FK solver below were
// fischerkolmogrov_prostate_common.h (namespace gismo::fkp).
// =============================================================================
#include <gismo.h>
#include <array>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace gismo {
namespace fkp {

typedef gsVector3d<real_t> Vec3;
typedef gsVector<real_t, 2> Vec2;

inline real_t wrapPi(real_t a) { return math::atan2(math::sin(a), math::cos(a)); }
inline real_t wrap01(real_t a) { return a - math::floor(a); }

inline gsMatrix<> toMatrix(const std::vector<Vec3> &pts)
{
    gsMatrix<> m(3, pts.size());
    for (size_t i = 0; i < pts.size(); ++i)
        m.col(i) = pts[i];
    return m;
}

// =============================================================================
// TRIANGULATED SURFACE
// =============================================================================
/// Closest point to p on the triangle (a,b,c) (Ericson, Real-Time Collision Detection)
inline Vec3 closestOnTriangle(const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c)
{
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const real_t d1 = ab.dot(ap), d2 = ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) return a;
    const Vec3 bp = p - b;
    const real_t d3 = ab.dot(bp), d4 = ac.dot(bp);
    if (d3 >= 0 && d4 <= d3) return b;
    const real_t vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + (d1 / (d1 - d3)) * ab;
    const Vec3 cp = p - c;
    const real_t d5 = ab.dot(cp), d6 = ac.dot(cp);
    if (d6 >= 0 && d5 <= d6) return c;
    const real_t vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + (d2 / (d2 - d6)) * ac;
    const real_t va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);
    const real_t den = 1 / (va + vb + vc);
    return a + (vb * den) * ab + (vc * den) * ac;
}

struct TriMesh
{
    std::vector<Vec3> verts;                    // unique vertices
    std::vector<std::array<index_t, 3> > faces; // outward oriented
    std::vector<Vec3> tris;                     // triangle soup (3 per face)
    std::vector<Vec3> nrm, ctr;                 // unit normal, centroid per face
    std::vector<real_t> rad;                    // bounding radius around ctr
    real_t vol = 0, diag = 0;
    Vec3 bbMin, bbMax;

    size_t nTri() const { return faces.size(); }

    void read(const std::string &file)
    {
        std::vector<Vec3> soup;
        const std::string ext = gsFileManager::getExtension(file);
        const bool ok = (ext == "stl" || ext == "STL") ? readStl(file, soup) : readObj(file, soup);
        GISMO_ENSURE(ok && !soup.empty() && soup.size() % 3 == 0, "Could not read triangles from " << file);

        // merge identical vertices
        std::map<std::array<real_t, 3>, index_t> id;
        faces.resize(soup.size() / 3);
        for (size_t i = 0; i < soup.size(); ++i)
        {
            const std::array<real_t, 3> key = {{soup[i][0], soup[i][1], soup[i][2]}};
            auto it = id.find(key);
            index_t k;
            if (it == id.end()) { k = verts.size(); id[key] = k; verts.push_back(soup[i]); }
            else k = it->second;
            faces[i / 3][i % 3] = k;
        }
        vol = 0;
        for (size_t f = 0; f < faces.size(); ++f)
            vol += verts[faces[f][0]].dot(verts[faces[f][1]].cross(verts[faces[f][2]])) / 6;
        if (vol < 0) // make the normals point outwards
        {
            for (size_t f = 0; f < faces.size(); ++f) std::swap(faces[f][1], faces[f][2]);
            vol = -vol;
        }
        bbMin = bbMax = verts[0];
        for (size_t i = 0; i < verts.size(); ++i) { bbMin = bbMin.cwiseMin(verts[i]); bbMax = bbMax.cwiseMax(verts[i]); }
        diag = (bbMax - bbMin).norm();
        tris.clear(); nrm.clear(); ctr.clear(); rad.clear();
        for (size_t f = 0; f < faces.size(); ++f)
        {
            const Vec3 &a = verts[faces[f][0]], &b = verts[faces[f][1]], &c = verts[faces[f][2]];
            tris.push_back(a); tris.push_back(b); tris.push_back(c);
            Vec3 n = (b - a).cross(c - a);
            const real_t nn = n.norm();
            nrm.push_back(nn > 0 ? Vec3(n / nn) : Vec3(Vec3::Zero()));
            const Vec3 g = (a + b + c) / 3.0;
            ctr.push_back(g);
            rad.push_back(math::max((a - g).norm(), math::max((b - g).norm(), (c - g).norm())));
        }
    }

    /// Closest point q on the surface; returns the SIGNED distance (> 0 outside)
    real_t closest(const Vec3 &p, Vec3 &q) const
    {
        real_t best = std::numeric_limits<real_t>::max();
        index_t fb = 0;
        for (size_t f = 0; f < faces.size(); ++f)
        {
            if ((p - ctr[f]).norm() - rad[f] > best) continue; // bounding-sphere pruning
            const Vec3 c = closestOnTriangle(p, tris[3 * f], tris[3 * f + 1], tris[3 * f + 2]);
            const real_t d = (p - c).norm();
            if (d < best) { best = d; q = c; fb = f; }
        }
        return (p - q).dot(nrm[fb]) >= 0 ? best : -best;
    }

    /// Writes the triangulation with a per-vertex field
    void writeVtp(const std::string &fn, const std::vector<real_t> &field = std::vector<real_t>(),
                  const std::string &name = "field") const
    {
        std::ofstream f((fn + ".vtp").c_str());
        f << std::setprecision(10)
          << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\">\n<PolyData>\n"
          << "<Piece NumberOfPoints=\"" << verts.size() << "\" NumberOfPolys=\"" << faces.size() << "\">\n";
        if (!field.empty())
        {
            f << "<PointData Scalars=\"" << name << "\">\n<DataArray type=\"Float64\" Name=\"" << name
              << "\" format=\"ascii\">\n";
            for (size_t i = 0; i < field.size(); ++i) f << field[i] << " ";
            f << "\n</DataArray>\n</PointData>\n";
        }
        f << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < verts.size(); ++i) f << verts[i][0] << " " << verts[i][1] << " " << verts[i][2] << " ";
        f << "\n</DataArray>\n</Points>\n<Polys>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < faces.size(); ++i) f << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << " ";
        f << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t i = 0; i < faces.size(); ++i) f << 3 * (i + 1) << " ";
        f << "\n</DataArray>\n</Polys>\n</Piece>\n</PolyData>\n</VTKFile>\n";
    }

private:
    static bool readObj(const std::string &filename, std::vector<Vec3> &tris)
    {
        std::ifstream in(filename.c_str());
        if (!in) return false;
        std::vector<Vec3> verts;
        std::vector<std::vector<long> > faces;
        std::string line, tok;
        while (std::getline(in, line))
        {
            if (line.size() < 2 || (line[1] != ' ' && line[1] != '\t')) continue;
            std::istringstream iss(line.substr(2));
            if (line[0] == 'v')
            {
                real_t x, y, z;
                if (iss >> x >> y >> z) verts.push_back(Vec3(x, y, z));
            }
            else if (line[0] == 'f')
            {
                std::vector<long> ids;
                while (iss >> tok) ids.push_back(std::atol(tok.c_str())); // stops at the first '/'
                faces.push_back(ids);
            }
        }
        tris.clear();
        const long nv = static_cast<long>(verts.size());
        for (size_t f = 0; f < faces.size(); ++f)
        {
            std::vector<long> ids = faces[f];
            for (size_t k = 0; k < ids.size(); ++k)
            {
                if (ids[k] < 0) ids[k] += nv + 1;
                if (ids[k] < 1 || ids[k] > nv) return false;
            }
            for (size_t k = 1; k + 1 < ids.size(); ++k) // fan triangulation
            {
                tris.push_back(verts[ids[0] - 1]);
                tris.push_back(verts[ids[k] - 1]);
                tris.push_back(verts[ids[k + 1] - 1]);
            }
        }
        return !tris.empty();
    }

    static bool readStl(const std::string &filename, std::vector<Vec3> &tris)
    {
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (!in) return false;
        in.seekg(0, std::ios::end);
        const std::streamoff size = in.tellg();
        in.seekg(0, std::ios::beg);
        if (size < 84) return false;
        char header[80];
        in.read(header, 80);
        uint32_t triCount = 0;
        in.read(reinterpret_cast<char *>(&triCount), 4);
        tris.clear();
        if (84 + static_cast<std::streamoff>(triCount) * 50 == size) // binary
        {
            for (uint32_t t = 0; t < triCount; ++t)
            {
                float buf[12];
                uint16_t attr;
                in.read(reinterpret_cast<char *>(buf), 12 * sizeof(float));
                in.read(reinterpret_cast<char *>(&attr), sizeof(uint16_t));
                if (!in) break;
                for (int k = 1; k < 4; ++k) tris.push_back(Vec3(buf[3 * k], buf[3 * k + 1], buf[3 * k + 2]));
            }
            return !tris.empty();
        }
        in.close();
        std::ifstream txt(filename.c_str());
        std::string line;
        const std::string key = "vertex ";
        while (std::getline(txt, line))
        {
            const std::size_t pos = line.find(key);
            if (pos == std::string::npos) continue;
            std::istringstream iss(line.substr(pos + key.size()));
            real_t x, y, z;
            if (iss >> x >> y >> z) tris.push_back(Vec3(x, y, z));
        }
        return !tris.empty() && tris.size() % 3 == 0;
    }
};

// Nearest hit t > 0 of the ray o + t*d with the triangle soup (Moeller-Trumbore)
inline real_t rayFirstHit(const std::vector<Vec3> &tris, const Vec3 &o, const Vec3 &d, index_t &hitTri)
{
    real_t tBest = -1;
    hitTri = -1;
    for (size_t i = 0; i + 2 < tris.size(); i += 3)
    {
        const Vec3 e1 = tris[i + 1] - tris[i], e2 = tris[i + 2] - tris[i];
        const Vec3 pv = d.cross(e2);
        const real_t det = e1.dot(pv);
        if (math::abs(det) < 1e-14) continue;
        const Vec3 tv = o - tris[i];
        const real_t u = tv.dot(pv) / det;
        if (u < -1e-9 || u > 1 + 1e-9) continue;
        const Vec3 qv = tv.cross(e1);
        const real_t v = d.dot(qv) / det;
        if (v < -1e-9 || u + v > 1 + 1e-9) continue;
        const real_t t = e2.dot(qv) / det;
        if (t > 0 && (tBest < 0 || t < tBest)) { tBest = t; hitTri = static_cast<index_t>(i / 3); }
    }
    return tBest;
}

// =============================================================================
// CENTRELINE AND CROSS-SECTIONS
// =============================================================================
// Sampled centreline with arc length and a rotation-minimizing frame.
struct Centerline
{
    std::vector<Vec3> P, T, N, B;
    std::vector<real_t> s;

    void frame(real_t sq, Vec3 &p, Vec3 &n, Vec3 &b) const
    {
        const size_t K = s.size();
        size_t k = std::upper_bound(s.begin(), s.end(), sq) - s.begin();
        k = std::min(std::max<size_t>(k, 1), K - 1);
        const real_t a = (sq - s[k - 1]) / (s[k] - s[k - 1]); // may extrapolate
        p = (1 - a) * P[k - 1] + a * P[k];
        n = ((1 - a) * N[k - 1] + a * N[k]).normalized();
        b = ((1 - a) * B[k - 1] + a * B[k]).normalized();
    }

    Vec3 map(real_t sq, real_t r, real_t phi) const
    {
        Vec3 p, n, b;
        frame(sq, p, n, b);
        return p + r * (math::cos(phi) * n + math::sin(phi) * b);
    }

    void project(const Vec3 &x, real_t &sq, real_t &r, real_t &phi) const
    {
        size_t kb = 0;
        real_t db = std::numeric_limits<real_t>::max();
        for (size_t k = 0; k < P.size(); ++k)
        {
            const real_t d = (x - P[k]).squaredNorm();
            if (d < db) { db = d; kb = k; }
        }
        sq = s[kb];
        for (size_t k = (kb > 0 ? kb - 1 : 0); k < std::min(kb + 1, P.size() - 1); ++k)
        {
            const Vec3 e = P[k + 1] - P[k];
            const real_t a = math::min((real_t)1, math::max((real_t)0, (x - P[k]).dot(e) / e.squaredNorm()));
            const real_t d = (x - (P[k] + a * e)).squaredNorm();
            if (d <= db) { db = d; sq = s[k] + a * (s[k + 1] - s[k]); }
        }
        Vec3 p, n, b;
        frame(sq, p, n, b);
        const Vec3 rv = x - p;
        r = math::sqrt(math::pow(rv.dot(n), 2) + math::pow(rv.dot(b), 2));
        phi = math::atan2(rv.dot(b), rv.dot(n));
    }
};

typedef std::vector<gsVector<real_t, 4> > Segments; // (s0 r0 s1 r1) in a half-plane

struct ProstateSurface
{
    TriMesh m;
    Centerline cl;
    std::vector<Vec3> lumenPts;
    std::vector<real_t> vs, vr, vphi; // straightened coordinates of the soup vertices
    real_t kappaMax = 0, rMax = 0;

    void build(const std::string &file, index_t clDegree = 3, real_t lumenMax = 0)
    {
        m.read(file);
        const std::vector<Vec3> &tris = m.tris;
        const size_t nTri = m.nTri();
        const real_t diag = m.diag;
        if (lumenMax <= 0) lumenMax = 0.1 * diag;
        gsInfo << "Surface: " << file << ", " << nTri << " triangles, " << m.verts.size()
               << " vertices, enclosed volume " << m.vol << " mm^3, bbox diagonal " << diag << " mm\n";

        // Lumen triangles: the outward normal hits the opposite wall (facing back)
        // within one lumen diameter; their midpoints sample the urethra.
        lumenPts.clear();
        const size_t triStride = std::max<size_t>(1, nTri / 3000);
        for (size_t t = 0; t < nTri; t += triStride)
        {
            const Vec3 n = m.nrm[t];
            if (n.norm() < 0.5) continue;
            const Vec3 g = m.ctr[t];
            index_t hit;
            const real_t d = rayFirstHit(tris, g + 1e-7 * diag * n, n, hit);
            if (d <= 0 || d > lumenMax) continue;
            if (m.nrm[hit].dot(n) < 0) lumenPts.push_back(g + 0.5 * d * n);
        }
        GISMO_ENSURE(lumenPts.size() >= 10, "No urethra found (is the surface genus 1? try --lumenMax).");

        Vec3 m0 = Vec3::Zero();
        for (size_t i = 0; i < lumenPts.size(); ++i) m0 += lumenPts[i] / lumenPts.size();
        gsMatrix<real_t, 3, 3> cov = gsMatrix<real_t, 3, 3>::Zero();
        for (size_t i = 0; i < lumenPts.size(); ++i) cov += (lumenPts[i] - m0) * (lumenPts[i] - m0).transpose();
        Vec3 axis(1, 1, 1);
        for (int it = 0; it < 100; ++it) axis = (cov * axis).normalized();

        // Least-squares polynomial offset from the axis, twice (outlier removal)
        std::vector<bool> use(lumenPts.size(), true);
        real_t tL0 = 0, tL1 = 0;
        gsMatrix<> coef;
        for (int pass = 0; pass < 2; ++pass)
        {
            tL0 = std::numeric_limits<real_t>::max(); tL1 = -tL0;
            for (size_t i = 0; i < lumenPts.size(); ++i)
                if (use[i]) { const real_t t = axis.dot(lumenPts[i] - m0); tL0 = math::min(tL0, t); tL1 = math::max(tL1, t); }
            index_t nUse = 0;
            for (size_t i = 0; i < use.size(); ++i) nUse += use[i];
            gsMatrix<> A(nUse, clDegree + 1), Y(nUse, 3);
            for (size_t i = 0, row = 0; i < lumenPts.size(); ++i)
            {
                if (!use[i]) continue;
                const Vec3 q = lumenPts[i] - m0;
                const real_t t = axis.dot(q);
                const real_t tau = 2 * (t - tL0) / (tL1 - tL0) - 1;
                for (index_t k = 0; k <= clDegree; ++k) A(row, k) = math::pow(tau, (real_t)k);
                Y.row(row) = (q - t * axis).transpose();
                ++row;
            }
            coef = (A.transpose() * A).ldlt().solve(A.transpose() * Y);
            gsMatrix<> res = (A * coef - Y).rowwise().norm();
            std::vector<real_t> rs(res.data(), res.data() + res.size());
            std::nth_element(rs.begin(), rs.begin() + rs.size() / 2, rs.end());
            const real_t cut = 3 * rs[rs.size() / 2] + 1e-12;
            for (size_t i = 0, row = 0; i < lumenPts.size(); ++i)
                if (use[i]) use[i] = (res(row++) <= cut);
        }
        auto curvePoint = [&](real_t t) -> Vec3
        {
            const real_t tc = math::min(tL1, math::max(tL0, t));
            const real_t tau = 2 * (tc - tL0) / (tL1 - tL0) - 1, dtau = 2 / (tL1 - tL0);
            Vec3 q = Vec3::Zero(), dq = Vec3::Zero();
            for (index_t k = 0; k <= clDegree; ++k)
            {
                q += math::pow(tau, (real_t)k) * coef.row(k).transpose();
                if (k > 0) dq += k * math::pow(tau, (real_t)(k - 1)) * dtau * coef.row(k).transpose();
            }
            const Vec3 p = m0 + tc * axis + q, tan = axis + dq;
            return p + (t - tc) * tan;
        };
        real_t tV0 = std::numeric_limits<real_t>::max(), tV1 = -tV0;
        for (size_t i = 0; i < tris.size(); ++i)
        {
            const real_t t = axis.dot(tris[i] - m0);
            tV0 = math::min(tV0, t); tV1 = math::max(tV1, t);
        }
        const real_t margin = 0.25 * (tV1 - tV0);
        const size_t K = 512;
        cl = Centerline();
        for (size_t k = 0; k < K; ++k)
            cl.P.push_back(curvePoint(tV0 - margin + (tV1 - tV0 + 2 * margin) * k / (K - 1)));
        cl.s.assign(K, 0);
        for (size_t k = 1; k < K; ++k) cl.s[k] = cl.s[k - 1] + (cl.P[k] - cl.P[k - 1]).norm();
        cl.T.resize(K); cl.N.resize(K); cl.B.resize(K);
        for (size_t k = 0; k < K; ++k)
            cl.T[k] = (cl.P[std::min(k + 1, K - 1)] - cl.P[k > 0 ? k - 1 : 0]).normalized();
        {
            Vec3 e = Vec3::Zero();
            index_t im;
            cl.T[0].cwiseAbs().minCoeff(&im);
            e[im] = 1;
            cl.N[0] = (e - e.dot(cl.T[0]) * cl.T[0]).normalized();
            for (size_t k = 1; k < K; ++k)
                cl.N[k] = (cl.N[k - 1] - cl.N[k - 1].dot(cl.T[k]) * cl.T[k]).normalized();
            for (size_t k = 0; k < K; ++k) cl.B[k] = cl.T[k].cross(cl.N[k]);
        }
        kappaMax = 0;
        for (size_t k = 1; k + 1 < K; ++k)
            kappaMax = math::max(kappaMax, (cl.T[k + 1] - cl.T[k - 1]).norm() / (cl.s[k + 1] - cl.s[k - 1]));

        vs.resize(tris.size()); vr.resize(tris.size()); vphi.resize(tris.size());
        rMax = 0;
        for (size_t i = 0; i < tris.size(); ++i)
        {
            cl.project(tris[i], vs[i], vr[i], vphi[i]);
            rMax = math::max(rMax, vr[i]);
        }
        gsInfo << "Urethra centreline: " << lumenPts.size() << " lumen points, length " << cl.s.back()
               << " mm, max r*curvature " << kappaMax * rMax << "\n";
        if (kappaMax * rMax > 1)
            gsWarn << "The urethra centreline bends too much for the prostate size (max r * curvature = "
                   << kappaMax * rMax << "); the volume may fold. Try --clDegree 2.\n";
    }

    Vec3 map(const Vec2 &sr, real_t phi) const { return cl.map(sr(0), sr(1), phi); }

    /// Cut of the surface with the half-plane {phi = const}: segments in (s, r)
    Segments cut(real_t phi) const
    {
        Segments seg;
        for (size_t t = 0; t < m.nTri(); ++t)
        {
            real_t dl[3];
            bool skip = false;
            for (int k = 0; k < 3; ++k)
            {
                dl[k] = wrapPi(vphi[3 * t + k] - phi);
                skip |= math::abs(dl[k]) > 0.5 * EIGEN_PI; // other side of the axis
            }
            if (skip) continue;
            gsVector<real_t, 4> sg;
            int np = 0;
            for (int k = 0; k < 3 && np < 2; ++k)
            {
                const int k1 = (k + 1) % 3;
                if ((dl[k] >= 0) == (dl[k1] >= 0)) continue;
                const real_t a = dl[k] / (dl[k] - dl[k1]);
                sg(2 * np)     = (1 - a) * vs[3 * t + k] + a * vs[3 * t + k1];
                sg(2 * np + 1) = (1 - a) * vr[3 * t + k] + a * vr[3 * t + k1];
                ++np;
            }
            if (np == 2) seg.push_back(sg);
        }
        GISMO_ENSURE(seg.size() >= 3, "Empty cross-section at phi=" << phi);
        return seg;
    }

    /// Inside test (odd number of crossings in +s) and distance to the boundary
    static bool inside(const Segments &seg, const Vec2 &p, real_t &dmin)
    {
        int crossings = 0;
        dmin = std::numeric_limits<real_t>::max();
        for (size_t k = 0; k < seg.size(); ++k)
        {
            const real_t s0 = seg[k](0), r0 = seg[k](1), s1 = seg[k](2), r1 = seg[k](3);
            if ((r0 > p(1)) != (r1 > p(1)) && p(0) < s0 + (p(1) - r0) * (s1 - s0) / (r1 - r0)) ++crossings;
            const real_t ex = s1 - s0, ey = r1 - r0, l2 = ex * ex + ey * ey;
            const real_t h = l2 > 0 ? math::min((real_t)1, math::max((real_t)0, ((p(0) - s0) * ex + (p(1) - r0) * ey) / l2)) : 0;
            dmin = math::min(dmin, math::sqrt(math::pow(p(0) - s0 - h * ex, 2) + math::pow(p(1) - r0 - h * ey, 2)));
        }
        return crossings % 2 == 1;
    }

    /// Deep interior point: max distance to the boundary over a grid of inside points
    static Vec2 deepPoint(const Segments &seg, index_t G = 48)
    {
        Vec2 lo, hi, best;
        lo << seg[0](0), seg[0](1); hi = lo;
        for (size_t k = 0; k < seg.size(); ++k)
            for (int e = 0; e < 2; ++e)
            {
                lo(0) = math::min(lo(0), seg[k](2 * e)); hi(0) = math::max(hi(0), seg[k](2 * e));
                lo(1) = math::min(lo(1), seg[k](2 * e + 1)); hi(1) = math::max(hi(1), seg[k](2 * e + 1));
            }
        real_t bd = -1, d;
        for (index_t a = 1; a < G; ++a)
            for (index_t b = 1; b < G; ++b)
            {
                Vec2 p;
                p << lo(0) + (hi(0) - lo(0)) * a / G, lo(1) + (hi(1) - lo(1)) * b / G;
                if (inside(seg, p, d) && d > bd) { bd = d; best = p; }
            }
        GISMO_ENSURE(bd > 0, "No interior point found in a cross-section.");
        return best;
    }

    /// Nearest boundary crossing of the ray c + t (cos th, sin th); hits = #crossings
    static Vec2 ray(const Segments &seg, const Vec2 &c, real_t th, int &hits)
    {
        const real_t dx = math::cos(th), dy = math::sin(th);
        real_t tNear = -1;
        hits = 0;
        for (size_t k = 0; k < seg.size(); ++k)
        {
            const real_t ex = seg[k](2) - seg[k](0), ey = seg[k](3) - seg[k](1);
            const real_t den = dx * ey - dy * ex;
            if (math::abs(den) < 1e-14) continue;
            const real_t wx = seg[k](0) - c(0), wy = seg[k](1) - c(1);
            const real_t t = (wx * ey - wy * ex) / den, h = (wx * dy - wy * dx) / den;
            if (t <= 0 || h < -1e-9 || h > 1 + 1e-9) continue;
            ++hits;
            if (tNear < 0 || t < tNear) tNear = t;
        }
        GISMO_ENSURE(tNear > 0, "Cross-section is not closed.");
        Vec2 b;
        b << c(0) + tNear * dx, c(1) + tNear * dy;
        return b;
    }
};

// =============================================================================
// TENSOR-PRODUCT VOLUMES
// =============================================================================
/// Reverses the v-direction of the coefficients (orientation fix, keeps the geometry)
inline void flipV(gsMatrix<> &C, index_t nu, index_t nv, index_t nw)
{
    gsMatrix<> F = C;
    for (index_t k = 0; k < nw; ++k)
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
                F.row(i + nu * (j + nv * k)) = C.row(i + nu * ((nv - 1 - j) + nv * k));
    C.swap(F);
}

inline real_t detAtMid(const gsGeometry<> &g)
{
    gsMatrix<> J, mid(3, 1);
    mid << 0.5, 0.5, 0.75;
    g.jacobian_into(mid, J);
    return J.determinant();
}

/// Legacy construction (as in torus_prostate_parametrization_example.cpp):
/// interpolation at the Greville points, raw cross-section centres.
inline gsGeometry<>::uPtr buildTensorInterp(const ProstateSurface &S, index_t degree, index_t nCirc,
                                           index_t nRad, index_t &nonStar)
{
    gsKnotVector<> kc(0, 1, nCirc - 1, degree + 1), kr(0, 1, nRad - 1, degree + 1);
    const gsMatrix<> gC = kc.greville(), gR = kr.greville();
    const index_t nu = gC.size(), nw = gR.size();
    std::vector<Vec2> cen(nu);
    std::vector<std::vector<Vec2> > bnd(nu, std::vector<Vec2>(nu));
    nonStar = 0;
    for (index_t i = 0; i + 1 < nu; ++i)
    {
        const Segments seg = S.cut(2 * EIGEN_PI * gC(i));
        cen[i] = ProstateSurface::deepPoint(seg);
        for (index_t j = 0; j + 1 < nu; ++j)
        {
            int hits;
            bnd[i][j] = ProstateSurface::ray(seg, cen[i], 2 * EIGEN_PI * gC(j), hits);
            nonStar += (hits > 2);
        }
        bnd[i][nu - 1] = bnd[i][0];
    }
    cen[nu - 1] = cen[0];
    bnd[nu - 1] = bnd[0];

    gsTensorBSplineBasis<3, real_t> basis(kc, kc, kr);
    gsMatrix<> vals(3, basis.size());
    for (index_t k = 0; k < nw; ++k)
        for (index_t j = 0; j < nu; ++j)
            for (index_t i = 0; i < nu; ++i)
            {
                const Vec2 x = cen[i] + gR(k) * (bnd[i][j] - cen[i]);
                const real_t phi = (i + 1 == nu) ? 0 : 2 * EIGEN_PI * gC(i);
                vals.col(i + nu * (j + nu * k)) = S.map(x, phi);
            }
    gsGeometry<>::uPtr vol = basis.interpolateAtAnchors(vals);
    if (detAtMid(*vol) < 0)
    {
        gsMatrix<> C = vol->coefs();
        flipV(C, nu, nu, nw);
        vol = basis.makeGeometry(give(C));
    }
    return vol;
}

/// Constraint matrix of a clamped uniform 1D basis (n functions, degree p) made
/// periodic: C1 (p >= 2): c_0 = c_{n-1} = (c_1 + c_{n-2})/2, free c_1..c_{n-2};
/// C0: c_{n-1} = c_0, free c_0..c_{n-2}.
inline gsMatrix<> periodicZ(index_t n, index_t p, bool C1)
{
    if (C1 && p >= 2)
    {
        gsMatrix<> Z = gsMatrix<>::Zero(n, n - 2);
        for (index_t r = 1; r < n - 1; ++r) Z(r, r - 1) = 1;
        Z(0, 0) = Z(0, n - 3) = Z(n - 1, 0) = Z(n - 1, n - 3) = 0.5;
        return Z;
    }
    gsMatrix<> Z = gsMatrix<>::Zero(n, n - 1);
    for (index_t r = 0; r < n - 1; ++r) Z(r, r) = 1;
    Z(n - 1, 0) = 1;
    return Z;
}

/// Sparse Z_u (x) Z_v (x) Z_w for the flat tensor index i + nu (j + nv k).
/// Z_w drops the first layer (k = 0: the collapsed core is not changed).
inline gsSparseMatrix<> tensorZ(const gsTensorBSplineBasis<3, real_t> &b, bool C1)
{
    const gsMatrix<> Zu = periodicZ(b.size(0), b.degree(0), C1), Zv = periodicZ(b.size(1), b.degree(1), C1);
    const index_t nu = b.size(0), nv = b.size(1), nw = b.size(2);
    const index_t fu = Zu.cols(), fv = Zv.cols();
    gsSparseEntries<> E;
    for (index_t k = 1; k < nw; ++k)
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
                for (index_t a = 0; a < fu; ++a)
                {
                    if (Zu(i, a) == 0) continue;
                    for (index_t bb = 0; bb < fv; ++bb)
                        if (Zv(j, bb) != 0)
                            E.add(i + nu * (j + nv * k), a + fu * (bb + fv * (k - 1)), Zu(i, a) * Zv(j, bb));
                }
    gsSparseMatrix<> Z(b.size(), fu * fv * (nw - 1));
    Z.setFrom(E);
    return Z;
}

/// Mode-d product of a 3-tensor (dims d[0..2], flat index i0 + d0 (i1 + d1 i2),
/// columns = components) with the matrix M (new_d x d[mode]).
inline gsMatrix<> modeProduct(const gsMatrix<> &X, index_t d[3], int mode, const gsMatrix<> &M)
{
    index_t e[3] = {d[0], d[1], d[2]};
    e[mode] = M.rows();
    gsMatrix<> Y = gsMatrix<>::Zero(e[0] * e[1] * e[2], X.cols());
    for (index_t i2 = 0; i2 < d[2]; ++i2)
        for (index_t i1 = 0; i1 < d[1]; ++i1)
            for (index_t i0 = 0; i0 < d[0]; ++i0)
            {
                const index_t src = i0 + d[0] * (i1 + d[1] * i2);
                const index_t ii[3] = {i0, i1, i2};
                for (index_t a = 0; a < M.rows(); ++a)
                {
                    const real_t m = M(a, ii[mode]);
                    if (m == 0) continue;
                    index_t jj[3] = {i0, i1, i2};
                    jj[mode] = a;
                    Y.row(jj[0] + e[0] * (jj[1] + e[1] * jj[2])) += m * X.row(src);
                }
            }
    for (int k = 0; k < 3; ++k) d[k] = e[k];
    return Y;
}

/// Collocation matrix (m x n) of a 1D B-spline basis at the points t
inline gsMatrix<> collocation(const gsBSplineBasis<> &b, const gsMatrix<> &t)
{
    gsMatrix<> vals;
    gsMatrix<index_t> act;
    b.eval_into(t, vals);
    b.active_into(t, act);
    gsMatrix<> B = gsMatrix<>::Zero(t.cols(), b.size());
    for (index_t c = 0; c < t.cols(); ++c)
        for (index_t r = 0; r < act.rows(); ++r) B(c, act(r, c)) += vals(r, c);
    return B;
}

/// Least-squares pseudo-inverse of A = B Z (rows of the result: free coefficients)
inline gsMatrix<> lsOperator(const gsMatrix<> &B, const gsMatrix<> &Z)
{
    const gsMatrix<> A = B * Z;
    const gsMatrix<> AtA = A.transpose() * A;
    return AtA.ldlt().solve(A.transpose());
}

struct lsOptions
{
    index_t samples      = 4;    // samples per element in each direction
    index_t centreModes  = 4;    // Fourier modes of the smoothed centre c(phi) (0: raw centres)
    bool    C1           = true; // C1 across the u- and v-seams (else C0)
    real_t  curvWeight   = 0.0;  // v-parametrization: < 0 polar angle; >= 0 arc length
                                 // + curvWeight * (turning angle / 2 pi) * perimeter
};

/// Boundary point of a section at v in [0,1): polar angle (curvWeight < 0) or
/// curvature-weighted arc length of the contour, starting at the ray theta = 0.
/// The contour is resolved by nRay polar rays. Returns the points for all vs.
inline std::vector<Vec2> sectionBoundary(const Segments &seg, const Vec2 &c, const std::vector<real_t> &vs,
                                         real_t curvWeight, index_t &nonStar, index_t nRay = 1024)
{
    std::vector<Vec2> out(vs.size());
    int hits;
    if (curvWeight < 0)
    {
        for (size_t j = 0; j < vs.size(); ++j)
        {
            out[j] = ProstateSurface::ray(seg, c, 2 * EIGEN_PI * vs[j], hits);
            nonStar += (hits > 2);
        }
        return out;
    }
    std::vector<Vec2> b(nRay + 1);
    for (index_t k = 0; k < nRay; ++k)
    {
        b[k] = ProstateSurface::ray(seg, c, 2 * EIGEN_PI * k / nRay, hits);
        nonStar += (hits > 2) && (k % (nRay / 32) == 0); // count like 32 coarse rays
    }
    b[nRay] = b[0];
    // segment lengths and turning angles at the vertices
    std::vector<real_t> len(nRay), turn(nRay);
    real_t perim = 0;
    for (index_t k = 0; k < nRay; ++k) { len[k] = (b[k + 1] - b[k]).norm(); perim += len[k]; }
    for (index_t k = 0; k < nRay; ++k)
    {
        const Vec2 e0 = b[k] - b[(k + nRay - 1) % nRay], e1 = b[k + 1] - b[k];
        const real_t n0 = e0.norm(), n1 = e1.norm();
        turn[k] = (n0 > 0 && n1 > 0) ? math::abs(math::atan2(e0(0) * e1(1) - e0(1) * e1(0), e0.dot(e1))) : 0;
    }
    // smooth the turning angle (the cut contour is piecewise linear and noisy)
    const index_t win = math::max<index_t>(1, nRay / 64);
    std::vector<real_t> ts(nRay, 0);
    for (index_t k = 0; k < nRay; ++k)
    {
        for (index_t l = -win; l <= win; ++l) ts[k] += turn[(k + l + nRay) % nRay];
        ts[k] /= (2 * win + 1);
    }
    // cumulative measure
    std::vector<real_t> cum(nRay + 1, 0);
    for (index_t k = 0; k < nRay; ++k)
        cum[k + 1] = cum[k] + len[k] + curvWeight * perim / (2 * EIGEN_PI) * 0.5 * (ts[k] + ts[(k + 1) % nRay]);
    for (size_t j = 0; j < vs.size(); ++j)
    {
        const real_t target = vs[j] * cum[nRay];
        const index_t k = math::min<index_t>(nRay - 1, std::upper_bound(cum.begin(), cum.end(), target) - cum.begin() - 1);
        const real_t a = (cum[k + 1] > cum[k]) ? (target - cum[k]) / (cum[k + 1] - cum[k]) : 0;
        out[j] = (1 - a) * b[k] + a * b[k + 1];
    }
    return out;
}

/// Least-squares tensor volume: densely sampled cross-sections, smoothed
/// centres, separable (Kronecker) LS with periodic seam constraints, and an
/// exactly collapsed core line (x = C(u) + D(u,v,w), D(u,v,0) = 0).
inline gsGeometry<>::uPtr buildTensorLS(const ProstateSurface &S, index_t degree, index_t nCirc,
                                       index_t nRad, const lsOptions &o, index_t &nonStar)
{
    gsKnotVector<> kc(0, 1, nCirc - 1, degree + 1), kr(0, 1, nRad - 1, degree + 1);
    gsBSplineBasis<> bc(kc), br(kr);
    const index_t n = bc.size(), nw = br.size();
    const index_t mu = math::max(o.samples * nCirc, n + 2), mw = math::max(o.samples * nRad, nw) + 1;

    gsMatrix<> tu(1, mu), tw(1, mw);
    for (index_t i = 0; i < mu; ++i) tu(0, i) = (i + 0.5) / mu;   // periodic: avoid duplicates
    for (index_t k = 0; k < mw; ++k) tw(0, k) = real_t(k) / (mw - 1);

    // Cross-sections and their deep interior points
    std::vector<Segments> segs(mu);
    std::vector<Vec2> raw(mu), cen(mu);
    for (index_t i = 0; i < mu; ++i)
    {
        segs[i] = S.cut(2 * EIGEN_PI * tu(0, i));
        raw[i] = ProstateSurface::deepPoint(segs[i]);
    }
    // Smoothed centre: least-squares Fourier series in phi (the grid argmax jumps)
    index_t nRaw = 0;
    if (o.centreModes > 0)
    {
        const index_t K = math::min(o.centreModes, (mu - 1) / 2);
        gsMatrix<> F(mu, 2 * K + 1), Y(mu, 2);
        for (index_t i = 0; i < mu; ++i)
        {
            const real_t ph = 2 * EIGEN_PI * tu(0, i);
            F(i, 0) = 1;
            for (index_t k = 1; k <= K; ++k) { F(i, 2 * k - 1) = math::cos(k * ph); F(i, 2 * k) = math::sin(k * ph); }
            Y.row(i) = raw[i].transpose();
        }
        const gsMatrix<> a = (F.transpose() * F).ldlt().solve(F.transpose() * Y);
        const gsMatrix<> Ys = F * a;
        for (index_t i = 0; i < mu; ++i)
        {
            cen[i] = Ys.row(i).transpose();
            real_t d;
            // keep the smoothed centre only if it is well inside the section
            if (!ProstateSurface::inside(segs[i], cen[i], d) || d < 0.25 * S.m.diag / nCirc)
            { cen[i] = raw[i]; ++nRaw; }
        }
    }
    else
        cen = raw;
    if (nRaw > 0)
        gsWarn << nRaw << " of " << mu << " smoothed cross-section centres fell outside/too close "
               << "to the boundary; the raw centre is used there.\n";

    // Samples: core C(u_i) and D(u_i, v_j, w_k) = x - C
    gsMatrix<> Cd(mu, 3), Dd(mu * mu * mw, 3);
    nonStar = 0;
    std::vector<real_t> vs(mu);
    for (index_t j = 0; j < mu; ++j) vs[j] = tu(0, j);
    for (index_t i = 0; i < mu; ++i)
    {
        const real_t phi = 2 * EIGEN_PI * tu(0, i);
        const Vec3 c3 = S.map(cen[i], phi);
        Cd.row(i) = c3.transpose();
        const std::vector<Vec2> bv = sectionBoundary(segs[i], cen[i], vs, o.curvWeight, nonStar);
        for (index_t j = 0; j < mu; ++j)
        {
            const Vec2 &b = bv[j];
            for (index_t k = 0; k < mw; ++k)
                Dd.row(i + mu * (j + mu * k)) = (S.map(cen[i] + tw(0, k) * (b - cen[i]), phi) - c3).transpose();
        }
    }

    const gsMatrix<> Zc = periodicZ(n, degree, o.C1);
    gsMatrix<> Zw = gsMatrix<>::Zero(nw, nw - 1); // drop k = 0: D = 0 on the core
    for (index_t k = 1; k < nw; ++k) Zw(k, k - 1) = 1;
    const gsMatrix<> Bu = collocation(bc, tu), Bw = collocation(br, tw);
    const gsMatrix<> Pu = lsOperator(Bu, Zc), Pw = lsOperator(Bw, Zw);

    // core curve
    const gsMatrix<> Cc = Zc * (Pu * Cd); // n x 3
    // volume offsets
    index_t d[3] = {mu, mu, mw};
    gsMatrix<> T = modeProduct(Dd, d, 0, Pu);
    T = modeProduct(T, d, 1, Pu);
    T = modeProduct(T, d, 2, Pw);
    T = modeProduct(T, d, 0, Zc);
    T = modeProduct(T, d, 1, Zc);
    T = modeProduct(T, d, 2, Zw); // n x n x nw
    for (index_t k = 0; k < nw; ++k)
        for (index_t j = 0; j < n; ++j)
            for (index_t i = 0; i < n; ++i)
                T.row(i + n * (j + n * k)) += Cc.row(i);

    gsTensorBSplineBasis<3, real_t> basis(kc, kc, kr);
    gsGeometry<>::uPtr vol = basis.makeGeometry(give(T));
    if (detAtMid(*vol) < 0)
    {
        gsMatrix<> C = vol->coefs();
        flipV(C, n, n, nw);
        vol = basis.makeGeometry(give(C));
    }
    gsInfo << "LS volume: " << mu << " cross-sections x " << mu << " rays x " << mw
           << " radial samples, centre Fourier modes " << o.centreModes
           << ", seams " << (o.C1 ? "C1" : "C0") << ", v-parametrization "
           << (o.curvWeight < 0 ? std::string("polar angle")
                                : "arc length + " + util::to_string(o.curvWeight) + " x turning") << "\n";
    return vol;
}

/// Adds the seam interfaces and the boundary labels of the torus topology
inline void setTorusTopology(gsMultiPatch<> &mp)
{
    mp.addInterface(0, boundary::west,  0, boundary::east);  // toroidal seam u=0 <-> u=1
    mp.addInterface(0, boundary::south, 0, boundary::north); // poloidal seam v=0 <-> v=1
    mp.addBoundary(0, boundary::front, "core");              // w = 0, collapsed polar axis
    mp.addBoundary(0, boundary::back,  "surface");           // w = 1, capsule + urethra wall
}

/// Indices of the functions that do not vanish on the face w = 0, grouped by
/// their u-anchor (functions of one group sit at the same point of the core line).
inline std::vector<std::vector<index_t> > coreGroups(const gsBasis<> &b)
{
    const gsMatrix<index_t> front = b.boundary(boundary::front);
    gsMatrix<> anc;
    b.anchors_into(anc);
    std::map<long long, std::vector<index_t> > g;
    for (index_t r = 0; r < front.size(); ++r)
    {
        const index_t i = front(r);
        g[static_cast<long long>(math::floor(anc(0, i) * 1e9 + 0.5))].push_back(i);
    }
    std::vector<std::vector<index_t> > out;
    for (auto &kv : g) out.push_back(kv.second);
    return out;
}

/// Dof mapper of a single-patch torus basis: C0 seams (patch self-interfaces);
/// optionally the core functions merged per u-anchor (polar) or eliminated.
inline gsDofMapper torusMapper(const gsMultiBasis<> &mb, bool polarMerge, bool eliminateCore)
{
    gsDofMapper dm(mb);
    for (gsBoxTopology::const_iiterator it = mb.topology().iBegin(); it != mb.topology().iEnd(); ++it)
        mb.matchInterface(*it, dm);
    if (eliminateCore)
        dm.markBoundary(0, mb.basis(0).boundary(boundary::front));
    else if (polarMerge)
    {
        const std::vector<std::vector<index_t> > g = coreGroups(mb.basis(0));
        for (size_t a = 0; a < g.size(); ++a)
            for (size_t k = 1; k < g[a].size(); ++k)
                dm.matchDof(0, g[a][0], 0, g[a][k]);
    }
    dm.finalize();
    return dm;
}

/// Sparse Z (full -> free) from a mapper (eliminated dofs give zero rows)
inline gsSparseMatrix<> mapperZ(const gsDofMapper &dm, index_t size)
{
    gsSparseEntries<> E;
    for (index_t i = 0; i < size; ++i)
        if (dm.is_free(i, 0)) E.add(i, dm.index(i, 0), 1.0);
    gsSparseMatrix<> Z(size, dm.freeSize());
    Z.setFrom(E);
    return Z;
}

// =============================================================================
// QUALITY: scaled Jacobian, surface deviation
// =============================================================================
struct jacobianStats { real_t sjMin, sjMax, detMin, vol; };

inline jacobianStats jacobianReport(const gsMultiPatch<> &mp, const gsMultiBasis<> &mb,
                                    const std::string &label, const std::string &vtkFile = "",
                                    index_t npts = 20000, bool quiet = false)
{
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    gsExprEvaluator<>::geometryMap Gm = ev.getMap(mp);
    auto scaledJac = jac(Gm).det() / (jac(Gm)[0].norm() * jac(Gm)[1].norm() * jac(Gm)[2].norm());
    jacobianStats s;
    s.sjMin  = ev.min(scaledJac);
    s.sjMax  = ev.max(scaledJac);
    s.detMin = ev.min(jac(Gm).det());
    s.vol    = ev.integral(meas(Gm));
    if (!quiet)
        gsInfo << label << ": " << mb.basis(0).numElements() << " elements, scaled Jacobian at the "
               << "quadrature points min " << s.sjMin << ", max " << s.sjMax << ", min det J "
               << s.detMin << (s.sjMin <= 0 ? "  (WARNING: FOLDED, not bijective)" : "  (bijective)")
               << ", volume " << s.vol << " mm^3\n";
    if (!vtkFile.empty())
    {
        ev.options().setInt("plot.npts", npts);
        ev.writeParaview(scaledJac, Gm, vtkFile);
    }
    return s;
}

/// Minimum scaled Jacobian per element of \a b (domain-iterator order, as the
/// gismo marking routines expect), sampled at n^3 points per element INCLUDING
/// the element faces and corners (so folds on the surface w = 1 are seen,
/// unlike at Gauss points). Points on the collapsed face w = 0 are moved
/// 10% into the element (det J = 0 there by construction).
inline std::vector<real_t> elementQuality(const gsGeometry<> &G, const gsBasis<> &b, index_t n = 3)
{
    std::vector<gsVector<real_t, 3> > lo, up;
    typename gsBasis<>::domainIter it = b.domain()->beginAll(), end = b.domain()->endAll();
    for (; it < end; ++it) { lo.push_back(it.lowerCorner()); up.push_back(it.upperCorner()); }
    const index_t nE = lo.size(), nP = n * n * n, chunk = 4096;
    std::vector<real_t> q(nE);
#pragma omp parallel
    {
    gsGeometry<>::uPtr Gt = G.clone(); // per-thread copy (evaluation caches are not thread-safe)
#pragma omp for schedule(dynamic)
    for (index_t e0 = 0; e0 < nE; e0 += chunk)
    {
        const index_t e1 = math::min(nE, e0 + chunk);
        gsMatrix<> P(3, (e1 - e0) * nP), D;
        for (index_t e = e0; e < e1; ++e)
            for (index_t k = 0; k < nP; ++k)
            {
                const index_t t[3] = {k % n, (k / n) % n, k / (n * n)};
                for (int dd = 0; dd < 3; ++dd)
                {
                    real_t s = real_t(t[dd]) / (n - 1);
                    if (dd == 2 && lo[e][2] <= 1e-12 && t[dd] == 0) s = 0.1;
                    P(dd, (e - e0) * nP + k) = lo[e][dd] + s * (up[e][dd] - lo[e][dd]);
                }
            }
        Gt->deriv_into(P, D);
        for (index_t e = e0; e < e1; ++e)
        {
            real_t m = std::numeric_limits<real_t>::max();
            for (index_t k = 0; k < nP; ++k)
            {
                const index_t c = (e - e0) * nP + k;
                gsMatrix<real_t, 3, 3> J;
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s) J(r, s) = D(3 * r + s, c);
                const real_t sj = J.determinant() / (J.col(0).norm() * J.col(1).norm() * J.col(2).norm());
                m = math::min(m, sj);
            }
            q[e] = m;
        }
    }
    }
    return q;
}

inline real_t minOf(const std::vector<real_t> &v) { return *std::min_element(v.begin(), v.end()); }

/// "a: n1, b: n2, ..." counts of the values below the thresholds
inline std::string qualityHistogram(const std::vector<real_t> &q)
{
    const real_t th[5] = {0, 0.1, 0.25, 0.5, 0.75};
    std::ostringstream os;
    for (int k = 0; k < 5; ++k)
        os << (k ? ", " : "") << "<" << th[k] << ": "
           << std::count_if(q.begin(), q.end(), [&](real_t x) { return x < th[k]; });
    return os.str();
}

/// Writes the closed boundary of the volume (the face w = 1; the seams are
/// interior and w = 0 is collapsed) as a watertight OBJ triangle mesh:
/// n x n samples, periodic in u and v (no duplicated seam vertices).
/// Returns the enclosed volume of that triangulation.
inline real_t writeSurfaceObj(const gsGeometry<> &G, const std::string &fn, index_t n = 256)
{
    gsMatrix<> P(3, n * n);
    for (index_t j = 0; j < n; ++j)
        for (index_t i = 0; i < n; ++i)
            P.col(i + n * j) << real_t(i) / n, real_t(j) / n, 1.0;
    const gsMatrix<> X = G.eval(P);
    std::ofstream f(fn.c_str());
    f << std::setprecision(10) << "# fitted prostate surface (w = 1), " << n << " x " << n << " samples\n";
    for (index_t c = 0; c < X.cols(); ++c) f << "v " << X(0, c) << " " << X(1, c) << " " << X(2, c) << "\n";
    real_t vol = 0;
    auto tri = [&](index_t a, index_t b, index_t c)
    {
        f << "f " << a + 1 << " " << b + 1 << " " << c + 1 << "\n";
        vol += Vec3(X.col(a)).dot(Vec3(X.col(b)).cross(Vec3(X.col(c)))) / 6;
    };
    for (index_t j = 0; j < n; ++j)
        for (index_t i = 0; i < n; ++i)
        {
            const index_t a = i + n * j, b = (i + 1) % n + n * j,
                          c = (i + 1) % n + n * ((j + 1) % n), d = i + n * ((j + 1) % n);
            tri(a, b, c);
            tri(a, c, d);
        }
    return math::abs(vol);
}

/// Sample grid of the torus-parametrized volume for ParaView, with the
/// topology of the PHYSICAL domain: the seams u = 0/1 and v = 0/1 share their
/// points and the collapsed face w = 0 is one point per u. Written as an
/// unstructured grid (hexahedra, wedges next to the core line), so ParaView
/// shows only the real boundary w = 1 (capsule + urethra wall). A structured
/// .vts of [0,1]^3 instead draws the seam faces as sheets inside the prostate
/// and the core line as a thin tube.
struct torusPlotGrid
{
    index_t nu, nv, nw;
    gsMatrix<> par; // 3 x nPoints parameter points

    torusPlotGrid(index_t nu_, index_t nv_, index_t nw_) : nu(nu_), nv(nv_), nw(math::max<index_t>(nw_, 2))
    {
        par.resize(3, nu + nu * nv * (nw - 1));
        for (index_t i = 0; i < nu; ++i) par.col(i) << real_t(i) / nu, 0, 0;
        for (index_t k = 1; k < nw; ++k)
            for (index_t j = 0; j < nv; ++j)
                for (index_t i = 0; i < nu; ++i)
                    par.col(id(i, j, k)) << real_t(i) / nu, real_t(j) / nv, real_t(k) / (nw - 1);
    }
    index_t id(index_t i, index_t j, index_t k) const
    {
        i %= nu; j %= nv;
        return k == 0 ? i : nu + i + nu * (j + nv * (k - 1));
    }
    index_t numPoints() const { return par.cols(); }

    /// Writes <fn>.vtu: points X (3 x numPoints) and optional point fields
    void write(const std::string &fn, const gsMatrix<> &X,
               const std::vector<std::pair<std::string, gsMatrix<> > > &fields) const
    {
        std::vector<index_t> conn, offs;
        std::vector<int> types;
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
            {   // wedge between the core line and the first layer
                const index_t w[6] = {id(i, 0, 0), id(i, j, 1), id(i, j + 1, 1),
                                      id(i + 1, 0, 0), id(i + 1, j, 1), id(i + 1, j + 1, 1)};
                conn.insert(conn.end(), w, w + 6); offs.push_back(conn.size()); types.push_back(13);
                for (index_t k = 1; k + 1 < nw; ++k)
                {
                    const index_t h[8] = {id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
                                          id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)};
                    conn.insert(conn.end(), h, h + 8); offs.push_back(conn.size()); types.push_back(12);
                }
            }
        std::ofstream f((fn + ".vtu").c_str());
        f << std::setprecision(9)
          << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" << X.cols() << "\" NumberOfCells=\"" << types.size() << "\">\n";
        if (!fields.empty())
        {
            f << "<PointData Scalars=\"" << fields[0].first << "\">\n";
            for (size_t q = 0; q < fields.size(); ++q)
            {
                f << "<DataArray type=\"Float64\" Name=\"" << fields[q].first << "\" format=\"ascii\">\n";
                for (index_t c = 0; c < fields[q].second.cols(); ++c) f << fields[q].second(0, c) << " ";
                f << "\n</DataArray>\n";
            }
            f << "</PointData>\n";
        }
        f << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (index_t c = 0; c < X.cols(); ++c) f << X(0, c) << " " << X(1, c) << " " << X(2, c) << "\n";
        f << "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t c = 0; c < conn.size(); ++c) f << conn[c] << " ";
        f << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t c = 0; c < offs.size(); ++c) f << offs[c] << " ";
        f << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t c = 0; c < types.size(); ++c) f << types[c] << " ";
        f << "\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }
};

/// The whole volume as <fn>.vtu (merged seams), with the scaled Jacobian as field
inline void writeTorusVolume(const gsGeometry<> &G, const std::string &fn, index_t n = 64)
{
    const torusPlotGrid tg(n, n, n / 2 + 1);
    const gsMatrix<> X = G.eval(tg.par);
    gsMatrix<> P = tg.par, D, sj(1, P.cols());
    for (index_t c = 0; c < P.cols(); ++c) if (P(2, c) == 0) P(2, c) = 1e-3; // core: det J = 0
    G.deriv_into(P, D);
    for (index_t c = 0; c < P.cols(); ++c)
    {
        gsMatrix<real_t, 3, 3> J;
        for (int r = 0; r < 3; ++r) for (int s = 0; s < 3; ++s) J(r, s) = D(3 * r + s, c);
        sj(0, c) = J.determinant() / (J.col(0).norm() * J.col(1).norm() * J.col(2).norm());
    }
    std::vector<std::pair<std::string, gsMatrix<> > > fields;
    fields.push_back(std::make_pair(std::string("scaledJacobian"), sj));
    tg.write(fn, X, fields);
}

/// Smallest scaled Jacobian on a parametric grid (n^3 points, w > 0) and where it is
inline real_t worstScaledJacobian(const gsGeometry<> &G, gsVector<> &where, gsVector<> &x, index_t n = 64)
{
    gsVector<> lo(3), hi(3);
    lo << 0, 0, 0.5 / n; hi << 1, 1, 1;
    gsVector<unsigned> np(3);
    np.setConstant(n);
    const gsMatrix<> P = gsPointGrid<real_t>(lo, hi, np);
    real_t best = std::numeric_limits<real_t>::max();
    index_t ib = 0;
    gsMatrix<> J;
    for (index_t c = 0; c < P.cols(); ++c)
    {
        G.jacobian_into(P.col(c), J);
        const real_t s = J.determinant() / (J.col(0).norm() * J.col(1).norm() * J.col(2).norm());
        if (s < best) { best = s; ib = c; }
    }
    where = P.col(ib);
    x = G.eval(P.col(ib));
    return best;
}

/// Closest points (foot points) on the surface S(u,v) = G(u,v,1) by a nearest
/// grid sample followed by damped Gauss-Newton; u and v are periodic.
inline void footPoints(const gsGeometry<> &G, const std::vector<Vec3> &pts, gsMatrix<> &uv,
                       std::vector<real_t> &dist, index_t nGrid = 160)
{
    gsVector<> lo(2), hi(2);
    lo << 0, 0; hi << 1, 1;
    gsVector<unsigned> np(2);
    np.setConstant(nGrid);
    const gsMatrix<> g2 = gsPointGrid<real_t>(lo, hi, np);
    gsMatrix<> g3(3, g2.cols());
    g3.topRows(2) = g2;
    g3.row(2).setOnes();
    const gsMatrix<> X = G.eval(g3);

    const index_t N = pts.size();
    uv.resize(2, N);
    dist.resize(N);
#pragma omp parallel
    {
    gsGeometry<>::uPtr Gt = G.clone(); // per-thread copy
    const gsGeometry<> &G = *Gt;
#pragma omp for schedule(dynamic, 16)
    for (index_t p = 0; p < N; ++p)
    {
        index_t ib = 0;
        real_t db = std::numeric_limits<real_t>::max();
        for (index_t c = 0; c < X.cols(); ++c)
        {
            const real_t d = (X.col(c) - pts[p]).squaredNorm();
            if (d < db) { db = d; ib = c; }
        }
        gsMatrix<> q(3, 1), x, J;
        q << g2(0, ib), g2(1, ib), 1;
        for (int it = 0; it < 12; ++it)
        {
            G.eval_into(q, x);
            G.jacobian_into(q, J);
            const gsMatrix<> Js = J.leftCols(2);
            const gsVector<> r = pts[p] - x.col(0);
            gsMatrix<> H = Js.transpose() * Js;
            H(0, 0) *= 1.0 + 1e-6; H(1, 1) *= 1.0 + 1e-6;
            gsVector<> du = H.ldlt().solve(Js.transpose() * r);
            const real_t step = du.norm();
            if (step > 0.05) du *= 0.05 / step; // trust region in parameter space
            gsMatrix<> qn = q;
            qn(0, 0) = wrap01(q(0, 0) + du(0));
            qn(1, 0) = wrap01(q(1, 0) + du(1));
            gsMatrix<> xn;
            G.eval_into(qn, xn);
            if ((xn.col(0) - pts[p]).squaredNorm() <= r.squaredNorm()) q = qn;
            else break;
            if (step < 1e-10) break;
        }
        G.eval_into(q, x);
        uv(0, p) = q(0, 0);
        uv(1, p) = q(1, 0);
        dist[p] = (x.col(0) - pts[p]).norm();
    }
    }
}

struct deviationStats
{
    real_t fwdMean = 0, fwdRms = 0, fwdMax = 0; // spline samples -> triangulation
    real_t bwdMean = 0, bwdRms = 0, bwdMax = 0; // triangulation vertices -> spline
    real_t hausdorff() const { return math::max(fwdMax, bwdMax); }
};

/// Two-sided distance between the w = 1 face of G and the triangulation.
/// Optionally writes <prefix>_spline.vts (signed distance on the spline
/// surface) and <prefix>_input.vtp (distance of the input vertices).
inline deviationStats surfaceDeviation(const gsGeometry<> &G, const TriMesh &m, const std::string &label,
                                       const std::string &prefix = "", index_t nGrid = 200)
{
    deviationStats s;
    gsVector<> lo(2), hi(2);
    lo << 0, 0; hi << 1, 1;
    gsVector<unsigned> np(2);
    np.setConstant(nGrid);
    const gsMatrix<> g2 = gsPointGrid<real_t>(lo, hi, np);
    gsMatrix<> g3(3, g2.cols());
    g3.topRows(2) = g2;
    g3.row(2).setOnes();
    const gsMatrix<> X = G.eval(g3);
    std::vector<real_t> sd(X.cols());
#pragma omp parallel for
    for (index_t c = 0; c < X.cols(); ++c)
    {
        Vec3 q;
        sd[c] = m.closest(X.col(c), q);
    }
    for (size_t c = 0; c < sd.size(); ++c)
    {
        const real_t d = math::abs(sd[c]);
        s.fwdMean += d; s.fwdRms += d * d; s.fwdMax = math::max(s.fwdMax, d);
    }
    s.fwdMean /= sd.size(); s.fwdRms = math::sqrt(s.fwdRms / sd.size());

    gsMatrix<> uv;
    std::vector<real_t> bd;
    footPoints(G, m.verts, uv, bd);
    for (size_t p = 0; p < bd.size(); ++p)
    {
        s.bwdMean += bd[p]; s.bwdRms += bd[p] * bd[p]; s.bwdMax = math::max(s.bwdMax, bd[p]);
    }
    s.bwdMean /= bd.size(); s.bwdRms = math::sqrt(s.bwdRms / bd.size());

    gsInfo << label << ": surface -> triangulation mean " << s.fwdMean << ", RMS " << s.fwdRms
           << ", max " << s.fwdMax << " mm; triangulation -> surface mean " << s.bwdMean
           << ", RMS " << s.bwdRms << ", max " << s.bwdMax << " mm; Hausdorff " << s.hausdorff() << " mm\n";

    if (!prefix.empty())
    {
        std::ofstream f((prefix + "_spline.vts").c_str());
        f << std::setprecision(10) << "<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n"
          << "<StructuredGrid WholeExtent=\"0 " << nGrid - 1 << " 0 " << nGrid - 1 << " 0 0\">\n"
          << "<Piece Extent=\"0 " << nGrid - 1 << " 0 " << nGrid - 1 << " 0 0\">\n"
          << "<PointData Scalars=\"signed_distance\">\n<DataArray type=\"Float64\" Name=\"signed_distance\" format=\"ascii\">\n";
        for (size_t c = 0; c < sd.size(); ++c) f << sd[c] << " ";
        f << "\n</DataArray>\n</PointData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (index_t c = 0; c < X.cols(); ++c) f << X(0, c) << " " << X(1, c) << " " << X(2, c) << " ";
        f << "\n</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        m.writeVtp(prefix + "_input", bd, "distance_to_spline");
    }
    return s;
}

// =============================================================================
// CLOSEST-POINT SURFACE CORRECTION
// =============================================================================
struct correctionOptions
{
    index_t its      = 3;     // closest-point iterations
    real_t  lambda   = 1e-3;  // relative weight of the parametric Dirichlet energy
    index_t nFwd     = 3;     // forward samples per element and direction on w = 1
    bool    verbose  = true;
};

/// Parametric Laplacian (stiffness) of a basis on [0,1]^3, identity numbering
inline gsSparseMatrix<> parametricLaplacian(const gsBasis<> &b)
{
    gsMultiBasis<> mb(b);
    gsExprAssembler<> A(1, 1);
    A.setIntegrationElements(mb);
    gsExprAssembler<>::space u = A.getSpace(mb);
    gsDofMapper dm(mb);
    dm.finalize();
    u.setupMapper(dm);
    A.initSystem();
    A.assemble(grad(u) * grad(u).tr());
    return A.matrix();
}

/// Minimizes, over the coefficient update d = Z y,
///   sum_p |G(u_p,v_p,1) + d(u_p,v_p,1) - y_p|^2  +  sum_q |G(q,1) + d(q,1) - z_q|^2
///   + lambda * a(d,d),
/// y_p: triangulation vertices with their foot points (u_p,v_p) on the spline
/// surface, z_q: closest points on the triangulation of spline samples q, a(.,.)
/// the parametric Dirichlet energy over the volume (smooths the correction and
/// extends it harmonically into the interior). Z encodes the seam constraints
/// and fixes the core layer w = 0. The update is halved until the volume is
/// unfolded. Returns false if no unfolded update was found in an iteration.
template <class Check>
bool correctSurface(gsGeometry<> &G, const gsSparseMatrix<> &Z, const TriMesh &m,
                    const correctionOptions &o, Check isValid)
{
    const gsBasis<> &basis = G.basis();
    const gsSparseMatrix<> K = parametricLaplacian(basis);
    // forward sample parameters: nFwd x nFwd points in every element of the
    // surface w = 1 (so locally refined functions always have data)
    gsMatrix<> g2;
    {
        gsBasis<>::uPtr bb = G.basis().boundaryBasis(boundary::back); // (u,v)
        std::vector<Vec2> pts;
        typename gsBasis<>::domainIter it = bb->domain()->beginAll(), end = bb->domain()->endAll();
        for (; it < end; ++it)
        {
            const gsVector<> lo = it.lowerCorner(), up = it.upperCorner();
            for (index_t a = 0; a < o.nFwd; ++a)
                for (index_t b = 0; b < o.nFwd; ++b)
                {
                    Vec2 p;
                    p << lo(0) + (up(0) - lo(0)) * (a + 0.5) / o.nFwd, lo(1) + (up(1) - lo(1)) * (b + 0.5) / o.nFwd;
                    pts.push_back(p);
                }
        }
        g2.resize(2, pts.size());
        for (size_t k = 0; k < pts.size(); ++k) g2.col(k) = pts[k];
    }
    bool allOk = true;
    for (index_t it = 0; it < o.its; ++it)
    {
        // backward data: input vertices at their foot points
        gsMatrix<> uv;
        std::vector<real_t> bd;
        footPoints(G, m.verts, uv, bd);
        // forward data: spline samples -> closest points on the triangulation

        const index_t nB = uv.cols(), nF = g2.cols();
        gsMatrix<> P(3, nB + nF), Y(nB + nF, 3);
        P.topLeftCorner(2, nB) = uv;
        P.topRightCorner(2, nF) = g2;
        P.row(2).setOnes();
        for (index_t p = 0; p < nB; ++p) Y.row(p) = m.verts[p].transpose();
        const gsMatrix<> X = G.eval(P);
#pragma omp parallel for
        for (index_t q = 0; q < nF; ++q)
        {
            Vec3 c;
            m.closest(X.col(nB + q), c);
            Y.row(nB + q) = c.transpose();
        }
        const gsMatrix<> R = Y - X.transpose(); // residuals

        gsMatrix<> vals;
        gsMatrix<index_t> act;
        basis.eval_into(P, vals);
        basis.active_into(P, act);
        gsSparseEntries<> E;
        for (index_t c = 0; c < P.cols(); ++c)
            for (index_t r = 0; r < act.rows(); ++r)
                if (vals(r, c) != 0) E.add(c, act(r, c), vals(r, c));
        gsSparseMatrix<> B(P.cols(), basis.size());
        B.setFrom(E);

        const gsSparseMatrix<> BtB = B.transpose() * B;
        const real_t sB = BtB.diagonal().sum() / BtB.rows(), sK = K.diagonal().sum() / K.rows();
        const real_t lam = o.lambda * sB / sK;
        const gsSparseMatrix<> M = Z.transpose() * (BtB + lam * K) * Z;
        const gsMatrix<> rhs = Z.transpose() * (B.transpose() * R);
        gsSparseSolver<>::SimplicialLDLT solver(M);
        GISMO_ENSURE(solver.info() == gsEigen::Success, "Surface correction: factorization failed.");
        const gsMatrix<> dC = Z * solver.solve(rhs);

        const gsMatrix<> C0 = G.coefs();
        real_t a = 1;
        bool ok = false;
        for (int bt = 0; bt < 6; ++bt, a *= 0.5)
        {
            G.coefs() = C0 + a * dC;
            if (isValid(G)) { ok = true; break; }
        }
        if (!ok) { G.coefs() = C0; allOk = false; }
        if (o.verbose)
        {
            real_t mb = 0, xb = 0;
            for (size_t p = 0; p < bd.size(); ++p) { mb += bd[p]; xb = math::max(xb, bd[p]); }
            gsInfo << "  correction it " << it + 1 << ": before it, vertices -> surface mean "
                   << mb / bd.size() << " max " << xb << " mm; |dC|_max " << dC.rowwise().norm().maxCoeff()
                   << " mm, step " << (ok ? a : 0) << (ok ? "" : " (rejected: would fold)") << "\n";
        }
        if (!ok) break;
    }
    return allOk;
}

// =============================================================================
// FISHER-KOLMOGOROV SOLVER
// =============================================================================
template <class Solver, class T>
index_t fkSolveCG(Solver &solver, const char *name, const gsSparseMatrix<T> &K, const gsMatrix<T> &b,
                  gsMatrix<T> &x, index_t maxIt, T tol)
{
    solver.setMaxIterations(maxIt);
    solver.setTolerance(tol);
    solver.compute(K);
    GISMO_ENSURE(solver.info() == gsEigen::Success, name << ": the preconditioner setup failed.");
    x = solver.solve(b);
    GISMO_ENSURE(solver.info() == gsEigen::Success,
                 name << " did not converge: " << solver.iterations() << "/" << maxIt << " iteration(s), "
                 "error = " << solver.error() << " > " << tol << ". Raise --cgMaxIt or --cgTol, or use "
                 "--solver cgic / direct.");
    return solver.iterations();
}

template <class T>
index_t fkSolveLinear(const std::string &linSolver, const gsSparseMatrix<T> &K, const gsMatrix<T> &b,
                      gsMatrix<T> &x, index_t maxIt, T tol)
{
    typedef gsEigen::SparseMatrix<T, 0, index_t> EigenMat;
    if (linSolver == "cgic")
    {
        gsEigen::ConjugateGradient<EigenMat, gsEigen::Lower | gsEigen::Upper,
                                   gsEigen::IncompleteCholesky<T, gsEigen::Lower, gsEigen::AMDOrdering<index_t> > > solver;
        return fkSolveCG(solver, "CG (incomplete Cholesky)", K, b, x, maxIt, tol);
    }
    else if (linSolver == "cgdiag")
    {
        gsEigen::ConjugateGradient<EigenMat, gsEigen::Lower | gsEigen::Upper,
                                   gsEigen::DiagonalPreconditioner<T> > solver;
        return fkSolveCG(solver, "CG (diagonal)", K, b, x, maxIt, tol);
    }
    GISMO_ENSURE(linSolver == "direct", "Unknown linear solver \"" << linSolver << "\" (cgdiag|cgic|direct)");
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::SimplicialLDLT solver;
#endif
    solver.compute(K);
    GISMO_ENSURE(solver.info() == gsEigen::Success, "Factorisation of the FK matrix failed.");
    x = solver.solve(b);
    return 1;
}

// Model parameters, PHYSICAL units (mm, day)
template <class T>
struct fkParams
{
    // defaults = the torus paper runs (run_fk_torus_params.sh, set 'median')
    T D     = 1.26e-3; // tumor cell diffusivity         [mm^2/day]
    T rho   = 1.86e-3; // net proliferation rate         [1/day]
    T R     = 2.5;     // radius of the initial tumor    [mm]
    T sigma = 1.0;     // interface width of the IC      [mm]
    T amp   = 0.5;    // peak density of the IC, N0 <= amp <= 1 [-]
    T Nth   = 0.15;   // threshold defining Omega_T      [-]
    T xc = std::numeric_limits<T>::quiet_NaN(); // seed centre [mm]; NaN: on the core line
    T yc = std::numeric_limits<T>::quiet_NaN();
    T zc = std::numeric_limits<T>::quiet_NaN();
};

template <class T>
gsFunctionExpr<T> tumorSeed(const fkParams<T> &par)
{
    std::ostringstream expr;
    expr.precision(16);
    expr << par.amp << "*0.5*(1-tanh((sqrt((x-(" << par.xc << "))^2+(y-(" << par.yc
         << "))^2+(z-(" << par.zc << "))^2)-(" << par.R << "))/(" << par.sigma << ")))";
    gsInfo << "Initial condition: N(x,0) = " << expr.str() << "\n";
    return gsFunctionExpr<T>(expr.str(), 3);
}

inline std::set<index_t> fkSaveSteps(const std::string &s, real_t dt)
{
    std::set<index_t> steps;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (tok.find_first_not_of(" \t") == std::string::npos) continue;
        const real_t t = std::stod(tok);
        const index_t n = static_cast<index_t>(math::floor(t / dt + 0.5));
        if (math::abs(n * dt - t) > 1e-8 * math::max((real_t)1, t))
            gsWarn << "--saveAt " << t << " is not a multiple of dt; saving at t = " << n * dt << " instead.\n";
        steps.insert(n);
    }
    return steps;
}

template <class T, class Expr>
T computeInt(const gsMultiBasis<T> &integrationBasis, const gsFunctionSet<T> &geometryMap,
             const Expr &f, const gsOptionList &options)
{
    gsExprAssembler<T> A(1, 1);
    A.options().update(options, gsOptionList::addIfUnknown);
    A.setIntegrationElements(integrationBasis);
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);
    gsExprEvaluator<T> ev(A);
    ev.options().update(A.options(), gsOptionList::addIfUnknown);
    ev.options().setSwitch("SameElement", false);
    return ev.integral(f * meas(G));
}

struct fkSolverSettings
{
    real_t  dt        = 1.0;   // [day]   (torus paper runs: dt = 1 day,
    index_t maxSteps  = 730;   //          T = 730 days)
    std::string out   = "output";
    std::string method = "bodyfitted_prostate";
    bool    plot      = false;
    index_t plotPts   = 50000;
    bool    polarMatch= true;   // merge the DoFs of the collapsed face w = 0
    std::string linSolver = "cgdiag";
    index_t cgMaxIt   = 20000;
    real_t  cgTol     = 1e-8;
    std::string saveAt;
};

/// Generalized-alpha / Newton Fisher-Kolmogorov solver on the map mp with the
/// solution basis dbasis (tensor or THB). Seams glued C0, polar axis merged.
template <class T>
void solveFK(const gsMultiPatch<T> &mp, const gsMultiBasis<T> &dbasis, const gsOptionList &TIMEopt,
             const gsOptionList &Aopt, const fkParams<T> &par, const fkSolverSettings &st,
             index_t numRefine, real_t volRef)
{
    const real_t dt = st.dt;
    const index_t maxSteps = st.maxSteps;
    const std::string &out = st.out;
    const std::set<index_t> saveSteps = fkSaveSteps(st.saveAt, dt);
    gsStopwatch swTotal;
    const gsBasis<T> &sb = dbasis.basis(0);

    // 1. Element sizes [mm] (longest mapped edge chord per element)
    real_t hMax = 0, hMin = std::numeric_limits<real_t>::max(), hSum = 0, hTorMax = 0, hCrossMax = 0;
    {
        typename gsBasis<T>::domainIter domIt = sb.domain()->beginAll(), domItEnd = sb.domain()->endAll();
        gsMatrix<T> c(3, 8), x;
        index_t nE = 0;
        for (; domIt < domItEnd; ++domIt)
        {
            const gsVector<T> lo = domIt.lowerCorner(), up = domIt.upperCorner();
            for (index_t k = 0; k != 8; ++k)
                for (index_t d = 0; d != 3; ++d) c(d, k) = ((k >> d) & 1) ? up[d] : lo[d];
            mp.patch(0).eval_into(c, x);
            real_t he = 0;
            for (index_t k = 0; k != 8; ++k)
                for (index_t d = 0; d != 3; ++d)
                    if (!((k >> d) & 1))
                    {
                        const real_t len = (x.col(k) - x.col(k | (1 << d))).norm();
                        he = math::max(he, len);
                        if (d == 0) hTorMax = math::max(hTorMax, len);
                        else        hCrossMax = math::max(hCrossMax, len);
                    }
            hMax = math::max(hMax, he); hMin = math::min(hMin, he); hSum += he; ++nE;
        }
        hSum /= nE;
    }
    gsInfo << "Element size (longest edge per element): h_max = " << hMax << " mm, h_min = " << hMin
           << " mm, h_mean = " << hSum << " mm; longest toroidal edge " << hTorMax
           << " mm, longest cross-section edge " << hCrossMax << " mm\n";
    if (hMax > 8 * math::sqrt(par.D / par.rho))
        gsWarn << "h_max = " << hMax << " mm is much larger than the front width ~ 4 sqrt(D/rho) = "
               << 4 * math::sqrt(par.D / par.rho) << " mm: refine (-r).\n";

    // 2. Generalized-alpha
    const real_t rho_inf = TIMEopt.askReal("rho_inf", 0.5);
    const real_t alpha_m = 0.5 * (3.0 - rho_inf) / (1.0 + rho_inf);
    const real_t alpha_f = 1.0 / (1.0 + rho_inf);
    const real_t gamma   = 0.5 + alpha_m - alpha_f;
    const real_t D = par.D, rho = par.rho;
    gsInfo << "Model: D = " << D << " mm^2/day, rho = " << rho << " 1/day, R = " << par.R
           << " mm, sigma = " << par.sigma << " mm, amp = " << par.amp << ", N_th = " << par.Nth << "\n";

    real_t V_T, V_N, V_P, N_dash, N_T, A_p, int_N, int_H;
    std::ofstream csvFile, compFile;
    gsFileManager::mkdir(out);
    csvFile.open(out + "/biomarkers.csv");
    csvFile << std::setprecision(12)
            << "TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p,N_max,N_min,"
               "V_ref,V_P_relerr,V_T,V_N,growth,dintN_dt\n";
    compFile.open(out + "/computational.csv");
    compFile << std::setprecision(8)
             << "TimeStep,Time,NewtonIts,LinIts,ResidualNorm,StepSolve_s,StepTotal_s,Wall_s\n";
    const real_t k_steep = 200.0;
    gsConstantFunction<> alpha_f_func(alpha_f, 3), alpha_m_func(alpha_m, 3);

    // 3. Assembler
    gsExprAssembler<> A(1, 1);
    A.options().setSwitch("SameElement", Aopt.askSwitch("SameElement", true));
    A.setIntegrationElements(dbasis);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;
    geometryMap G = A.getMap(mp);
    space w = A.getSpace(dbasis);
    gsMatrix<> Nnew, dNnew, Nold, dNold, Nalpha, dNalpha, dNupdate, Q;
    gsMultiPatch<> mp_nold, mp_dnold;
    auto nold = A.getCoeff(mp_nold);
    auto dnold = A.getCoeff(mp_dnold);
    solution nnew_sol = A.getSolution(w, Nnew);
    solution dnnew_sol = A.getSolution(w, dNnew);
    auto af = A.getCoeff(alpha_f_func);
    auto am = A.getCoeff(alpha_m_func);
    auto N     = nold + af * (nnew_sol - nold);
    auto dN    = dnold + am * (dnnew_sol - dnold);
    auto gradN = igrad(nold, G) + af * (igrad(nnew_sol, G) - igrad(nold, G));
    auto residual = w * dN + D * igrad(w, G) * gradN.tr() + (-rho) * w * N.val() * (1.0 - N.val());

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(mp);
    w.setup(bc, dirichlet::l2Projection, 0);
    w.setupMapper(torusMapper(dbasis, st.polarMatch, false));
    A.initSystem();
    gsInfo << "Discretization: numRefine = " << numRefine << ", " << sb.numElements() << " elements, "
           << sb.size() << " basis functions, DoFs (unknowns) = " << A.numDofs() << " (seams glued C0"
           << (st.polarMatch ? ", polar axis merged" : "") << "), dt = " << dt << ", steps = " << maxSteps << "\n";
    if (const gsTensorBSplineBasis<3, T> *tb = dynamic_cast<const gsTensorBSplineBasis<3, T> *>(&sb))
    {
        const index_t nu = tb->size(0), nv = tb->size(1), nw = tb->size(2);
        const index_t expected = st.polarMatch ? (nu - 1) * ((nv - 1) * (nw - 1) + 1) : (nu - 1) * (nv - 1) * nw;
        GISMO_ENSURE(A.numDofs() == expected, "Unexpected number of DoFs " << A.numDofs() << " != " << expected
                     << ": seam/axis gluing failed.");
    }

    gsSparseMatrix<> K, K_const;
    gsInfo << "Linear solver: " << st.linSolver
           << (st.linSolver == "direct" ? "" : " (CG tol " + util::to_string(st.cgTol) + ")") << "\n";

    // 4. Initial condition: L2 projection of the seed at PHYSICAL points
    gsFunctionExpr<T> N0 = tumorSeed<T>(par);
    gsMatrix<T> N0coefs;
    real_t icRelErr = 0, tIC = 0;
    {
        gsStopwatch swIC;
        auto f0 = A.getCoeff(N0, G);
        A.initSystem();
        A.assemble(w * w.tr() * meas(G), w * f0 * meas(G));
        fkSolveLinear<T>(st.linSolver, A.matrix(), A.rhs(), Nold, st.cgMaxIt, st.cgTol);
        tIC = swIC.stop();
        solution n0_sol = A.getSolution(w, Nold);
        n0_sol.extract(N0coefs, 0);
        gsExprEvaluator<T> ev0(A);
        const real_t errIC = math::sqrt(ev0.integral((n0_sol - f0).sqNorm() * meas(G)));
        const real_t nrmIC = math::sqrt(ev0.integral(f0.sqNorm() * meas(G)));
        icRelErr = errIC / nrmIC;
        gsInfo << "IC L2-projection error (physical space): " << errIC << ", relative " << icRelErr << "\n";
        A.initSystem();
    }
    mp_nold.addPatch(sb.makeGeometry(N0coefs));
    mp_dnold.addPatch(sb.makeGeometry(gsMatrix<T>::Zero(sb.size(), 1)));
    dNold.setZero(A.numDofs(), 1);

    // 5. Constant part of the tangent
    A.assemble(meas(G) * (w * w.tr() * alpha_m + (alpha_f * gamma * dt) * (D * igrad(w, G) * igrad(w, G).tr())));
    K_const = A.giveMatrix();

    // 6. ParaView output of N: unstructured grid with merged seams and the
    //    collapsed core (see torusPlotGrid), one .vtu per frame
    const gsMatrix<T> parBox = mp.patch(0).support();
    const index_t npl = static_cast<index_t>(math::ceil(math::pow((T)st.plotPts, (T)1 / 3)));
    const torusPlotGrid tg(npl, npl, npl / 2 + 1);
    const gsMatrix<T> plotXYZ = st.plot ? mp.patch(0).eval(tg.par) : gsMatrix<T>();
    gsParaviewCollection collection(out + "/N");
    if (st.plot) gsFileManager::mkdir(out + "/N_pvd");
    index_t frame = 0;
    auto writeFrame = [&](T time)
    {
        std::vector<std::pair<std::string, gsMatrix<> > > fields;
        fields.push_back(std::make_pair(std::string("N"), mp_nold.patch(0).eval(tg.par)));
        const std::string rel = "N_pvd/N_" + util::to_string(frame++);
        tg.write(out + "/" + rel, plotXYZ, fields);
        collection.addPart(rel + ".vtu", time);
    };
    const gsMatrix<T> sampleGrid = gsPointGrid<T>(parBox, 40000);
    auto rangeN = [&](T &Nmax, T &Nmin)
    {
        gsMatrix<T> vals = mp_nold.patch(0).eval(sampleGrid);
        Nmax = vals.maxCoeff(); Nmin = vals.minCoeff();
    };
    T Nmax, Nmin;
    rangeN(Nmax, Nmin);
    gsInfo << "Initial condition: max N = " << Nmax << ", min N = " << Nmin << "\n";
    const T icNmax = Nmax, icNmin = Nmin;
    if (st.plot) writeFrame(0);

    real_t growth = 0, int_N_prev = 0;
    auto writeBio = [&](index_t n, real_t time)
    {
        auto H = 1.0 / (1.0 + (-(k_steep) * (nnew_sol.val() - par.Nth)).exp());
        int_N = computeInt(dbasis, mp, nnew_sol.val(), A.options());
        int_H = computeInt(dbasis, mp, H, A.options());
        V_T = int_H;
        V_N = computeInt(dbasis, mp, nnew_sol.val() * H, A.options());
        V_P = computeInt(dbasis, mp, 1.0, A.options());
        N_dash = V_N / V_T;
        N_T = N_dash * V_T / V_P;
        A_p = computeInt(dbasis, mp, rho * nnew_sol.val() * (1.0 - nnew_sol.val()) * H, A.options());
        growth = computeInt(dbasis, mp, rho * nnew_sol.val() * (1.0 - nnew_sol.val()), A.options());
        rangeN(Nmax, Nmin);
        gsInfo << "  int_N = " << int_N << " mm^3, V_T = " << V_T << " mm^3, max N = " << Nmax
               << ", min N = " << Nmin << "\n";
        csvFile << n << "," << time << "," << A.numDofs() << "," << int_N << "," << int_H << "," << V_P << ","
                << N_dash << "," << N_T << "," << A_p << "," << Nmax << "," << Nmin << "," << volRef << ","
                << math::abs(V_P - volRef) / volRef << "," << V_T << "," << V_N << "," << growth << ",";
        if (n > 0) csvFile << (int_N - int_N_prev) / dt;
        csvFile << "\n";
        csvFile.flush();
        int_N_prev = int_N;
    };
    auto saveN = [&](index_t n)
    {
        if (!saveSteps.count(n)) return;
        const std::string fn = out + "/N_t" + util::to_string(n * dt) + ".xml";
        gsWrite(mp_nold.patch(0), fn);
        gsInfo << "  Solution written: " << fn << "\n";
    };

    Nnew = Nold;
    writeBio(0, 0.0);
    saveN(0);
    const real_t V_P_relerr = math::abs(V_P - volRef) / volRef;
    const real_t tSetup = swTotal.stop();

    // 7. Time loop
    const index_t maxIt = 50;
    const real_t tol = TIMEopt.askReal("tol", 1e-4);
    index_t itsTotal = 0, linItsTotal = 0;
    real_t tSolveTotal = 0;
    gsStopwatch swLoop;
    for (index_t step = 0; step < maxSteps; ++step)
    {
        gsInfo << "Time Step " << step + 1 << " / " << maxSteps << "\n";
        gsStopwatch swStep;
        Nnew = Nold;
        dNnew = (gamma - 1.0) / gamma * dNold;
        index_t its = maxIt, linIts = 0;
        real_t resNorm = 0;
        for (index_t it = 0; it < maxIt; ++it)
        {
            Nalpha = Nold + alpha_f * (Nnew - Nold);
            dNalpha = dNold + alpha_m * (dNnew - dNold);
            A.initMatrix();
            A.clearRhs();
            A.assemble(residual * meas(G));
            Q = A.rhs();
            resNorm = Q.norm();
            if (it > 0 && resNorm < tol)
            {
                gsInfo << "  Converged in " << it << " iterations.\n";
                its = it;
                break;
            }
            A.initMatrix();
            A.assemble(meas(G) * (-alpha_f * gamma * dt) * (rho * w * (1.0 - 2.0 * N.val()) * w.tr()));
            K = A.giveMatrix();
            K += K_const;
            const index_t li = fkSolveLinear<T>(st.linSolver, K, -Q, dNupdate, st.cgMaxIt, st.cgTol);
            linIts += li; linItsTotal += li;
            dNnew += dNupdate;
            Nnew.noalias() += (gamma * dt) * dNupdate;
        }
        if (its == maxIt)
            gsWarn << "Newton did not converge in " << maxIt << " iterations (|Q| = " << resNorm << ").\n";
        const real_t tSolve = swStep.stop();
        nnew_sol.extract(mp_nold);
        dnnew_sol.extract(mp_dnold);
        Nold = Nnew;
        dNold = dNnew;
        if (st.plot) writeFrame((step + 1) * dt);
        writeBio(step + 1, (step + 1) * dt);
        saveN(step + 1);
        itsTotal += its;
        tSolveTotal += tSolve;
        compFile << step + 1 << "," << (step + 1) * dt << "," << its << "," << linIts << "," << resNorm << ","
                 << tSolve << "," << swStep.stop() << "," << swTotal.stop() << "\n";
        compFile.flush();
    }
    csvFile.close();
    compFile.close();
    if (st.plot) collection.save();

    const real_t tLoop = swLoop.stop(), tTotal = swTotal.stop();
    std::ofstream sumFile((out + "/computational_summary.csv").c_str());
    sumFile << std::setprecision(10) << "quantity,value\n"
            << "method," << st.method << "\nnumRefine," << numRefine << "\ndegree," << sb.degree(0) << "\n"
            << "seam,glued_C0\npolar_axis," << (st.polarMatch ? "merged" : "free") << "\n"
            << "elements," << sb.numElements() << "\n"
            << "h_max_mm," << hMax << "\nh_min_mm," << hMin << "\nh_mean_mm," << hSum << "\n"
            << "h_tor_max_mm," << hTorMax << "\nh_cross_max_mm," << hCrossMax << "\n"
            << "DoFs," << A.numDofs() << "\n"
            << "time_IC_projection_s," << tIC << "\nIC_L2_relerr," << icRelErr << "\n"
            << "IC_Nmax," << icNmax << "\nIC_Nmin," << icNmin << "\n"
            << "V_ref_mm3," << volRef << "\nV_P_relerr," << V_P_relerr << "\n"
            << "dt," << dt << "\nsteps," << maxSteps << "\n"
            << "newton_its_total," << itsTotal << "\n"
            << "newton_its_mean," << (maxSteps ? (real_t)itsTotal / maxSteps : 0) << "\n"
            << "lin_solver," << st.linSolver << "\ncg_tol," << st.cgTol << "\n"
            << "lin_its_total," << linItsTotal << "\n"
            << "lin_its_per_newton_mean," << (itsTotal ? (real_t)linItsTotal / itsTotal : 0) << "\n"
            << "time_setup_s," << tSetup << "\ntime_solve_total_s," << tSolveTotal << "\n"
            << "time_loop_s," << tLoop << "\ntime_per_step_mean_s," << (maxSteps ? tLoop / maxSteps : 0) << "\n"
            << "time_total_s," << tTotal << "\n";
    gsInfo << "Times: setup " << tSetup << " s, time loop " << tLoop << " s ("
           << (maxSteps ? tLoop / maxSteps : 0) << " s/step, "
           << (maxSteps ? (real_t)itsTotal / maxSteps : 0) << " Newton its/step), total " << tTotal << " s\n";
}

/// Default seed on the core line G(seedU,0,0), distance check, regime diagnostics
inline void prepareSeed(fkParams<real_t> &par, const gsGeometry<> &geo, const TriMesh &m,
                        real_t seedU, real_t T)
{
    if (math::isnan(par.xc) && math::isnan(par.yc) && math::isnan(par.zc))
    {
        gsMatrix<> u(3, 1), x;
        u << seedU, 0, 0;
        geo.eval_into(u, x);
        par.xc = x(0, 0); par.yc = x(1, 0); par.zc = x(2, 0);
    }
    GISMO_ENSURE(!math::isnan(par.xc) && !math::isnan(par.yc) && !math::isnan(par.zc),
                 "Give all three of --xc, --yc, --zc (or none).");
    Vec3 q;
    const real_t d = m.closest(Vec3(par.xc, par.yc, par.zc), q);
    gsInfo << "Seed centre (" << par.xc << ", " << par.yc << ", " << par.zc << ") mm, signed distance to the "
           << "surface " << d << " mm\n";
    if (d > 0) gsWarn << "The seed centre lies OUTSIDE the prostate surface.\n";
    else if (par.R + 2 * par.sigma > -d)
        gsWarn << "Tumour seed R + 2 sigma = " << par.R + 2 * par.sigma
               << " mm exceeds the distance to the surface: the seed is cut by the boundary.\n";
    gsInfo << "Regime over T = " << T << " days: rho*T = " << par.rho * T << ", D/R^2 = "
           << par.D / (par.R * par.R) << " 1/day vs rho = " << par.rho << " 1/day, front speed 2 sqrt(D rho) = "
           << 2 * math::sqrt(par.D * par.rho) << " mm/day, front width ~ " << 4 * math::sqrt(par.D / par.rho) << " mm\n";
    if (par.rho * T < 1) gsWarn << "rho*T < 1: almost no proliferation within the horizon.\n";
    if (par.D / (par.R * par.R) > par.rho) gsWarn << "D/R^2 > rho: diffusion-dominated, the tumour will be diluted.\n";
}

} // namespace fkp
} // namespace gismo

using namespace gismo;
using namespace gismo::fkp;

int main(int argc, char *argv[])
{
    std::string fn("pde/fisher_kolmogorov_cube_bvp.xml");
    for (int i = 1; i + 1 < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-f" || arg == "--file") fn = argv[i + 1];
    }
    gsFileData<> fd(fn);
    gsInfo << "Loaded configuration file: " << fd.lastPath() << "\n";
    gsOptionList TIMEopt, Aopt;
    fd.getId(4, TIMEopt); // rho_inf, tol
    fd.getId(5, Aopt);    // SameElement

    fkParams<real_t>   par;
    fkSolverSettings   st;
    lsOptions          lso;
    correctionOptions  co;
    std::string meshFile = "optional/gsGmsh/filedata/prostata_lorenzo.obj";
    std::string fitMode  = "ls";
    index_t degree = 2, nCirc = 32, nRad = 8, clDegree = 3, numRefine = 3, geoPlotPts = 20000;
    real_t  lumenMax = 0, seedU = 0, minSJ = 0.01;
    bool geometryOnly = false, noGeoPlot = false, noPolarMatch = false, C0fit = false;
    st.out = "00_output_FK_prostate_bf";
    st.method = "bodyfitted_prostate_tensor";

    gsCmdLine cmd("Fisher-Kolmogorov on a body-fitted (tensor B-spline) prostate mapped from a solid torus");
    cmd.addString("f", "file", "Input XML setup file (rho_inf, tol)", fn);
    cmd.addString("m", "mesh", "Prostate surface with urethra (OBJ or STL, genus 1, mm)", meshFile);
    cmd.addString("o", "output", "Output directory", st.out);
    cmd.addInt   ("d", "degree", "Spline degree", degree);
    cmd.addInt   ("c", "circ", "Coarse elements in u (toroidal) and v (poloidal)", nCirc);
    cmd.addInt   ("", "radial", "Coarse elements in w (radial)", nRad);
    cmd.addInt   ("r", "refine", "Uniform knot-insertion steps (3: 32x32x8 -> 256x256x64)", numRefine);
    cmd.addInt   ("", "clDegree", "Polynomial degree of the urethra centreline fit", clDegree);
    cmd.addReal  ("", "lumenMax", "Max. lumen diameter used to detect urethra triangles [mm] (0: auto)", lumenMax);
    cmd.addString("", "fit", "Volume construction: ls (least squares, default) | interp (legacy)", fitMode);
    cmd.addInt   ("", "fitSamples", "ls: samples per element in each direction", lso.samples);
    cmd.addInt   ("", "centreModes", "ls: Fourier modes of the section centre c(phi) (0: raw)", lso.centreModes);
    cmd.addReal  ("", "curvWeight", "ls: v-parametrization, <0 polar angle, >=0 arc length + weight x turning angle", lso.curvWeight);
    cmd.addSwitch("C0fit", "ls/correction: C0 instead of C1 across the u- and v-seams", C0fit);
    cmd.addInt   ("", "fitIts", "Closest-point surface correction iterations (0: none)", co.its);
    cmd.addReal  ("", "fitLambda", "Correction: relative weight of the parametric Dirichlet energy", co.lambda);
    cmd.addReal  ("", "minSJ", "Correction: smallest accepted scaled Jacobian (steps are halved below it)", minSJ);
    cmd.addSwitch("geometryOnly", "Build, refine, check and export the geometry; no time stepping", geometryOnly);
    cmd.addSwitch("noGeoPlot", "Do not write the geometry/mesh ParaView files", noGeoPlot);
    cmd.addInt   ("", "geoPlotPts", "Sample points of the scaled-Jacobian ParaView field", geoPlotPts);
    cmd.addSwitch("noPolarMatch", "Do not merge the DoFs on the collapsed face w = 0", noPolarMatch);
    cmd.addReal("t", "dt", "Time step size [day]", st.dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", st.maxSteps);
    cmd.addSwitch("plot", "Write ParaView files of N for every time step", st.plot);
    cmd.addInt("", "plotPts", "ParaView sample points per frame of N", st.plotPts);
    cmd.addString("", "solver", "Linear solver: cgdiag (default) | cgic | direct", st.linSolver);
    cmd.addInt   ("", "cgMaxIt", "Max CG iterations per linear solve", st.cgMaxIt);
    cmd.addReal  ("", "cgTol", "CG tolerance, relative to ||rhs||", st.cgTol);
    cmd.addString("", "saveAt", "Write the solution spline at these times [day], e.g. 365,730", st.saveAt);
    cmd.addReal("D", "diffusivity", "Tumor cell diffusivity D [mm^2/day]", par.D);
    cmd.addReal("p", "rho", "Net proliferation rate rho [1/day]", par.rho);
    cmd.addReal("R", "radius", "Radius R of the initial tumor [mm]", par.R);
    cmd.addReal("s", "sigma", "Interface width sigma of the initial tumor [mm]", par.sigma);
    cmd.addReal("a", "amp", "Peak density of the initial tumor (0..1]", par.amp);
    cmd.addReal("", "Nth", "Cell density threshold N_th defining Omega_T", par.Nth);
    cmd.addReal("", "xc", "x-coordinate of the tumor centre [mm] (default: on the core line)", par.xc);
    cmd.addReal("", "yc", "y-coordinate of the tumor centre [mm] (default: on the core line)", par.yc);
    cmd.addReal("", "zc", "z-coordinate of the tumor centre [mm] (default: on the core line)", par.zc);
    cmd.addReal("", "seedU", "Toroidal parameter u in [0,1) of the default seed centre G(u,0,0)", seedU);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree >= 1 && nCirc >= 3 && nRad >= 1 && numRefine >= 0, "Invalid discretization.");
    GISMO_ENSURE(fitMode == "ls" || fitMode == "interp", "--fit must be ls or interp.");
    GISMO_ENSURE(par.amp > 0 && par.amp <= 1, "The seed amplitude must lie in (0,1].");
    lso.C1 = !C0fit;
    st.polarMatch = !noPolarMatch;
    const std::string &out = st.out;
    gsFileManager::mkdir(out);
    // all geometry / mesh / fit-quality files go to <out>/geometry/; the solver
    // output (biomarkers.csv, computational*.csv, N.pvd) stays in <out>/
    const std::string geoDir = out + "/geometry";
    gsFileManager::mkdir(geoDir);

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------
    gsStopwatch swGeo;
    ProstateSurface S;
    S.build(meshFile, clDegree, lumenMax);

    index_t nonStar = 0;
    gsMultiPatch<> mp;
    mp.addPatch(fitMode == "ls" ? buildTensorLS(S, degree, nCirc, nRad, lso, nonStar)
                                : buildTensorInterp(S, degree, nCirc, nRad, nonStar));
    setTorusTopology(mp);
    if (nonStar > 0)
        gsWarn << nonStar << " polar rays cross the prostate boundary more than once (section not "
               << "star-shaped w.r.t. its centre): the ray fit is only approximate there.\n";
    gsMultiBasis<> mbCoarse(mp, true);
    const deviationStats dev0 = surfaceDeviation(mp.patch(0), S.m, "Deviation (" + fitMode + ")");

    // Closest-point correction of the surface (seam constraints as in the fit)
    deviationStats dev = dev0;
    if (co.its > 0)
    {
        const gsTensorBSplineBasis<3, real_t> &tb =
            static_cast<const gsTensorBSplineBasis<3, real_t> &>(mp.patch(0).basis());
        // the interp volume is only C0 at the seams: correct it with C0 constraints
        const gsSparseMatrix<> Z = tensorZ(tb, fitMode == "ls" && lso.C1);
        real_t qCur = minOf(elementQuality(mp.patch(0), mbCoarse.basis(0)));
        auto valid = [&](const gsGeometry<> &g)
        {
            const real_t qn = minOf(elementQuality(g, mbCoarse.basis(0)));
            if (qn > 0 && qn >= math::min(minSJ, qCur)) { qCur = qn; return true; }
            return false;
        };
        gsInfo << "Closest-point surface correction (" << co.its << " iterations):\n";
        if (!correctSurface(mp.patch(0), Z, S.m, co, valid))
            gsWarn << "A correction step was rejected (it would fold the volume); try a larger --fitLambda.\n";
        dev = surfaceDeviation(mp.patch(0), S.m, "Deviation (after correction)");
    }

    const jacobianStats jc = jacobianReport(mp, mbCoarse, "Coarse mesh",
                                            noGeoPlot ? "" : geoDir + "/geometry_coarse_ScaledJacobian", geoPlotPts);
    gsInfo << "Volume: IGA " << jc.vol << " vs. triangulation " << S.m.vol << " mm^3 ("
           << 100 * (jc.vol - S.m.vol) / S.m.vol << " %)\n";
    {
        const std::vector<real_t> q = elementQuality(mp.patch(0), mbCoarse.basis(0));
        gsVector<> wq, wx;
        const real_t sw = worstScaledJacobian(mp.patch(0), wq, wx);
        gsInfo << "Element quality (min scaled Jacobian incl. faces/corners): min " << minOf(q)
               << "; elements " << qualityHistogram(q) << " of " << q.size() << "\n"
               << "Worst scaled Jacobian on a 64^3 parameter grid: " << sw << " at (u,v,w) = ("
               << wq.transpose() << "), x = (" << wx.transpose() << ") mm\n";
        GISMO_ENSURE(minOf(q) > 0, "The volume is folded at element faces/corners (min scaled Jacobian "
                     << minOf(q) << "): lower --fitIts or raise --fitLambda / --minSJ.");
    }
    GISMO_ENSURE(jc.sjMin > 0, "The fitted volume is folded (min scaled Jacobian " << jc.sjMin
                 << " <= 0): try --clDegree 2, --fitLambda, another --circ, or check the surface.");

    // Global refinement by knot insertion; the map is unchanged, the solver
    // uses the coarse map with the refined basis.
    gsMultiBasis<> dbasis(mp, true);
    for (index_t r = 0; r != numRefine; ++r) dbasis.uniformRefine();
    const jacobianStats jr = jacobianReport(mp, dbasis, "Refined mesh (quadrature points of the solution basis)");
    GISMO_ENSURE(jr.sjMin > 0, "Folded at the refined quadrature points.");

    gsWrite(mp, geoDir + "/prostate_torus_coarse.xml");
    gsInfo << "Coarse geometry (seams as interfaces) written to " << geoDir << "/prostate_torus_coarse.xml\n";

    if (!noGeoPlot)
    {
        surfaceDeviation(mp.patch(0), S.m, "Deviation (written)", geoDir + "/deviation");
        S.m.writeVtp(geoDir + "/input_surface");
        // whole volume (solid) and its closed boundary as OBJ, for comparison with the input
        writeTorusVolume(mp.patch(0), geoDir + "/geometry_volume"); // merged seams, scaled Jacobian field
        const real_t vObj = writeSurfaceObj(mp.patch(0), geoDir + "/geometry_surface_fit.obj");
        gsInfo << "Whole volume: " << geoDir << "/geometry_volume.vtu; its boundary (w = 1) as OBJ: "
               << geoDir << "/geometry_surface_fit.obj (enclosed volume " << vObj << " mm^3 vs input "
               << S.m.vol << " mm^3)\n";
        gsWriteParaviewPoints(toMatrix(S.cl.P), geoDir + "/centreline");
        gsWriteParaviewPoints(toMatrix(S.lumenPts), geoDir + "/lumen_points");
        gsWriteParaview(*mp.patch(0).boundary(boundary::back), geoDir + "/geometry_surface", 40000);
        gsWriteParaview(*mp.patch(0).boundary(boundary::west), geoDir + "/geometry_crosssection", 10000);
        writeSingleCompMesh(mbCoarse.basis(0), mp.patch(0), geoDir + "/mesh_coarse_volume", 4);
        gsBasis<>::uPtr bs = dbasis.basis(0).boundaryBasis(boundary::back);
        gsBasis<>::uPtr bcs = dbasis.basis(0).boundaryBasis(boundary::west);
        writeSingleCompMesh(*bs, *mp.patch(0).boundary(boundary::back), geoDir + "/mesh_refined_surface", 2);
        writeSingleCompMesh(*bcs, *mp.patch(0).boundary(boundary::west), geoDir + "/mesh_refined_crosssection", 2);
        if (dbasis.basis(0).numElements() <= 300000)
            writeSingleCompMesh(dbasis.basis(0), mp.patch(0), geoDir + "/mesh_refined_volume", 2);
        else
            gsInfo << "Refined volume wireframe skipped (" << dbasis.basis(0).numElements() << " elements).\n";
        gsInfo << "ParaView files written to " << geoDir << "/ (deviation_spline: signed distance on the "
               << "spline surface, deviation_input: distance of the input vertices; input_surface, "
               << "centreline, geometry_*, mesh_*, geometry_coarse_ScaledJacobian)\n";
    }
    gsInfo << "Geometry built, checked and exported in " << swGeo.stop() << " s\n";

    {
        std::ofstream g((geoDir + "/geometry_summary.csv").c_str());
        g << std::setprecision(10) << "quantity,value\n"
          << "surface," << meshFile << "\nfit," << fitMode << "\nseams_fit," << (fitMode == "ls" && lso.C1 ? "C1" : "C0") << "\n"
          << "fitIts," << co.its << "\nfitLambda," << co.lambda << "\n"
          << "coarse_elements," << mbCoarse.basis(0).numElements() << "\n"
          << "refined_elements," << dbasis.basis(0).numElements() << "\n"
          << "numRefine," << numRefine << "\ndegree," << degree << "\n"
          << "scaledJac_min_coarse," << jc.sjMin << "\nscaledJac_max_coarse," << jc.sjMax << "\n"
          << "scaledJac_min_refined," << jr.sjMin << "\ndetJ_min_refined," << jr.detMin << "\n"
          << "V_IGA_mm3," << jc.vol << "\nV_ref_mm3," << S.m.vol << "\n"
          << "V_relerr," << math::abs(jc.vol - S.m.vol) / S.m.vol << "\n"
          << "hausdorff_fit_mm," << dev0.hausdorff() << "\n"
          << "hausdorff_mm," << dev.hausdorff() << "\n"
          << "dist_surf_to_mesh_mean_mm," << dev.fwdMean << "\ndist_surf_to_mesh_rms_mm," << dev.fwdRms << "\n"
          << "dist_surf_to_mesh_max_mm," << dev.fwdMax << "\n"
          << "dist_mesh_to_surf_mean_mm," << dev.bwdMean << "\ndist_mesh_to_surf_rms_mm," << dev.bwdRms << "\n"
          << "dist_mesh_to_surf_max_mm," << dev.bwdMax << "\n"
          << "non_star_rays," << nonStar << "\n";
    }
    if (geometryOnly) return EXIT_SUCCESS;

    prepareSeed(par, mp.patch(0), S.m, seedU, st.dt * st.maxSteps);
    solveFK<real_t>(mp, dbasis, TIMEopt, Aopt, par, st, numRefine, S.m.vol);
    return EXIT_SUCCESS;
}
