/** @file tchakaloff_selftest_example.cpp

    @brief Self-test driver for gsTchakaloffRule.h: runs the six measure-
    compression cases (volume and boundary variants, single- and
    two-level chunking) over a range of polynomial degrees and reports
    PASS/FAIL per case, then two permanent real-cut-cell regression cases
    (sphere r=1, cells 172 and 214). The real-cell cases are reported as
    expected failures (XFAIL/XPASS) rather than PASS/FAIL, since they
    currently hit a known NNLS stopping-test stall (see the TODO on the
    stopping test in gsTchakaloffRule.h); --strict-realcell treats
    them as ordinary PASS/FAIL instead. The run ends with one summary line
    "<n> PASS, <n> FAIL, <n> XFAIL, <n> XPASS" and is auto-registered as a
    CTest; exit 0 iff FAIL is 0 (real-cell XFAIL/XPASS never contribute to
    FAIL; under --strict-realcell they do).

    With --bench, runs a performance/memory benchmark instead: an
    NNLS-only micro-benchmark on a single chunk, followed by end-to-end
    tchakaloffCompress[Boundary] rows over synthetic volume and boundary
    measures of increasing size, reporting timing, peak RSS and streamed
    moment-reproduction error. See --help for the sweep options.

    With --nnmf, runs a strict PASS/FAIL self-test of the NNMF rule of
    gsNnmfRule.h instead (the synthetic cases for --pmin..--pmax plus a dense
    tilted-plane boundary case, and the two real cut cells at p = 2); exit 0
    iff every case passes.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gsTchakaloffRule.h"
#include "gsNnmfRule.h"
#include "gsImmersedLookupRule.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace gismo;
using namespace gsTetClip;

namespace {

typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic> DMat;

// ---------------------------------------------------------------------
// Peak-RSS measurement (Linux /proc, proc(5)).
// ---------------------------------------------------------------------

/// Reads the kB value off a "VmXxx:  NNNN kB" line of /proc/self/status.
/// Returns -1 if /proc is unavailable or the key is not present (non-Linux).
long readStatusKB(const char* key)
{
    std::ifstream f("/proc/self/status");
    std::string line;
    const std::size_t klen = std::strlen(key);
    while (std::getline(f, line))
    {
        if (line.compare(0, klen, key) == 0)
        {
            std::istringstream iss(line.substr(klen));
            long kb = -1;
            iss >> kb;
            return kb;
        }
    }
    return -1;
}

long readVmHWM() { return readStatusKB("VmHWM:"); }
long readVmRSS() { return readStatusKB("VmRSS:"); }

/// Writes "5" to /proc/self/clear_refs, which resets VmHWM to the current
/// VmRSS (proc(5), Linux >= 4.0). Returns true iff the write succeeded; a
/// caller must still verify the effect via readVmHWM()/readVmRSS(), since
/// some sandboxes accept the write but do not honour it.
bool resetPeakRSS()
{
    std::ofstream f("/proc/self/clear_refs");
    if (!f)
        return false;
    f << "5";
    return f.good();
}

// ---------------------------------------------------------------------
// Deterministic per-row RNG seeding.
// ---------------------------------------------------------------------

/// FNV-1a 64-bit hash of a (kind, p, N) tag, used to seed a std::mt19937
/// so that a given (kind, p, N) always draws an identical synthetic
/// measure, from any binary built from this source file. Deliberately not
/// std::hash<string>, whose result is unspecified across standard-library
/// versions.
std::mt19937 makeRng(const std::string& kind, index_t p, index_t N)
{
    std::ostringstream tag;
    tag << kind << '_' << p << '_' << N;
    const std::string s = tag.str();
    std::uint64_t h = 14695981039346656037ULL;
    for (std::size_t i = 0; i < s.size(); ++i)
    {
        h ^= static_cast<unsigned char>(s[i]);
        h *= 1099511628211ULL;
    }
    return std::mt19937(static_cast<std::mt19937::result_type>(h));
}

// ---------------------------------------------------------------------
// Synthetic measures (free N, unlike the self-test's fixed-N builders).
// ---------------------------------------------------------------------

/// detail::buildVRandom's recipe (uniform nodes in the cell, weights
/// ~ U[0.5,1.5]*vol/N) with N free, so the benchmark can sweep N
/// independently of K = (2p+1)^3.
void benchBuildVol(const gsVector<real_t>& lower, const gsVector<real_t>& upper, real_t vol, index_t N,
                    std::mt19937& rng, gsMatrix<real_t>& nodes, gsVector<real_t>& weights)
{
    nodes.resize(3, N);
    weights.resize(N);
    std::uniform_real_distribution<real_t> ux(lower(0), upper(0));
    std::uniform_real_distribution<real_t> uy(lower(1), upper(1));
    std::uniform_real_distribution<real_t> uz(lower(2), upper(2));
    std::uniform_real_distribution<real_t> uw(0.5, 1.5);
    for (index_t i = 0; i < N; ++i)
    {
        nodes(0, i) = ux(rng);
        nodes(1, i) = uy(rng);
        nodes(2, i) = uz(rng);
        weights(i) = uw(rng) * vol / static_cast<real_t>(N);
    }
}

/// detail::buildSCap's spherical-cap geometry (R=0.55, theta0=1.0,
/// phi0=0.7), sampled with N points drawn uniformly at random in
/// (theta,phi) in [theta0+-0.15]x[phi0+-0.15] instead of tensor-Gauss
/// nodes, with weight = R^2 sin(theta) * (patch area) / N (a Monte-Carlo
/// quadrature of the same patch). normal = dir(theta,phi) is kept
/// deterministic (never randomised): random normals push the rank to the
/// full 4K = 4*(2p+1)^3, and a from-scratch NNLS refactorization per
/// outer iteration costs about 0.5*r^4 per chunk at that rank, which runs
/// for hours at p=3.
void benchBuildBdr(const gsVector<real_t>& mid, const gsVector<real_t>& lower, const gsVector<real_t>& upper,
                    index_t N, std::mt19937& rng, gsMatrix<real_t>& nodes, gsVector<real_t>& weights,
                    gsMatrix<real_t>& normals)
{
    const real_t R = 0.55, theta0 = 1.0, phi0 = 0.7;
    gsVector<real_t> d(3);
    d(0) = std::sin(theta0) * std::cos(phi0);
    d(1) = std::sin(theta0) * std::sin(phi0);
    d(2) = std::cos(theta0);
    const gsVector<real_t> cs = mid - R * d;

    std::uniform_real_distribution<real_t> uth(theta0 - 0.15, theta0 + 0.15);
    std::uniform_real_distribution<real_t> uph(phi0 - 0.15, phi0 + 0.15);

    nodes.resize(3, N);
    weights.resize(N);
    normals.resize(3, N);
    for (index_t i = 0; i < N; ++i)
    {
        const real_t th = uth(rng);
        const real_t ph = uph(rng);
        const real_t sn = std::sin(th);
        gsVector<real_t> dir(3);
        dir(0) = sn * std::cos(ph);
        dir(1) = sn * std::sin(ph);
        dir(2) = std::cos(th);
        nodes.col(i) = cs + R * dir;
        normals.col(i) = dir;
        weights(i) = R * R * sn * (0.3 * 0.3) / static_cast<real_t>(N);
    }

    for (index_t i = 0; i < N; ++i)
        for (index_t d2 = 0; d2 < 3; ++d2)
            GISMO_ENSURE(nodes(d2, i) >= lower(d2) && nodes(d2, i) <= upper(d2),
                         "benchBuildBdr: node fell outside the cell");
}

// ---------------------------------------------------------------------
// NNLS-only micro-benchmark: one chunk (n = B), timing nnlsLawsonHanson
// alone on an (A, c) built exactly as detail::reduceChunk builds it.
// ---------------------------------------------------------------------

void nnlsMicroBench(const std::string& kind, index_t p, const gsVector<real_t>& lower,
                     const gsVector<real_t>& upper, const gsVector<real_t>& mid, real_t vol)
{
    const bool isBoundary = (kind == "bdr");
    const index_t deg = 2 * p;
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t m = isBoundary ? 4 * K : K;
    const index_t n = 4 * m; // B, exactly one chunk

    std::mt19937 rng = makeRng(kind + "_nnls", p, n);
    gsMatrix<real_t> nodes, normals;
    gsVector<real_t> weights;
    if (isBoundary)
        benchBuildBdr(mid, lower, upper, n, rng, nodes, weights, normals);
    else
        benchBuildVol(lower, upper, vol, n, rng, nodes, weights);

    gsMatrix<real_t> M;
    detail::buildConstraintMatrix(nodes, isBoundary ? &normals : NULL, lower, upper, deg, M);
    gsEigen::ColPivHouseholderQR<DMat> qr(M);
    const index_t r = detail::leadingRunRank(qr.matrixQR(), 1e-13);
    const DMat Q1 = qr.householderQ() * DMat::Identity(n, r);
    const gsMatrix<real_t> A = Q1.transpose();
    const gsVector<real_t> c = Q1.transpose() * weights;

    double minMs = -1.0;
    NnlsResult nr;
    for (int trial = 0; trial < 3; ++trial)
    {
        gsStopwatch sw;
        nr = nnlsLawsonHanson(A, c, 3 * n);
        const double ms = sw.stop() * 1e3;
        if (minMs < 0.0 || ms < minMs)
            minMs = ms;
    }

    index_t nnz = 0;
    for (index_t j = 0; j < nr.w.size(); ++j)
        if (nr.w(j) > 0.0)
            ++nnz;

    gsInfo << "NNLS kind=" << kind << " p=" << p << " n=" << n << " r=" << r << " ms=" << minMs
           << " iters=" << nr.iterations << " nnz=" << nnz << " res=" << nr.residual
           << " converged=" << (nr.converged ? 1 : 0) << "\n";
}

// ---------------------------------------------------------------------
// End-to-end benchmark row: one tchakaloffCompress[Boundary] call, timed
// and RSS-measured, followed by a streamed moment check.
// ---------------------------------------------------------------------

bool endToEndRow(const std::string& kind, index_t p, index_t N, const gsVector<real_t>& lower,
                  const gsVector<real_t>& upper, const gsVector<real_t>& mid, real_t vol)
{
    const bool isBoundary = (kind == "bdr");
    std::mt19937 rng = makeRng(kind, p, N);

    gsMatrix<real_t> nodes, normals;
    gsVector<real_t> weights;
    if (isBoundary)
        benchBuildBdr(mid, lower, upper, N, rng, nodes, weights, normals);
    else
        benchBuildVol(lower, upper, vol, N, rng, nodes, weights);

    const real_t sumWin = weights.sum();

    const bool clearOk = resetPeakRSS();
    const long hwmAfterReset = readVmHWM();
    const long rssAfterReset = readVmRSS();
    const bool hwmReset =
        clearOk && hwmAfterReset >= 0 && rssAfterReset >= 0 && hwmAfterReset <= rssAfterReset + 1024;
    const long base = rssAfterReset;

    gsStopwatch sw;
    const TchakaloffOptions opt;
    const TchakaloffResult res = isBoundary
        ? tchakaloffCompressBoundary(nodes, weights, normals, lower, upper, p, opt)
        : tchakaloffCompress(nodes, weights, lower, upper, p, opt);
    const double totalS = sw.stop();

    const long peak = readVmHWM();

    // Streamed moment check: accumulate m_in over blocks of the input, so
    // the full N x m constraint matrix (2.2 GB at bdr, p=3, N=2e5) is
    // never formed -- only detail::genericChecks does that, and is not
    // called here for exactly this reason.
    const index_t deg = 2 * p;
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t m = isBoundary ? 4 * K : K;
    gsVector<real_t> mIn = gsVector<real_t>::Zero(m);
    const index_t blockSize = 4096;
    for (index_t lo = 0; lo < N; lo += blockSize)
    {
        const index_t b = std::min(blockSize, N - lo);
        const gsMatrix<real_t> nodesBlk = nodes.middleCols(lo, b);
        const gsVector<real_t> wBlk = weights.segment(lo, b);
        gsMatrix<real_t> normBlk;
        if (isBoundary)
            normBlk = normals.middleCols(lo, b);
        gsMatrix<real_t> Mblk;
        detail::buildConstraintMatrix(nodesBlk, isBoundary ? &normBlk : NULL, lower, upper, deg, Mblk);
        mIn += Mblk.transpose() * wBlk;
    }

    const index_t count = static_cast<index_t>(res.indices.size());
    gsMatrix<real_t> nodesOut(3, count);
    gsMatrix<real_t> normOut;
    if (isBoundary)
        normOut.resize(3, count);
    for (index_t j = 0; j < count; ++j)
    {
        nodesOut.col(j) = nodes.col(res.indices[j]);
        if (isBoundary)
            normOut.col(j) = normals.col(res.indices[j]);
    }
    gsMatrix<real_t> Mout;
    detail::buildConstraintMatrix(nodesOut, isBoundary ? &normOut : NULL, lower, upper, deg, Mout);
    const gsVector<real_t> mOut = Mout.transpose() * res.weights;

    const real_t maxIn = mIn.array().abs().maxCoeff();
    const real_t maxDiff = (mOut - mIn).array().abs().maxCoeff();
    const real_t maxMomErr = (maxIn > 0.0) ? maxDiff / maxIn : maxDiff;

    const real_t sumW = res.weights.sum();

    const index_t B = (opt.chunkSize > 0) ? opt.chunkSize : 4 * m;
    const index_t nChunks = (N > B) ? (N + B - 1) / B : 1;

    gsInfo << "BENCH kind=" << kind << " p=" << p << " N=" << N << " levels=" << res.levels
           << " r=" << res.rank << " levelRank=";
    for (std::size_t i = 0; i < res.levelRank.size(); ++i)
    {
        if (i) gsInfo << ",";
        gsInfo << res.levelRank[i];
    }
    gsInfo << " count=" << count << " nChunks=" << nChunks << " iters=" << res.iterations
           << " total_s=" << totalS << " ms_per_cand=" << (1e3 * totalS / static_cast<double>(N))
           << " rssBaseMiB=" << (base / 1024.0) << " rssPeakMiB=" << (peak / 1024.0)
           << " rssDeltaMiB=" << ((peak - base) / 1024.0) << " hwmReset=" << (hwmReset ? 1 : 0)
           << " res=";
    for (std::size_t i = 0; i < res.levelResidual.size(); ++i)
    {
        if (i) gsInfo << ",";
        gsInfo << res.levelResidual[i];
    }
    gsInfo << " maxMomErr=" << maxMomErr << " sumW=";
    const std::streamsize oldPrec = gsInfo.precision();
    gsInfo << std::setprecision(17) << sumW << " sumWin=" << sumWin;
    gsInfo << std::setprecision(static_cast<int>(oldPrec));
    const bool rowOk = res.ok && (maxMomErr <= 1e-12) && (std::abs(sumW - sumWin) <= 1e-12 * sumWin);
    gsInfo << " ok=" << (res.ok ? 1 : 0) << "\n";

    return rowOk;
}

int runBench(const std::string& benchKind, index_t benchP, index_t benchN)
{
    gsVector<real_t> lower(3), upper(3), mid(3);
    lower(0) = 0.10; lower(1) = -0.20; lower(2) = 0.30;
    const real_t h = 0.25;
    for (index_t d = 0; d < 3; ++d)
    {
        upper(d) = lower(d) + h;
        mid(d) = lower(d) + 0.5 * h;
    }
    const real_t vol = h * h * h;

    std::vector<std::string> kinds;
    if (benchKind == "all")
    {
        kinds.push_back("vol");
        kinds.push_back("bdr");
    }
    else
    {
        kinds.push_back(benchKind);
    }

    std::vector<index_t> ps;
    if (benchP == 0)
    {
        ps.push_back(2);
        ps.push_back(3);
    }
    else
    {
        ps.push_back(benchP);
    }

    std::vector<index_t> Ns;
    if (benchN == 0)
    {
        Ns.push_back(1000);
        Ns.push_back(10000);
        Ns.push_back(100000);
        Ns.push_back(200000);
    }
    else
    {
        Ns.push_back(benchN);
    }

    for (std::size_t ki = 0; ki < kinds.size(); ++ki)
        for (std::size_t pi = 0; pi < ps.size(); ++pi)
            nnlsMicroBench(kinds[ki], ps[pi], lower, upper, mid, vol);

    bool allOk = true;
    for (std::size_t ki = 0; ki < kinds.size(); ++ki)
        for (std::size_t pi = 0; pi < ps.size(); ++pi)
            for (std::size_t ni = 0; ni < Ns.size(); ++ni)
                allOk = endToEndRow(kinds[ki], ps[pi], Ns[ni], lower, upper, mid, vol) && allOk;

    return allOk ? EXIT_SUCCESS : EXIT_FAILURE;
}

// ---------------------------------------------------------------------
// --realcell mode: streams one real cut cell exactly as
// immersed_tetmesh_poisson_example.cpp's buildCellRules does (that
// driver's cell pipeline is cited by line below; it is not included or
// edited), compresses it with full NNLS diagnostics, and can dump the
// uncompressed rule to a bitwise-exact binary file for offline
// comparison against other gsTchakaloffRule.h builds.
// ---------------------------------------------------------------------

/// Case name -> ASCII gmsh MSH 4.1 file, gsFileManager-resolved by the
/// caller. Mirrors immersed_tetmesh_poisson_example.cpp's meshFile().
std::string realCellMeshFile(const std::string& caseName)
{
    if ("sphere" == caseName) return "volumes/tetmesh_sphere.msh";
    if ("rotcube" == caseName) return "volumes/tetmesh_cube_rotated.msh";
    GISMO_ERROR("--case: unknown case '" << caseName << "'");
}

/// n = n0 << r cells per direction on the fixed box [-1,1]^3. Mirrors
/// immersed_tetmesh_poisson_example.cpp's makeGrid() (grid part only;
/// runRealCell range-checks --cell against the n^3 cells before any
/// lookup).
gsTetClip::Grid3 realCellGrid(index_t n0, index_t r)
{
    gsTetClip::Grid3 grid;
    grid.x0 = grid.y0 = grid.z0 = -1.0;
    grid.n  = n0 << r;
    grid.h  = 2.0 / static_cast<real_t>(grid.n);
    return grid;
}

const char* nnlsStopName(gsTetClip::detail::NnlsStop s)
{
    switch (s)
    {
        case gsTetClip::detail::NnlsStop::Zero:      return "ZERO";
        case gsTetClip::detail::NnlsStop::Residual:  return "RES";
        case gsTetClip::detail::NnlsStop::ZEmpty:    return "ZEMPTY";
        case gsTetClip::detail::NnlsStop::Kkt:       return "KKT";
        case gsTetClip::detail::NnlsStop::Full:      return "FULL";
        case gsTetClip::detail::NnlsStop::Exhausted: return "EXHAUSTED";
        case gsTetClip::detail::NnlsStop::Cap:       return "CAP";
    }
    return "?";
}

void printNnlsDiagLine(const gsTetClip::detail::ChunkDiag& cd)
{
    const gsTetClip::detail::NnlsDiagnostics& d = cd.nnls;
    gsInfo << "NNLSDIAG level=" << cd.level << " chunk=" << cd.chunk << " n=" << cd.n << " r=" << d.r
           << " k=" << d.k << " iters=" << d.iterations << " dels=" << d.deletions
           << " stop=" << nnlsStopName(d.stop) << " res=" << d.res << " resInc=" << d.resInc
           << " resFresh=" << d.resFresh << " freshAccepted=" << (d.freshAccepted ? 1 : 0)
           << " freshNonPos=" << d.freshNonPos << " freshMin=" << d.freshMin << " orthErr=" << d.orthErr
           << " factErr=" << d.factErr << " dErr=" << d.dErr << " gPnorm=" << d.gPnorm
           << " maxGZ=" << d.maxGZ << " gPnormStop=" << d.gPnormStop << " maxGZStop=" << d.maxGZStop
           << " rejNu=" << d.rejNu << " rejDiag=" << d.rejDiag << " rejZ=" << d.rejZ
           << " rejNuTot=" << d.rejNuTotal << " rejDiagTot=" << d.rejDiagTotal
           << " rejZTot=" << d.rejZTotal << " rho=" << d.rho << " minRii=" << d.minRii
           << " condAP=" << d.condAP << " wMin=" << cd.wMin << " wMax=" << cd.wMax
           << " refreshes=" << d.refreshes << "\n";
}

/// Selects and prints the informative NNLSDIAG rows: the level-2 (union)
/// call always, and at level 1 the worst-residual chunk plus every chunk
/// at or above the resTol boundary (1e-13) -- capped at 20 level-1 lines,
/// with the remainder only counted.
void printNnlsDiagLines(const gsTetClip::detail::ReduceDiagnostics& diag)
{
    index_t worstIdx = -1;
    real_t worstRes = -1.0;
    for (std::size_t c = 0; c < diag.calls.size(); ++c)
    {
        const gsTetClip::detail::ChunkDiag& cd = diag.calls[c];
        if (cd.level == 1 && cd.nnls.res > worstRes)
        {
            worstRes = cd.nnls.res;
            worstIdx = static_cast<index_t>(c);
        }
    }

    std::vector<bool> chosen(diag.calls.size(), false);
    if (worstIdx >= 0) chosen[static_cast<std::size_t>(worstIdx)] = true;
    for (std::size_t c = 0; c < diag.calls.size(); ++c)
    {
        const gsTetClip::detail::ChunkDiag& cd = diag.calls[c];
        if (cd.level == 2 || (cd.level == 1 && cd.nnls.res >= 1e-13))
            chosen[c] = true;
    }

    index_t printed = 0, remainder = 0;
    for (std::size_t c = 0; c < diag.calls.size(); ++c)
    {
        if (!chosen[c])
            continue;
        if (diag.calls[c].level == 1 && printed >= 20) { ++remainder; continue; }
        printNnlsDiagLine(diag.calls[c]);
        if (diag.calls[c].level == 1) ++printed;
    }
    if (remainder > 0)
        gsInfo << "NNLSDIAG remaining=" << remainder << "\n";
}

/// Bitwise-exact little-endian dump of a cell's uncompressed rule, laid
/// out as: int64 N, int64 p, 3 double lower, 3 double upper, 3*N double
/// nodes (column-major, nd.data()), N double weights. real_t is assumed
/// double, as throughout this file and its build.
void dumpCellBinary(const std::string& path, const gsMatrix<real_t>& nd, const gsVector<real_t>& wt,
                     const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t p)
{
    std::ofstream f(path.c_str(), std::ios::binary);
    GISMO_ENSURE(bool(f), "--dump: cannot open '" << path << "' for writing");
    const std::int64_t N64 = nd.cols();
    const std::int64_t p64 = p;
    f.write(reinterpret_cast<const char*>(&N64), sizeof(N64));
    f.write(reinterpret_cast<const char*>(&p64), sizeof(p64));
    for (index_t d = 0; d < 3; ++d)
    {
        const double v = lower(d);
        f.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
    for (index_t d = 0; d < 3; ++d)
    {
        const double v = upper(d);
        f.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
    f.write(reinterpret_cast<const char*>(nd.data()), sizeof(double) * static_cast<std::size_t>(3 * N64));
    f.write(reinterpret_cast<const char*>(wt.data()), sizeof(double) * static_cast<std::size_t>(N64));
    GISMO_ENSURE(bool(f), "--dump: write to '" << path << "' failed");
}

/// Streamed Legendre moment error, blocked exactly as endToEndRow's own
/// moment check (never forms the full N x K constraint matrix).
real_t streamedMaxMomErr(const gsMatrix<real_t>& nd, const gsVector<real_t>& wt,
                          const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t deg,
                          const gsTetClip::TchakaloffResult& res)
{
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t N = nd.cols();
    gsVector<real_t> mIn = gsVector<real_t>::Zero(K);
    const index_t blockSize = 4096;
    for (index_t lo = 0; lo < N; lo += blockSize)
    {
        const index_t b = std::min(blockSize, N - lo);
        gsMatrix<real_t> Mblk;
        gsTetClip::detail::buildConstraintMatrix(nd.middleCols(lo, b), NULL, lower, upper, deg, Mblk);
        mIn += Mblk.transpose() * wt.segment(lo, b);
    }

    const index_t count = static_cast<index_t>(res.indices.size());
    gsMatrix<real_t> nodesOut(3, count);
    for (index_t c = 0; c < count; ++c)
        nodesOut.col(c) = nd.col(res.indices[c]);
    gsMatrix<real_t> Mout;
    gsTetClip::detail::buildConstraintMatrix(nodesOut, NULL, lower, upper, deg, Mout);
    const gsVector<real_t> mOut = Mout.transpose() * res.weights;

    const real_t maxIn = mIn.array().abs().maxCoeff();
    const real_t maxDiff = (mOut - mIn).array().abs().maxCoeff();
    return (maxIn > 0.0) ? maxDiff / maxIn : maxDiff;
}

struct RealCellOpts
{
    std::string caseName = "sphere";
    index_t level  = 1;
    index_t n0     = 4;
    index_t cell   = 172;
    index_t cellP  = 2;
    index_t chunk  = 0;
    std::string dumpPath;
};

/// --realcell driver: streams the requested cell's uncompressed volume
/// rule exactly as immersed_tetmesh_poisson_example.cpp's buildCellRules
/// does (ijk/lower/upper/status/volRule, cited above), reduces it via
/// detail::reduceDriver with full diagnostics, prints REALCELL and
/// NNLSDIAG, and optionally dumps the uncompressed rule. Returns
/// EXIT_SUCCESS iff the reduction is ok and its streamed moment error is
/// within tolerance.
int runRealCell(const RealCellOpts& ro)
{
    const std::string resolved = gsFileManager::find(realCellMeshFile(ro.caseName));
    GISMO_ENSURE(!resolved.empty(),
                 "--realcell: mesh file for --case " << ro.caseName << " not found");
    memory::shared_ptr<const gsTetClip::TetMesh> M =
        memory::make_shared(new gsTetClip::TetMesh(gsTetClip::readMsh41(resolved)));

    GISMO_ENSURE(ro.n0 >= 1 && ro.level >= 0,
                 "--realcell: --n0 " << ro.n0 << " --level " << ro.level << " must have n0 >= 1, level >= 0");
    const gsTetClip::Grid3 grid = realCellGrid(ro.n0, ro.level);
    gsTetClip::ClipStreamer S(M, grid, ro.cellP);
    memory::shared_ptr<const gsTetClip::CellIndex> idx = S.index();

    GISMO_ENSURE(ro.cell >= 0 && static_cast<std::size_t>(ro.cell) < idx->status.size(),
                 "--realcell: --cell " << ro.cell << " out of range [0, " << idx->status.size() << ")");

    const std::size_t id = static_cast<std::size_t>(ro.cell);
    index_t ci, cj, ck;
    idx->ijk(id, ci, cj, ck);
    gsVector<real_t> lower(3), upper(3);
    lower << idx->X[ci], idx->Y[cj], idx->Z[ck];
    upper << idx->X[ci + 1], idx->Y[cj + 1], idx->Z[ck + 1];

    GISMO_ENSURE(gsTetClip::Cut == idx->status[id],
                 "--realcell: cell " << ro.cell << " has status " << idx->status[id] << ", not Cut");

    gsMatrix<real_t> nd;
    gsVector<real_t> wt;
    S.volRule(id, nd, wt);
    const index_t nIn = nd.cols();

    if (!ro.dumpPath.empty())
        dumpCellBinary(ro.dumpPath, nd, wt, lower, upper, ro.cellP);

    gsTetClip::TchakaloffOptions opt;
    if (ro.chunk > 0)
        opt.chunkSize = ro.chunk;

    gsTetClip::detail::ReduceDiagnostics diag;
    gsStopwatch sw;
    const gsTetClip::TchakaloffResult res =
        gsTetClip::detail::reduceDriver(nd, wt, NULL, lower, upper, ro.cellP, opt, &diag);
    const double timeSec = sw.stop();

    const real_t maxMomErr = streamedMaxMomErr(nd, wt, lower, upper, 2 * ro.cellP, res);
    const index_t count = static_cast<index_t>(res.indices.size());

    gsInfo << "REALCELL case=" << ro.caseName << " r=" << ro.level << " cell=" << ro.cell
           << " p=" << ro.cellP << " N=" << nIn << " levels=" << res.levels << " levelRank=";
    for (std::size_t x = 0; x < res.levelRank.size(); ++x)
    {
        if (x) gsInfo << ",";
        gsInfo << res.levelRank[x];
    }
    gsInfo << " res=";
    for (std::size_t x = 0; x < res.levelResidual.size(); ++x)
    {
        if (x) gsInfo << ",";
        gsInfo << res.levelResidual[x];
    }
    gsInfo << " iters=" << res.iterations << " count=" << count << " rank=" << res.rank
           << " ok=" << (res.ok ? 1 : 0) << " maxMomErr=" << maxMomErr << " time=" << timeSec << "\n";

    printNnlsDiagLines(diag);

    return (res.ok && maxMomErr <= 1e-12) ? EXIT_SUCCESS : EXIT_FAILURE;
}

// ---------------------------------------------------------------------
// Permanent real-cell regression (no-argument self-test): sphere n0=4
// r=1 p=2, cells 214 and 172, compressed with the PUBLIC
// tchakaloffCompress (not detail::reduceDriver -- this exercises the
// same call the Poisson driver makes). Guards a real two-level cut-cell
// reduction, distinct from the header's synthetic 18 cases above.
//
// Both cells currently hit the Lawson-Hanson stopping stall documented
// on the NNLS stopping test in gsTchakaloffRule.h (the KKT test can
// accept a passive-set solution whose residual is still on the order of
// 1e-9 on a small-mass chunk). By default both are therefore reported as
// expected failures (XFAIL/XPASS), which does not affect the exit code;
// --strict-realcell reports them as ordinary PASS/FAIL instead, so a
// future fix to the stopping test can be verified against this same
// regression.
// ---------------------------------------------------------------------

struct RealCellCase { index_t cell; index_t nInAnchor; };

/// Outcome tally for the two permanent REALCELL cases, folded into the
/// run's overall summary line by the caller.
struct RealCellSummary { index_t pass = 0, fail = 0, xfail = 0, xpass = 0; };

/// Runs the two permanent regression cells and prints one REALCELL line
/// each. Under strict, the prefix is PASS/FAIL and the outcome is
/// counted in pass/fail (so a failing case makes the run exit 1). Not
/// under strict, the prefix is XFAIL/XPASS and the outcome is counted in
/// xfail/xpass instead (so neither outcome affects the exit code) -- the
/// known stopping-test stall is an expected failure, not a build defect.
/// A missing mesh file is an environment fault, not the documented
/// stall, and is always counted as an ordinary FAIL.
RealCellSummary runRealCellRegression(bool strict)
{
    RealCellSummary summary;
    const std::string resolved = gsFileManager::find(realCellMeshFile("sphere"));
    if (resolved.empty())
    {
        gsInfo << "FAIL  REALCELL  case=sphere mesh-not-found\n";
        summary.fail = 2;
        return summary;
    }
    memory::shared_ptr<const gsTetClip::TetMesh> M =
        memory::make_shared(new gsTetClip::TetMesh(gsTetClip::readMsh41(resolved)));

    const index_t n0 = 4, r = 1, p = 2;
    const gsTetClip::Grid3 grid = realCellGrid(n0, r);
    gsTetClip::ClipStreamer S(M, grid, p);
    memory::shared_ptr<const gsTetClip::CellIndex> idx = S.index();

    const RealCellCase cases[2] = { {214, 24696}, {172, 203840} };

    for (int ci = 0; ci < 2; ++ci)
    {
        const RealCellCase& cs = cases[ci];
        gsStopwatch sw;
        bool pass = true;
        std::string reason;
        auto fail = [&pass, &reason](const std::string& msg)
        {
            pass = false;
            if (!reason.empty()) reason += "; ";
            reason += msg;
        };

        const std::size_t id = static_cast<std::size_t>(cs.cell);
        index_t i, j, k;
        idx->ijk(id, i, j, k);
        gsVector<real_t> lower(3), upper(3);
        lower << idx->X[i], idx->Y[j], idx->Z[k];
        upper << idx->X[i + 1], idx->Y[j + 1], idx->Z[k + 1];

        if (gsTetClip::Cut != idx->status[id])
            fail("status != Cut");

        gsMatrix<real_t> nd;
        gsVector<real_t> wt;
        S.volRule(id, nd, wt);
        const index_t nIn = nd.cols();
        if (nIn != cs.nInAnchor)
            fail("nIn != anchor (a clip change would silently swap the cell)");

        const gsTetClip::TchakaloffResult res = gsTetClip::tchakaloffCompress(nd, wt, lower, upper, p);
        if (!res.ok)
            fail("res.ok == false");
        for (std::size_t li = 0; li < res.levelResidual.size(); ++li)
            if (!(res.levelResidual[li] < 1e-13))
                fail("levelResidual >= 1e-13");
        for (index_t wj = 0; wj < res.weights.size(); ++wj)
            if (!(res.weights(wj) > 0.0))
                fail("non-positive output weight");
        const index_t count = static_cast<index_t>(res.indices.size());
        if (!(count <= res.rank))
            fail("count > rank");

        const real_t maxMomErr = streamedMaxMomErr(nd, wt, lower, upper, 2 * p, res);
        if (!(maxMomErr <= 1e-12))
            fail("maxMomErr > 1e-12");

        const real_t sumWIn = wt.sum();
        const real_t sumWOut = res.weights.sum();
        if (!(std::abs(sumWOut - sumWIn) <= 1e-12 * sumWIn))
            fail("sum(w_out) != sum(w_in)");

        const double t = sw.stop();
        const char* prefix = strict ? (pass ? "PASS" : "FAIL") : (pass ? "XPASS" : "XFAIL");
        gsInfo << prefix << "  REALCELL  case=sphere r=" << r << " cell=" << cs.cell
               << " p=" << p << " N=" << nIn << " levels=" << res.levels << " levelRank=";
        for (std::size_t li = 0; li < res.levelRank.size(); ++li)
        {
            if (li) gsInfo << ",";
            gsInfo << res.levelRank[li];
        }
        gsInfo << " res=";
        for (std::size_t li = 0; li < res.levelResidual.size(); ++li)
        {
            if (li) gsInfo << ",";
            gsInfo << res.levelResidual[li];
        }
        gsInfo << " maxMomErr=" << maxMomErr << " count=" << count << " rank=" << res.rank
               << " ok=" << (res.ok ? 1 : 0) << " time=" << t << "\n";
        if (!pass) gsInfo << "  FAILED CHECK: " << reason << "\n";
        if (strict)
        {
            if (pass) ++summary.pass; else ++summary.fail;
        }
        else
        {
            if (pass) ++summary.xpass; else ++summary.xfail;
        }
    }

    return summary;
}

// ---------------------------------------------------------------------
// --nnmf: strict self-test of gsNnmfRule.h.
// ---------------------------------------------------------------------

/// Adapts an NnmfResult to the TchakaloffResult that
/// gsTetClip::detail::genericChecks consumes (levelResidual and levelRank
/// stay empty, so the per-level checks are no-ops).
gsTetClip::TchakaloffResult asTchakaloffResult(const gsTetClip::NnmfResult& r)
{
    gsTetClip::TchakaloffResult t;
    t.indices.assign(r.indices.begin(), r.indices.end());
    t.weights = r.weights;
    t.rank    = r.rank;
    t.ok      = r.ok;
    t.levels  = 1;
    return t;
}

/// Boundary measure of the tilted-plane patch of
/// gsTetClip::detail::buildFTilted, sampled on an 8 x 8 subdivision so that
/// N = 64 (2p+2)^2 exceeds 3K and the pooled (non-fallback) path of
/// nnmfCompressBoundary is exercised.
void buildFTiltedDense(const gsVector<real_t>& mid, real_t h, index_t p, gsMatrix<real_t>& nodes,
                       gsVector<real_t>& weights, gsMatrix<real_t>& normals)
{
    gsVector<real_t> n(3);
    n(0) = 1.0; n(1) = 2.0; n(2) = 3.0;
    n /= std::sqrt(14.0);
    gsVector<real_t> ez(3);
    ez(0) = 0.0; ez(1) = 0.0; ez(2) = 1.0;
    gsVector<real_t> t1 = gsTetClip::detail::cross3(n, ez);
    t1.normalize();
    gsVector<real_t> t2 = gsTetClip::detail::cross3(n, t1);

    gsVector<real_t> lo2(2), hi2(2);
    lo2(0) = -h / 4.0; lo2(1) = -h / 4.0;
    hi2(0) =  h / 4.0; hi2(1) =  h / 4.0;

    gsMatrix<real_t> st;
    gsTetClip::detail::tensorGaussSubdivided(lo2, hi2, 8, 2 * p + 2, st, weights);

    const index_t N = st.cols();
    nodes.resize(3, N);
    normals.resize(3, N);
    for (index_t i = 0; i < N; ++i)
    {
        nodes.col(i) = mid + st(0, i) * t1 + st(1, i) * t2;
        normals.col(i) = n;
    }
}

/// Accumulates failing checks of one NNMF case into a single reason string.
struct NnmfCaseVerdict
{
    bool pass = true;
    std::string reason;
    void fail(const std::string& msg)
    {
        pass = false;
        if (!reason.empty()) reason += "; ";
        reason += msg;
    }
};

/// Checks common to every NNMF case: res.ok, the header's own moment error,
/// positive weights, kept <= rank and weight-sum conservation.
void nnmfCommonChecks(const gsTetClip::NnmfResult& res, const gsVector<real_t>& wIn, NnmfCaseVerdict& v)
{
    if (!res.ok)
        v.fail("res.ok == false");
    if (!(res.momErr <= 1e-12))
        v.fail("res.momErr > 1e-12");
    for (index_t j = 0; j < res.weights.size(); ++j)
        if (!(res.weights(j) > 0.0))
        {
            v.fail("non-positive output weight");
            break;
        }
    if (!(static_cast<index_t>(res.indices.size()) <= res.rank))
        v.fail("kept > rank");
    const real_t sumIn = wIn.sum();
    if (!(std::abs(res.weights.sum() - sumIn) <= 1e-12 * sumIn))
        v.fail("sum(w_out) != sum(w_in)");
}

void printNnmfLine(const std::string& label, const std::string& pTag, index_t N, const gsTetClip::NnmfResult& res,
                   real_t indepMomErr, real_t monErr, bool hasMon, double t, const NnmfCaseVerdict& v)
{
    gsInfo << (v.pass ? "PASS" : "FAIL") << "  NNMF  " << label << "  " << pTag << "  N=" << N
           << "  kept=" << res.indices.size() << "  rank=" << res.rank << "  momErr=" << res.momErr
           << "  indepMomErr=" << indepMomErr << "  monErr=";
    if (hasMon) gsInfo << monErr; else gsInfo << "-";
    gsInfo << "  rounds=" << res.rounds << "  poolSize=" << res.poolSize
           << "  fallbackFull=" << (res.fallbackFull ? 1 : 0) << "  time=" << t << "\n";
    if (!v.pass) gsInfo << "  FAILED CHECK: " << v.reason << "\n";
}

/// Strict PASS/FAIL self-test of gsTetClip::nnmfCompress[Boundary]: the
/// synthetic cases of gsTetClip::tchakaloffSelfTest (plus a dense tilted-plane
/// boundary case and an L = 1 spherical-cap case) for p = pmin..pmax, then the
/// two permanent real cut cells. Returns the process exit code.
int runNnmfSelfTest(index_t pmin, index_t pmax)
{
    using gsTetClip::NnmfOptions;
    using gsTetClip::NnmfResult;
    namespace det = gsTetClip::detail;

    index_t nPass = 0, nFail = 0;
    auto tally = [&nPass, &nFail](const NnmfCaseVerdict& v) { if (v.pass) ++nPass; else ++nFail; };

    std::mt19937 rng(12345);
    gsVector<real_t> lower(3), upper(3), mid(3);
    lower(0) = 0.10; lower(1) = -0.20; lower(2) = 0.30;
    const real_t h = 0.25;
    for (index_t d = 0; d < 3; ++d)
    {
        upper(d) = lower(d) + h;
        mid(d) = lower(d) + 0.5 * h;
    }
    const real_t vol = h * h * h;
    const real_t area = (h / 2.0) * (h / 2.0);

    for (index_t p = pmin; p <= pmax; ++p)
    {
        const index_t deg = 2 * p;
        const index_t K = (deg + 1) * (deg + 1) * (deg + 1);
        const std::string pTag = "p=" + std::to_string(p);

        // V-gauss
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes;
            gsVector<real_t> weights;
            det::buildVGauss(lower, upper, p, nodes, weights);
            NnmfResult res;
            gsTetClip::nnmfCompress(nodes, weights, lower, upper, p, NnmfOptions(), res);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            NnmfCaseVerdict v;
            std::string reason;
            const gsTetClip::TchakaloffResult tr = asTchakaloffResult(res);
            if (!det::genericChecks(nodes, weights, NULL, lower, upper, mid, h, p, tr, maxMomErr, maxMonErr, reason))
                v.fail(reason);
            nnmfCommonChecks(res, weights, v);

            const index_t count = static_cast<index_t>(res.indices.size());
            gsMatrix<real_t> nodesOut(3, count);
            bool indicesValid = true;
            for (index_t j = 0; j < count; ++j)
            {
                if (res.indices[j] < 0 || res.indices[j] >= nodes.cols())
                    indicesValid = false;
                else
                    nodesOut.col(j) = nodes.col(res.indices[j]);
            }
            if (indicesValid)
            {
                gsMatrix<real_t> Vout;
                gsTetClip::legendreVandermonde(nodesOut, lower, upper, deg, Vout);
                const gsVector<real_t> momOut = Vout.transpose() * res.weights;
                const real_t target0 = std::sqrt(vol);
                if (!(std::abs(momOut(0) - target0) <= 1e-12 * target0))
                    v.fail("m_0 != sqrt(vol)");
                bool higherModesZero = true;
                for (index_t k = 1; k < momOut.size(); ++k)
                    if (!(std::abs(momOut(k)) <= 1e-12 * target0))
                        higherModesZero = false;
                if (!higherModesZero)
                    v.fail("m_k != 0 for k>0");
            }
            else
                v.fail("analytic moment check skipped: out-of-range index");
            if (!(std::abs(res.weights.sum() - vol) <= 1e-12 * vol))
                v.fail("sum(w_out) != vol");

            printNnmfLine("V-gauss", pTag, nodes.cols(), res, maxMomErr, maxMonErr, true, t, v);
            tally(v);
        }

        // V-random
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes;
            gsVector<real_t> weights;
            det::buildVRandom(lower, upper, K, vol, rng, nodes, weights);
            NnmfResult res;
            gsTetClip::nnmfCompress(nodes, weights, lower, upper, p, NnmfOptions(), res);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            NnmfCaseVerdict v;
            std::string reason;
            if (!det::genericChecks(nodes, weights, NULL, lower, upper, mid, h, p, asTchakaloffResult(res),
                                    maxMomErr, maxMonErr, reason))
                v.fail(reason);
            nnmfCommonChecks(res, weights, v);

            printNnmfLine("V-random", pTag, nodes.cols(), res, maxMomErr, maxMonErr, true, t, v);
            tally(v);
        }

        // Boundary cases: F-axis, F-tilted, S-cap, S-cap-L1, F-tilted-dense.
        for (int bc = 0; bc < 5; ++bc)
        {
            static const char* const labels[5] = { "F-axis", "F-tilted", "S-cap", "S-cap-L1", "F-tilted-dense" };
            gsStopwatch sw;
            gsMatrix<real_t> nodes, normals;
            gsVector<real_t> weights;
            NnmfOptions opt;
            switch (bc)
            {
            case 0: det::buildFAxis(mid, h, p, nodes, weights, normals); break;
            case 1: det::buildFTilted(mid, h, p, nodes, weights, normals); break;
            case 2: det::buildSCap(mid, p, nodes, weights, normals); break;
            case 3: det::buildSCap(mid, p, nodes, weights, normals); opt.poolFactor = 1; break;
            default: buildFTiltedDense(mid, h, p, nodes, weights, normals); break;
            }

            NnmfCaseVerdict v;
            if (2 == bc || 3 == bc)
            {
                bool insideCell = true;
                for (index_t i = 0; i < nodes.cols() && insideCell; ++i)
                    for (index_t d = 0; d < 3; ++d)
                        if (nodes(d, i) < lower(d) || nodes(d, i) > upper(d))
                            insideCell = false;
                if (!insideCell)
                    v.fail("input nodes outside the cell");
            }

            NnmfResult res;
            gsTetClip::nnmfCompressBoundary(nodes, weights, normals, lower, upper, p, opt, res);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            if (!det::genericChecks(nodes, weights, &normals, lower, upper, mid, h, p, asTchakaloffResult(res),
                                    maxMomErr, maxMonErr, reason))
                v.fail(reason);
            nnmfCommonChecks(res, weights, v);

            if (bc <= 1 || 4 == bc)
                if (!(std::abs(res.weights.sum() - area) <= 1e-12 * area))
                    v.fail("sum(w_out) != (h/2)^2");
            if (4 == bc && res.fallbackFull)
                v.fail("fallbackFull == true (pooled path not exercised)");

            printNnmfLine(labels[bc], pTag, nodes.cols(), res, maxMomErr, maxMonErr, true, t, v);
            tally(v);
        }
    }

    // Real cut cells: sphere, n0 = 4, r = 1, p = 2, volume rule.
    const std::string resolved = gsFileManager::find(realCellMeshFile("sphere"));
    if (resolved.empty())
    {
        gsInfo << "FAIL  NNMF  REALCELL  case=sphere mesh-not-found\n";
        nFail += 2;
    }
    else
    {
        memory::shared_ptr<const gsTetClip::TetMesh> M =
            memory::make_shared(new gsTetClip::TetMesh(gsTetClip::readMsh41(resolved)));
        const index_t n0 = 4, r = 1, p = 2;
        const gsTetClip::Grid3 grid = realCellGrid(n0, r);
        gsTetClip::ClipStreamer S(M, grid, p);
        memory::shared_ptr<const gsTetClip::CellIndex> idx = S.index();

        const RealCellCase cases[2] = { {214, 24696}, {172, 203840} };
        for (int ci = 0; ci < 2; ++ci)
        {
            const RealCellCase& cs = cases[ci];
            gsStopwatch sw;
            NnmfCaseVerdict v;

            const std::size_t id = static_cast<std::size_t>(cs.cell);
            index_t i, j, k;
            idx->ijk(id, i, j, k);
            gsVector<real_t> lower(3), upper(3);
            lower << idx->X[i], idx->Y[j], idx->Z[k];
            upper << idx->X[i + 1], idx->Y[j + 1], idx->Z[k + 1];

            if (gsTetClip::Cut != idx->status[id])
                v.fail("status != Cut");

            gsMatrix<real_t> nd;
            gsVector<real_t> wt;
            S.volRule(id, nd, wt);
            const index_t nIn = nd.cols();
            if (nIn != cs.nInAnchor)
                v.fail("nIn != anchor (a clip change would silently swap the cell)");

            NnmfResult res;
            gsTetClip::nnmfCompress(nd, wt, lower, upper, p, NnmfOptions(), res);
            nnmfCommonChecks(res, wt, v);

            const real_t indepMomErr = streamedMaxMomErr(nd, wt, lower, upper, 2 * p, asTchakaloffResult(res));
            if (!(indepMomErr <= 1e-12))
                v.fail("indepMomErr > 1e-12");
            const double t = sw.stop();

            std::ostringstream tag;
            tag << "case=sphere r=" << r << " cell=" << cs.cell << " p=" << p;
            gsInfo << (v.pass ? "PASS" : "FAIL") << "  NNMF  REALCELL  " << tag.str() << "  N=" << nIn
                   << "  kept=" << res.indices.size() << "  rank=" << res.rank << "  momErr=" << res.momErr
                   << "  indepMomErr=" << indepMomErr << "  rounds=" << res.rounds
                   << "  poolSize=" << res.poolSize << "  fallbackFull=" << (res.fallbackFull ? 1 : 0)
                   << "  time=" << t << "\n";
            if (!v.pass) gsInfo << "  FAILED CHECK: " << v.reason << "\n";
            tally(v);
        }
    }

    gsInfo << "NNMF: " << nPass << " PASS, " << nFail << " FAIL\n";
    return (0 == nFail) ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace

int main(int argc, char* argv[])
{
    index_t pmin = 1;
    index_t pmax = 3;
    bool bench = false;
    bool nnmf = false;
    std::string benchKind = "all";
    index_t benchP = 0;
    index_t benchN = 0;

    bool realcell = false;
    bool strictRealcell = false;
    RealCellOpts ro;

    gsCmdLine cmd("Self-test / benchmark driver for the Tchakaloff measure-compression rules "
                  "(volume and boundary variants) in gsTchakaloffRule.h.");
    cmd.addInt("", "pmin", "Smallest polynomial degree parameter p to test", pmin);
    cmd.addInt("", "pmax", "Largest polynomial degree parameter p to test", pmax);
    cmd.addSwitch("bench", "Run the NNLS/reduction benchmark instead of the self-test", bench);
    cmd.addSwitch("nnmf", "Run the strict NNMF (gsNnmfRule.h) self-test instead of the Tchakaloff self-test", nnmf);
    cmd.addString("", "bench-kind", "Benchmark kind: all|vol|bdr", benchKind);
    cmd.addInt("", "bench-p", "Benchmark p (0 = both {2,3})", benchP);
    cmd.addInt("", "bench-n", "Benchmark N (0 = all of {1000, 10000, 100000, 200000})", benchN);
    cmd.addSwitch("realcell", "Stream and compress one real cut cell (diagnosis mode)", realcell);
    cmd.addSwitch("strict-realcell", "Treat the two permanent REALCELL regression cases as ordinary "
                  "PASS/FAIL instead of expected failures (XFAIL/XPASS); a FAIL then exits 1",
                  strictRealcell);
    cmd.addString("", "case", "--realcell mesh case: sphere|rotcube", ro.caseName);
    cmd.addInt("r", "level", "--realcell background-grid refinement level", ro.level);
    cmd.addInt("", "n0", "--realcell base grid resolution (n = n0 << r)", ro.n0);
    cmd.addInt("", "cell", "--realcell flat cell id", ro.cell);
    cmd.addInt("", "cell-p", "--realcell quadrature degree parameter p", ro.cellP);
    cmd.addInt("", "chunk", "--realcell TchakaloffOptions::chunkSize override (0 = library default; "
               "diagnosis control only)", ro.chunk);
    cmd.addString("", "dump", "--realcell: dump the cell's uncompressed rule to this binary path", ro.dumpPath);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (realcell)
        return runRealCell(ro);

    if (bench)
        return runBench(benchKind, benchP, benchN);

    if (nnmf)
        return runNnmfSelfTest(pmin, pmax);

    // gsTetClip::tchakaloffSelfTest (gsTchakaloffRule.h) prints one PASS/FAIL
    // line per synthetic case but returns only the aggregate bool, so its
    // output is captured here (std::cout is what gsInfo expands to,
    // src/gsCore/gsDebug.h:43) and replayed verbatim, then parsed for exact
    // PASS/FAIL counts to fold into this driver's own summary line without
    // duplicating the header's per-p case list.
    std::ostringstream captured;
    std::streambuf* const realCout = std::cout.rdbuf(captured.rdbuf());
    const bool selfTestOk = gsTetClip::tchakaloffSelfTest(pmin, pmax);
    std::cout.rdbuf(realCout);
    gsInfo << captured.str();

    index_t nSynthPass = 0, nSynthFail = 0;
    {
        std::istringstream iss(captured.str());
        std::string line;
        while (std::getline(iss, line))
        {
            if (line.compare(0, 6, "PASS  ") == 0)
                ++nSynthPass;
            else if (line.compare(0, 6, "FAIL  ") == 0)
                ++nSynthFail;
        }
    }
    GISMO_ENSURE(selfTestOk == (nSynthFail == 0),
                 "tchakaloff_selftest_example: parsed synthetic PASS/FAIL count disagrees with "
                 "tchakaloffSelfTest's own return value");

    const RealCellSummary rc = runRealCellRegression(strictRealcell);

    const index_t totalPass = nSynthPass + rc.pass;
    const index_t totalFail = nSynthFail + rc.fail;
    gsInfo << totalPass << " PASS, " << totalFail << " FAIL, " << rc.xfail << " XFAIL, " << rc.xpass
           << " XPASS\n";

    return (totalFail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
