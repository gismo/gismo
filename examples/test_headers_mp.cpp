/** @file test_headers_mp.cpp

    @brief Self-test for the multi-patch AS-G1 header API
    (`examples/asg1_mp/*.h`).

    By default it builds an N-patch strip of unit squares (no interior
    vertex) and verifies AS-G1 smoothness across every interface -- this
    exercises the Interior + Edge spaces (Phase B).  A geometry file may
    also be supplied with -f.

    Usage:
        ./bin/test_headers_mp [-N 4] [-r refinements] [-d minDeg]
        ./bin/test_headers_mp -f planar/lshape2d_3patches_tens.xml
*/

#include <gismo.h>
#include "asg1_mp/asg1_multipatch_mp.h"

using namespace gismo;

// Build a 1xN strip of unit squares: patch k = [k,k+1] x [0,1].
static gsMultiPatch<real_t> makeStrip(index_t N)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 0, 2);   // degree 1
    gsMultiPatch<real_t> mp;
    for (index_t k = 0; k < N; ++k)
    {
        gsMatrix<real_t> coefs(4, 2);
        coefs << real_t(k),   0.0,
                 real_t(k+1), 0.0,
                 real_t(k),   1.0,
                 real_t(k+1), 1.0;
        mp.addPatch(gsTensorBSpline<2,real_t>(kv, kv, coefs));
    }
    mp.computeTopology();
    return mp;
}

int main(int argc, char* argv[])
{
    std::string fn;
    index_t N = 4, refs = 1, minDeg = 3, nCheck = 21;
    bool noTrunc = false, noVtx = false;
    index_t vWin = 4, vSpan = 2, vSamp = 24;
    real_t vTol = 1e-7;
    bool dumpEdges = false, noCheck = false;

    gsCmdLine cmd("Test multi-patch AS-G1 header API.");
    cmd.addString("f", "file",   "Geometry XML file (overrides strip)", fn);
    cmd.addInt   ("N", "npatch", "Number of strip patches (if no file)", N);
    cmd.addInt   ("r", "refs",   "Uniform refinements", refs);
    cmd.addInt   ("d", "degree", "Minimum spline degree", minDeg);
    cmd.addInt   ("n", "ncheck", "Samples per interface", nCheck);
    cmd.addSwitch("notrunc", "Disable Phase D edge-end truncation", noTrunc);
    cmd.addSwitch("novtx",   "Disable Phase E vertex space", noVtx);
    cmd.addSwitch("dump",    "Dump per-edge gluing data", dumpEdges);
    cmd.addSwitch("nocheck", "Skip smoothness + reproduction (time build only)", noCheck);
    cmd.addInt   ("w", "vwin",  "Vertex corner window size", vWin);
    cmd.addInt   ("s", "vspan", "Vertex tangent spans covered", vSpan);
    cmd.addInt   ("m", "vsamp", "Vertex samples per edge", vSamp);
    cmd.addReal  ("t", "vtol",  "Vertex null-space rel. tolerance", vTol);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> mp;
    if (fn.empty())
    {
        mp = makeStrip(N);
        gsInfo << "Geometry    : " << N << "-patch unit-square strip\n";
    }
    else
    {
        gsReadFile<real_t>(fn, mp);
        if (mp.nPatches() == 0) { gsInfo << "Cannot read " << fn << "\n"; return 1; }
        gsInfo << "Geometry    : " << fn << "\n";
    }

    asg1mp::AsG1MpOptions<real_t> opts;
    opts.refinements = refs;
    opts.minDegree   = static_cast<short_t>(minDeg);
    opts.truncate    = !noTrunc;
    opts.vertex      = !noVtx;
    opts.vtxWindow   = vWin;
    opts.vtxSpanCover = vSpan;
    opts.vtxSamples  = vSamp;
    opts.vtxNullTol  = vTol;
    opts.verbose     = false;

    gsStopwatch clock;
    asg1mp::AsG1Multipatch<real_t> B =
        asg1mp::buildAsG1Multipatch<real_t>(mp, opts);
    const double buildTime = clock.stop();

    gsInfo << "nPatches    : " << B.mp.nPatches() << "\n";
    gsInfo << "nInterfaces : " << B.mp.nInterfaces() << "\n";
    gsInfo << "buildTime   : " << buildTime << " s\n";
    gsInfo << "nGlobal     : " << B.nGlobal
           << "  (interior=" << B.edgeBlockOff
           << ", edgeBlock=" << (B.vertexBlockOff - B.edgeBlockOff)
           << ", vertexBlock=" << (B.nGlobal - B.vertexBlockOff) << ")\n";
    for (size_t v = 0; v < B.vertexFuncs.size(); ++v)
        gsInfo << "  vertex block " << v << " : id=" << B.vertexFuncs[v].vertexId
               << "  funcs=" << B.vertexFuncs[v].W.cols() << "\n";
    for (size_t p = 0; p < B.G.size(); ++p)
        gsInfo << "  G[" << p << "] : " << B.G[p].rows()
               << " x " << B.G[p].cols()
               << "  interiorDofs=" << B.interiorDofs[p].size() << "\n";

    for (size_t e = 0; e < B.edges.size(); ++e)
    {
        const asg1mp::EdgeInfo<real_t>& ed = B.edges[e];
        if (ed.isBoundary) continue;
        if (!dumpEdges) continue;
        gsInfo << "  edge " << ed.id << " p" << ed.side1.patch << "s" << ed.side1.side()
               << " <-> p" << ed.side2.patch << "s" << ed.side2.side()
               << (ed.flipped ? " [flip] " : " ") << "gd=[" << ed.gluingData.transpose() << "]\n";
    }

    clock.restart();
    asg1mp::G1ReportMP<real_t> R;
    double checkTime = 0;
    if (!noCheck)
    {
        R = asg1mp::smoothnessCheckMP<real_t>(B, nCheck);
        checkTime = clock.stop();
    }

    gsInfo << "\n=== AS-G1 smoothness ===\n";
    gsInfo << "  checkTime  : " << checkTime << " s\n";
    gsInfo << "  maxValErr  : " << R.maxValErr  << "\n";
    gsInfo << "  maxGradErr : " << R.maxGradErr << "\n";
    for (size_t i = 0; i < R.perEdgeGradErr.size(); ++i)
        gsInfo << "    edge " << i << " gradErr : " << R.perEdgeGradErr[i] << "\n";
    gsInfo << "  STATUS     : " << (noCheck ? "skipped" : (R.pass ? "PASS" : "FAIL")) << "\n";

    // ---- Polynomial reproduction check ----
    // Can the assembled global space represent a given global polynomial
    // f(x,y) exactly?  Residual ~0 means the space contains f.  Uses a
    // sparse normal-equation solve so it scales to fine meshes.
    if (!noCheck)
    {
        const gsMultiPatch<real_t>& mp2 = B.mp;
        const size_t P = mp2.nPatches();
        index_t totRows = 0;
        std::vector<index_t> roffv(P+1, 0);
        for (size_t p = 0; p < P; ++p)
        {
            roffv[p] = totRows;
            totRows += mp2.patch(p).basis().size();
        }
        roffv[P] = totRows;

        // Build the stacked sparse G (totRows x nGlobal) once.
        gsSparseMatrix<real_t> Gst(totRows, B.nGlobal);
        {
            gsSparseEntries<real_t> ent;
            for (size_t p = 0; p < P; ++p)
                for (index_t cc = 0; cc < B.G[p].outerSize(); ++cc)
                    for (gsSparseMatrix<real_t>::InnerIterator it(B.G[p], cc); it; ++it)
                        ent.add(roffv[p] + it.row(), it.col(), it.value());
            Gst.setFrom(ent); Gst.makeCompressed();
        }
        gsSparseMatrix<real_t> A = Gst.transpose() * Gst;
        for (index_t i = 0; i < A.rows(); ++i) A.coeffRef(i,i) += 1e-12;
        A.makeCompressed();
        typename gsSparseSolver<real_t>::SimplicialLDLT solver(A);

        auto reproduce = [&](std::function<real_t(real_t,real_t)> f)->real_t
        {
            gsMatrix<real_t> c(totRows, 1);
            for (size_t p = 0; p < P; ++p)
            {
                const gsBasis<real_t>& tb = mp2.patch(p).basis();
                const index_t n = tb.size();
                gsMatrix<real_t> anchors = tb.anchors();
                gsMatrix<real_t> phys = mp2.patch(p).eval(anchors);
                gsMatrix<real_t> fval(1, n);
                for (index_t k = 0; k < n; ++k)
                    fval(0,k) = f(phys(0,k), phys(1,k));
                typename gsGeometry<real_t>::uPtr geo = tb.interpolateAtAnchors(fval);
                c.block(roffv[p],0,n,1) = geo->coefs();
            }
            gsMatrix<real_t> g = solver.solve(Gst.transpose() * c);
            const real_t den = std::max<real_t>(c.norm(), 1e-30);
            return (Gst*g - c).norm() / den;
        };

        gsInfo << "\n=== Polynomial reproduction (rel. residual) ===\n";
        gsInfo << "  1    : " << reproduce([](real_t  , real_t  ){ return 1.0; }) << "\n";
        gsInfo << "  x    : " << reproduce([](real_t x, real_t  ){ return x;   }) << "\n";
        gsInfo << "  y    : " << reproduce([](real_t  , real_t y){ return y;   }) << "\n";
        gsInfo << "  x^2  : " << reproduce([](real_t x, real_t  ){ return x*x; }) << "\n";
        gsInfo << "  x*y  : " << reproduce([](real_t x, real_t y){ return x*y; }) << "\n";
        gsInfo << "  y^2  : " << reproduce([](real_t  , real_t y){ return y*y; }) << "\n";
    }

    // ---- gsMappedBasis export ----
    {
        gsMappedBasis<2,real_t> mb;
        asg1mp::toMappedBasis<real_t>(B, mb);
        gsInfo << "\n=== gsMappedBasis export ===\n";
        gsInfo << "  globalSize : " << mb.globalSize()
               << "  (expected " << B.nGlobal << ")\n";
        gsInfo << "  localSize  : " << mb.localSize() << "\n";
        gsInfo << "  status     : "
               << (mb.globalSize() == B.nGlobal ? "OK" : "MISMATCH") << "\n";
    }


    return R.pass ? 0 : 2;
}
