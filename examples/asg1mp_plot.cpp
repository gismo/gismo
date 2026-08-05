/** @file asg1mp_plot.cpp

    @brief Plot AS-G1 multi-patch global basis functions to Paraview.

    For a chosen geometry it builds the AS-G1 basis and writes, for
    selected global DOFs, the corresponding field (the basis function
    evaluated on every patch) as a Paraview file.  Functions are grouped
    and named by category: vertex, edge (trace/deriv) and interior.

    Usage:
      ./bin/asg1mp_plot -f <file> -r <refs> -o <outDir> \
            --which both --maxedge 3 --npts 3000
        --which : vertex | edge | both | all
        --maxedge : how many trace and deriv functions per interface
*/

#include <gismo.h>
#include "asg1_mp/asg1_multipatch_mp.h"

using namespace gismo;

static void writeField(const asg1mp::AsG1Multipatch<real_t>& B, index_t idx,
                       const std::string& fname, index_t npts)
{
    gsVector<real_t> e = gsVector<real_t>::Zero(B.nGlobal);
    e(idx) = 1.0;
    gsMultiPatch<real_t> sol;
    for (size_t p = 0; p < B.mp.nPatches(); ++p)
    {
        const gsBasis<real_t>& tb = B.mp.patch(p).basis();
        gsVector<real_t> c = B.G[p] * e;
        gsMatrix<real_t> coefs = c;
        sol.addPatch(tb.makeGeometry(give(coefs)));
    }
    gsField<real_t> field(B.mp, sol);
    gsWriteParaview<real_t>(field, fname, npts);
}

int main(int argc, char* argv[])
{
    std::string fn, outDir("plots"), which("both");
    index_t refs = 2, minDeg = 3, maxEdge = 3, npts = 3000;

    gsCmdLine cmd("Plot AS-G1 multi-patch basis functions to Paraview.");
    cmd.addString("f", "file",   "Geometry XML file", fn);
    cmd.addString("o", "outDir", "Output directory", outDir);
    cmd.addString("",  "which",  "vertex | edge | both | all", which);
    cmd.addInt   ("r", "refs",   "Uniform refinements", refs);
    cmd.addInt   ("d", "degree", "Minimum spline degree", minDeg);
    cmd.addInt   ("",  "maxedge","Trace & deriv funcs per interface", maxEdge);
    cmd.addInt   ("",  "npts",   "Paraview sample points per patch", npts);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>(fn, mp);
    if (mp.nPatches() == 0) { gsInfo << "Cannot read " << fn << "\n"; return 1; }
    mp.computeTopology();

    if (outDir.back() != '/') outDir += '/';
    gsFileManager::mkdir(outDir);

    asg1mp::AsG1MpOptions<real_t> opts;
    opts.refinements = refs;
    opts.minDegree   = static_cast<short_t>(minDeg);
    opts.verbose     = false;
    asg1mp::AsG1Multipatch<real_t> B = asg1mp::buildAsG1Multipatch<real_t>(mp, opts);

    const bool doVtx  = (which=="vertex"||which=="both"||which=="all");
    const bool doEdge = (which=="edge"  ||which=="both"||which=="all");
    const bool doInt  = (which=="all");

    index_t nWritten = 0;

    // ---- geometry mesh (for context) ----
    gsWriteParaview<real_t>(B.mp, outDir + "geometry", 1000, true);

    // ---- vertex functions ----
    if (doVtx)
        for (const asg1mp::VertexFuncs<real_t>& vf : B.vertexFuncs)
            for (index_t j = 0; j < vf.W.cols(); ++j)
            {
                std::string nm = outDir + "vtx" + std::to_string(vf.vertexId)
                               + "_bf" + std::to_string(j);
                writeField(B, vf.gColBase + j, nm, npts);
                ++nWritten;
            }

    // ---- edge functions (a spread of trace and deriv) ----
    if (doEdge)
        for (const asg1mp::EdgeEmbedding<real_t>& em : B.edgeEmb)
        {
            auto spread = [&](index_t base, index_t cnt, const std::string& tag)
            {
                if (cnt <= 0) return;
                const index_t m = std::min(maxEdge, cnt);
                for (index_t s = 0; s < m; ++s)
                {
                    const index_t k = (m==1) ? cnt/2 : (s*(cnt-1))/(m-1);
                    std::string nm = outDir + "edge" + std::to_string(em.edgeId)
                                   + "_" + tag + std::to_string(k);
                    writeField(B, base + k, nm, npts);
                    ++nWritten;
                }
            };
            spread(em.gOffTrace, em.nSmKept, "trace");
            spread(em.gOffDeriv, em.nLDKept, "deriv");
        }

    // ---- a couple of interior functions ----
    if (doInt && B.edgeBlockOff > 0)
    {
        const index_t step = std::max<index_t>(B.edgeBlockOff/4, 1);
        for (index_t idx = 0; idx < B.edgeBlockOff; idx += step)
        {
            writeField(B, idx, outDir + "interior" + std::to_string(idx), npts);
            ++nWritten;
        }
    }

    gsInfo << "Wrote " << nWritten << " basis-function fields + geometry to "
           << outDir << "\n";
    return 0;
}
