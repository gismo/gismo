/** @file asg1mp_report.cpp

    @brief Collect AS-G1 multi-patch results for one geometry and append
    them as a CSV row.  Used by the experiment driver to compile result
    tables across many geometries.

    Outputs (per invocation):
      * one data row appended to the --csv file (summary metrics)
      * one row per interior interface appended to the --edgecsv file

    Usage:
      ./bin/asg1mp_report -f <file> -r <refs> --name <label> \
            --csv results.csv --edgecsv per_edge.csv
*/

#include <gismo.h>
#include "asg1_mp/asg1_multipatch_mp.h"

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string fn, name("geom"), csv("results.csv"), edgecsv("per_edge.csv");
    index_t refs = 2, minDeg = 3, nCheck = 25;
    bool useTemplate = false;

    gsCmdLine cmd("AS-G1 multi-patch result collector (CSV).");
    cmd.addString("f", "file",   "Geometry XML file", fn);
    cmd.addString("",  "name",   "Geometry label for the CSV", name);
    cmd.addString("",  "csv",    "Summary CSV file (appended)", csv);
    cmd.addString("",  "edgecsv","Per-edge CSV file (appended)", edgecsv);
    cmd.addInt   ("r", "refs",   "Uniform refinements", refs);
    cmd.addInt   ("d", "degree", "Minimum spline degree", minDeg);
    cmd.addInt   ("n", "ncheck", "Samples per interface", nCheck);
    cmd.addSwitch("template", "Gluing data on the bilinear template (default: actual geometry)", useTemplate);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>(fn, mp);
    if (mp.nPatches() == 0) { gsInfo << "Cannot read " << fn << "\n"; return 1; }
    mp.computeTopology();

    asg1mp::AsG1MpOptions<real_t> opts;
    opts.refinements = refs;
    opts.minDegree   = static_cast<short_t>(minDeg);
    opts.gluingSource = useTemplate ? asg1mp::GluingSource::BilinearTemplate
                                    : asg1mp::GluingSource::Direct;
    opts.verbose     = false;

    gsStopwatch clock;
    asg1mp::AsG1Multipatch<real_t> B = asg1mp::buildAsG1Multipatch<real_t>(mp, opts);
    const double buildTime = clock.stop();

    clock.restart();
    asg1mp::G1ReportMP<real_t> R = asg1mp::smoothnessCheckMP<real_t>(B, nCheck);
    const double checkTime = clock.stop();

    // count interior (>= 2 interfaces) vertices
    index_t nXverts = 0;
    for (const asg1mp::VertexInfo<real_t>& vi : B.verts)
        if (asg1mp::vtxInterfaceCount<real_t>(vi, B.edges) >= 2) ++nXverts;

    // --- polynomial reproduction (sparse normal equations) ---
    const gsMultiPatch<real_t>& mp2 = B.mp;
    const size_t P = mp2.nPatches();
    std::vector<index_t> roffv(P+1, 0);
    index_t totRows = 0;
    for (size_t p = 0; p < P; ++p) { roffv[p]=totRows; totRows += mp2.patch(p).basis().size(); }
    roffv[P] = totRows;
    gsSparseMatrix<real_t> Gst(totRows, B.nGlobal);
    {
        gsSparseEntries<real_t> ent;
        for (size_t p = 0; p < P; ++p)
            for (index_t cc = 0; cc < B.G[p].outerSize(); ++cc)
                for (gsSparseMatrix<real_t>::InnerIterator it(B.G[p], cc); it; ++it)
                    ent.add(roffv[p]+it.row(), it.col(), it.value());
        Gst.setFrom(ent); Gst.makeCompressed();
    }
    gsSparseMatrix<real_t> A = Gst.transpose() * Gst;
    for (index_t i = 0; i < A.rows(); ++i) A.coeffRef(i,i) += 1e-12;
    A.makeCompressed();
    gsSparseSolver<real_t>::SimplicialLDLT solver(A);
    auto reproduce = [&](std::function<real_t(real_t,real_t)> f)->real_t
    {
        gsMatrix<real_t> c(totRows,1);
        for (size_t p = 0; p < P; ++p)
        {
            const gsBasis<real_t>& tb = mp2.patch(p).basis();
            const index_t n = tb.size();
            gsMatrix<real_t> anchors = tb.anchors();
            gsMatrix<real_t> phys = mp2.patch(p).eval(anchors);
            gsMatrix<real_t> fv(1,n);
            for (index_t k=0;k<n;++k) fv(0,k)=f(phys(0,k),phys(1,k));
            typename gsGeometry<real_t>::uPtr geo = tb.interpolateAtAnchors(fv);
            c.block(roffv[p],0,n,1) = geo->coefs();
        }
        gsMatrix<real_t> g = solver.solve(Gst.transpose()*c);
        return (Gst*g-c).norm()/std::max<real_t>(c.norm(),1e-30);
    };
    const real_t r1  = reproduce([](real_t  ,real_t  ){return 1.0;});
    const real_t rx  = reproduce([](real_t x,real_t  ){return x;  });
    const real_t ry  = reproduce([](real_t  ,real_t y){return y;  });
    const real_t rxx = reproduce([](real_t x,real_t  ){return x*x;});
    const real_t rxy = reproduce([](real_t x,real_t y){return x*y;});
    const real_t ryy = reproduce([](real_t  ,real_t y){return y*y;});

    // --- append summary CSV row ---
    std::ofstream os(csv.c_str(), std::ios::app);
    os << name << "," << refs << "," << B.mp.nPatches() << "," << B.mp.nInterfaces()
       << "," << nXverts << "," << B.nGlobal << "," << B.edgeBlockOff
       << "," << (B.vertexBlockOff - B.edgeBlockOff) << "," << (B.nGlobal - B.vertexBlockOff)
       << "," << std::scientific << std::setprecision(4)
       << buildTime << "," << checkTime << "," << R.maxValErr << "," << R.maxGradErr
       << "," << r1 << "," << rx << "," << ry << "," << rxx << "," << rxy << "," << ryy
       << "," << (R.pass ? "PASS" : "FAIL") << "\n";
    os.close();

    // --- append per-edge CSV rows ---
    std::ofstream oe(edgecsv.c_str(), std::ios::app);
    index_t ei = 0;
    for (const asg1mp::EdgeInfo<real_t>& ed : B.edges)
    {
        if (ed.isBoundary) continue;
        const real_t ge = (ei < (index_t)R.perEdgeGradErr.size()) ? R.perEdgeGradErr[ei] : real_t(-1);
        oe << name << "," << refs << "," << ed.id
           << ",p" << ed.side1.patch << "s" << ed.side1.side()
           << ",p" << ed.side2.patch << "s" << ed.side2.side()
           << "," << (ed.flipped ? 1 : 0)
           << "," << std::scientific << std::setprecision(4) << ge << "\n";
        ++ei;
    }
    oe.close();

    gsInfo << name << ": nGlobal=" << B.nGlobal << " maxGradErr=" << R.maxGradErr
           << " repro(1,x2)=" << r1 << "," << rxx
           << " build=" << buildTime << "s  " << (R.pass ? "PASS" : "FAIL") << "\n";
    return 0;
}
