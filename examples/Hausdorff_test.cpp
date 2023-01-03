/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsSpectra/gsSpectra.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

#include <gsKLShell/gsThinShellUtils.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot       = false;
    bool mesh       = false;
    index_t numRefine  = 2;
    index_t degree = 3;
    index_t smoothness = 2;

    std::string fn1,fn2;
    fn1 = "pde/2p_square_geom.xml";
    fn2 = "pde/2p_square_bvp.xml";
    std::string out = "ModalResults";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addString( "g", "ori","File containing the geometry",  fn1 );
    cmd.addString( "G", "def","File containing the geometry",  fn2 );
    cmd.addString( "o", "out", "Output directory",  out );
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("mesh", "mesh",mesh);

    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp, mp_def;

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");

    gsFileData<> fd;
    gsInfo<<"Reading geometry from "<<fn1<<"..."<<std::flush;
    gsReadFile<>(fn1, mp);
    if (mp.nInterfaces()==0 && mp.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it..."<<std::flush;
        mp.computeTopology();
    }
    if (mp.geoDim()==2)
        mp.embed(3);

    gsReadFile<>(fn2, mp_def);

    gsInfo<<"Finished\n";
    if (mp_def.nInterfaces()==0 && mp_def.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it..."<<std::flush;
        mp_def.computeTopology();
    }
    if (mp.geoDim()==2)
        mp.embed(3);

    // gsInfo<<"Setting degree and refinement..."<<std::flush;
    // GISMO_ENSURE(degree>=mp.patch(0).degree(0),"Degree must be larger than or equal to the degree of the initial geometry, but degree = "<<degree<<" and the original degree = "<<mp.patch(0).degree(0));
    // mp.degreeIncrease(degree-mp.patch(0).degree(0));

    // // h-refine each basis
    // for (int r =0; r < numRefine; ++r)
    //     mp.uniformRefine(1,degree-smoothness);
    // gsInfo<<"Finished\n";
    // gsInfo<<"Patch 0 has basis: "<<mp.basis(0)<<"\n";

    if (plot)
        gsWriteParaview(mp,out + "/" + "mp",200,true);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    gsInfo<<"Computing Hausdorff distance..."<<std::flush;
    std::vector<real_t> hausdorffs = mp_def.HausdorffDistance(mp,1000,1e-6,true);
    gsInfo<<"Finished.\n";
    gsPiecewiseFunction<real_t> distances(mp.nPatches());
    std::vector<gsConstantFunction<real_t>> funcs(mp.nPatches());
    real_t area;
    for (size_t p = 0; p!=hausdorffs.size(); p++)
    {
        gsExprEvaluator<> ev;
        gsMultiPatch<> mp_tmp(mp_def.patch(p));
        gsMultiBasis<> dbasis(mp_tmp);
        ev.setIntegrationElements(dbasis);
        geometryMap G = ev.getMap(mp_tmp);
        area = ev.integral(meas(G));
        funcs.at(p) = gsConstantFunction<>(hausdorffs.at(p) / std::sqrt(area),2);
        distances.addPiece(funcs.at(p));
        gsDebugVar(hausdorffs.at(p));
    }

    gsInfo<<"Plotting Hausdorff distance..."<<std::flush;
    gsWriteParaview(mp,distances,"distField",100);
    gsInfo<<"Finished.\n";

    return EXIT_SUCCESS;
}
