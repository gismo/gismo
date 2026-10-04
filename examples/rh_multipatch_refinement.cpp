/** @file rh_multipatch_refinement.cpp

    @brief Multipatch r-refinement from an analytic density function, using gsAdaptiveMultiPatchBuilder

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot          = false;
    double Intensity   = 9.;
    index_t numRefine  = 4;
    index_t numElevate = 0;
    index_t maxIter    = 50;
    std::string fn("pde/annulus2d_bvp.xml");

    gsCmdLine cmd("Multipatch r-refinement from an analytic density function.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps for the mapping basis", numElevate);
    cmd.addInt("u", "uniformRefine", "Number of uniform h-refinement loops for the mapping basis", numRefine);
    cmd.addReal("I", "intensity", "Intensity of the density function", Intensity);
    cmd.addString("f", "file", "Input XML file with the initial multipatch geometry (id 0)", fn);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the result", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    gsMultiPatch<> mpLeft;// Initial geometry
    fd.getId(0,mpLeft);
    mpLeft.computeTopology(); // the interface sides stored in some files do not match the geometry
    gsInfo<<"The domain is "<< mpLeft.detail() << "\n";

    // Analytical density function (defined on the physical domain)
    gsFunctionExpr<> f("1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2)))",2);
    gsInfo<<"Source function "<< f << "\n";

    //! [Adaptive mapping]
    gsAdaptiveMultiPatchBuilder MAE(mpLeft, numRefine, maxIter, Intensity, 0, numElevate);
    auto density = MAE.buildAnalyticDensity(f);   // one density patch per geometry patch
    MAE.buildMultiPatch(density, 1e-8);           // one adaptive mapping of the unit square per patch
    //! [Adaptive mapping]

    //! [Composition]
    gsMultiBasis<> dbasis(mpLeft, true);
    for (index_t r = 0; r < numRefine; ++r)
        dbasis.uniformRefine();
    gsMultiPatch<> newLeft = MAE.buildColCompMultiPatch(dbasis);
    //! [Composition]

    gsInfo << "Adaptive multipatch geometry:\n" << newLeft.detail() << "\n";

    // Maximal gap of the geometry across the interfaces (0 for a conforming multipatch)
    real_t gap = 0.;
    gsMatrix<> t = gsVector<>::LinSpaced(50, 0., 1.).transpose();
    auto sidePts = [&](const patchSide & ps)
    {
        gsMatrix<> p(2, t.cols());
        const index_t s = ps.index();
        p.row((s <= 2) ? 1 : 0) = t;
        p.row((s <= 2) ? 0 : 1).setConstant((s % 2 == 1) ? 0. : 1.);
        return p;
    };
    for (auto & ifc : newLeft.interfaces())
    {
        const gsMatrix<> a = newLeft.patch(ifc.first().patch).eval(sidePts(ifc.first()));
        const gsMatrix<> b = newLeft.patch(ifc.second().patch).eval(sidePts(ifc.second()));
        gap = std::max(gap, std::min((a - b).cwiseAbs().maxCoeff(), (a - b.rowwise().reverse()).cwiseAbs().maxCoeff()));
    }
    gsInfo << "Max interface gap: " << gap << "\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiBasis<> pbasis(newLeft, true);
        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(pbasis);
        gsExprEvaluator<> ev(A);
        gsExprAssembler<>::geometryMap PPFinal = A.getMap(newLeft);
        auto ff_Psi = A.getCoeff(f, PPFinal);
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection<> collection("ParaviewOutput/solution", ev);
        collection.options().setSwitch("elements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("elementResolution", 16);
        collection.options().setInt("numPoints", 100000);
        collection.newTimeStep(newLeft);
        collection.addField(ff_Psi, "density function");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;
}
