/** @file quadrature_example.cpp

    @brief Example of different quadrature rules

    Examples from the following paper can be created with this file.
    Adam, C., Hughes, T. J. R., Bouabdallah, S., Zarroug, M., & Maitournam, H. (2015). Selective and reduced numerical integrations for NURBS-based isogeometric analysis. Computer Methods in Applied Mechanics and Engineering, 284, 732–761. https://doi.org/10.1016/j.cma.2014.11.001

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================


    index_t order = 2;
    index_t regularity = 1;
    bool plot = false;
    bool verbose = false;
    bool overInt = false;

    gsCmdLine cmd("Quadrature rules in G+Smo.");
    cmd.addInt("P","deg","order of target space",order);
    cmd.addInt("R","reg","regularity of target space",regularity);
    cmd.addSwitch("plot","Plot with paraview",plot);
    cmd.addSwitch("verbose","Verbose points and weights",verbose);
    cmd.addSwitch("over","overintegrate",overInt);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsKnotVector<> kv1(0, 1.0, 3, 3, 1);

    gsTensorBSplineBasis<1,real_t> tbs(kv1);
    gsTHBSplineBasis<1,real_t> thb(tbs);

    thb.refineElements({1,2,6});
    thb.refineElements({2,6,10});

    if (plot)
    {
        gsWriteParaview(thb,"basis",1000,true);
        gsMesh<> mesh(thb);
        gsWriteParaview(mesh,"mesh");
    }

    // ======================================================================
    // some properties
    // ======================================================================
    if (verbose)
        gsInfo << thb <<"\n";

    // ======================================================================
    // Define quadrature rules
    // ======================================================================

    // Gauss Legendre
    gsOptionList legendreOpts;
    legendreOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::GaussLegendre);
    legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    legendreOpts.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    legendreOpts.addSwitch("overInt","Apply over-integration or not?",false);
    gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(thb, legendreOpts);

    // Mixed Quadrature
    gsOptionList lobattoOpts;
    lobattoOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::GaussLobatto);
    lobattoOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 0.0  );
    lobattoOpts.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    lobattoOpts.addReal("quAb", "Number of quadrature points: quA*deg + quB", 1.0  );
    lobattoOpts.addInt ("quBb", "Number of quadrature points: quA*deg + quB", 1    );
    lobattoOpts.addSwitch("overInt","Apply over-integration or not?",true);
    gsQuadRule<real_t>::uPtr mixedLobatto = gsQuadrature::getPtr(thb, lobattoOpts);

    // PatchRule
    gsOptionList patchOpts;
    patchOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::HPatchRule);
    patchOpts.addReal("quA", "Order of the target space", order  );
    patchOpts.addInt ("quB", "Regularity of the targed space", regularity    );
    patchOpts.addSwitch("overInt","Apply over-integration or not?",overInt);
    gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(thb, patchOpts);

    gsMatrix<> points;
    gsVector<> weights;

    // --------------------------------------------------------------------------------------

    typename gsBasis<>::domainIter domIt = thb.domain()->beginAll();
    typename gsBasis<>::domainIter domItEnd = thb.domain()->endAll();

    gsMatrix<> GaussRule(thb.dim(),0);
    gsMatrix<> MixedRule(thb.dim(),0);
    gsMatrix<> TensorPatch(thb.dim(),0);
    index_t start;
    for (; domIt<domItEnd; ++domIt )
    {
        if (verbose)
        {
            gsInfo<<"---------------------------------------------------------------------------\n";
            gsInfo  <<"Element with corners (lower) "
                    <<domIt.lowerCorner().transpose()<<" and (higher) "
                    <<domIt.upperCorner().transpose()<<" :\n";
        }

        //---------------------------------------------------------------------------
        // Gauss-Legendre rule (w/o over-integration)
        legendre->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t Gauss-Legendre\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = GaussRule.cols();
        GaussRule.conservativeResize(gsEigen::NoChange,GaussRule.cols()+points.cols());
        GaussRule.block(0,start,GaussRule.rows(),points.cols()) = points;

        //---------------------------------------------------------------------------
        // Gauss-Lobatto rule (w/ over-integration)
        mixedLobatto->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t Gauss-Lobatto (overintegrated)\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = MixedRule.cols();
        MixedRule.conservativeResize(gsEigen::NoChange,MixedRule.cols()+points.cols());
        MixedRule.block(0,start,MixedRule.rows(),points.cols()) = points;
        //---------------------------------------------------------------------------
        //  Patch-rule
        patchRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t PatchRule ("<<( overInt ? "" : "not " )<<"overintegrated)\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = TensorPatch.cols();
        TensorPatch.conservativeResize(gsEigen::NoChange,TensorPatch.cols()+points.cols());
        TensorPatch.block(0,start,TensorPatch.rows(),points.cols()) = points;
        //---------------------------------------------------------------------------
    }
    if (verbose)
        gsInfo<<"---------------------------------------------------------------------------\n";

    if (plot)
    {
        gsWriteParaviewPoints(GaussRule,"Points_Original");
        gsWriteParaviewPoints(MixedRule,"Points_Mixed");
        gsWriteParaviewPoints(TensorPatch,"Points_Patch");
    }
    else
    {
        gsInfo<<"No plot produced! Re-run with --plot to export points and basis to Paraview\n";
    }

    if (verbose)
        gsInfo<<"---------------------------------------------------------------------------\n";

    if (plot)
    {
        gsWriteParaviewPoints(GaussRule,"BoundaryPoints_Original");
        gsWriteParaviewPoints(MixedRule,"BoundaryPoints_Mixed");
        gsWriteParaviewPoints(TensorPatch,"BoundaryPoints_Patch");
    }
    else
    {
        gsInfo<<"No plot produced! Re-run with --plot to export points and basis to Paraview\n";
    }

    return 0;
}
