/** @file quadrature_example.cpp

    @brief Example of different quadratyre rules

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
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2,real_t> tbsb2(kv1,kv2);

    if (plot)
        gsWriteParaview(tbsb2,"basis",1000,true);

    // ======================================================================
    // some properties
    // ======================================================================
    if (verbose)
        gsInfo << tbsb2 <<"\n";

    // ======================================================================
    // Define quadrature rules
    // ======================================================================

    // Gauss Legendre
    gsOptionList options;
    options.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::rule::GaussLegendre);
    options.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    options.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    options.addSwitch("overInt","Apply over-integration or not?",false);

    gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(tbsb2, options);

    // Mixed Quadrature
    options.setInt   ("quRule" ,gsQuadrature::rule::GaussLobatto);
    options.setSwitch("overInt",true);
    options.setReal("quA", 0  ); // interior
    options.setInt ("quB", 1  );
    options.addReal("quAb", "Number of quadrature points: quA*deg + quB", 1.0  );
    options.addInt ("quBb", "Number of quadrature points: quA*deg + quB", 1    );
    gsQuadRule<real_t>::uPtr mixedLobatto = gsQuadrature::getPtr(tbsb2, options);

    // PatchRule
    options.setInt   ("quRule",gsQuadrature::rule::PatchRule);
    options.setReal  ("quA", order      );
    options.setInt   ("quB", regularity );
    options.setSwitch("overInt",overInt);
    gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(tbsb2, options);

    gsMatrix<> points;
    gsVector<> weights;

    // --------------------------------------------------------------------------------------

    typename gsBasis<real_t>::domainIter domIt = tbsb2.makeDomainIterator();

    gsMatrix<> GaussRule(tbsb2.dim(),0);
    gsMatrix<> MixedRule(tbsb2.dim(),0);
    gsMatrix<> TensorPatch(tbsb2.dim(),0);
    index_t start;
    for (; domIt->good(); domIt->next() )
    {
        if (verbose)
        {
            gsInfo<<"---------------------------------------------------------------------------\n";
            gsInfo  <<"Element with corners (lower) "
                    <<domIt->lowerCorner().transpose()<<" and (higher) "
                    <<domIt->upperCorner().transpose()<<" :\n";
        }

        //---------------------------------------------------------------------------
        // Gauss-Legendre rule (w/o over-integration)
        legendre->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t Gauss-Legendre\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = GaussRule.cols();
        GaussRule.conservativeResize(Eigen::NoChange,GaussRule.cols()+points.cols());
        GaussRule.block(0,start,GaussRule.rows(),points.cols()) = points;

        //---------------------------------------------------------------------------
        // Gauss-Lobatto rule (w/ over-integration)
        mixedLobatto->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t Gauss-Lobatto (overintegrated)\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = MixedRule.cols();
        MixedRule.conservativeResize(Eigen::NoChange,MixedRule.cols()+points.cols());
        MixedRule.block(0,start,MixedRule.rows(),points.cols()) = points;
        //---------------------------------------------------------------------------
        //  Patch-rule
        patchRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        if (verbose)
        {
            gsInfo  <<"* \t PatchRule ("<<( overInt ? "" : "not " )<<"overintegrated)\n"
                    <<"- points:\n"<<points<<"\n"
                    <<"- weights:\n"<<weights.transpose()<<"\n";
        }
        start = TensorPatch.cols();
        TensorPatch.conservativeResize(Eigen::NoChange,TensorPatch.cols()+points.cols());
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

    return 0;
}