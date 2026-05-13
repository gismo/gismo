/** @file quadrature_example.cpp

    @brief Example of different quadratyre rules

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

    gsKnotVector<> kv1(0, 1.0, 0, 3, 1);
    kv1.uniformRefine();
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2,real_t> tbsb2(kv1,kv2);

    gsTHBSplineBasis<2,real_t,true> thb2(tbsb2);

    std::vector<index_t> box(10);
    for (index_t i=0; i!=4; ++i)
    {
        box[0] = i+1; // level
        box[1] = 0; // lower x
        box[2] = 0; // lower y
        box[3] = 2; // upper x
        box[4] = 2; // upper y
        box[5] = i+1; // level
        box[6] = 4*math::pow(2,i)-2; // lower x
        box[7] = 4*math::pow(2,i)-2; // lower y
        box[8] = 4*math::pow(2,i);   // upper x
        box[9] = 4*math::pow(2,i);   // upper y
        thb2.refineElements(box);
    }

    if (plot)
        gsWriteParaview(thb2,"basis",1000,true);

    // Gauss Legendre
    gsOptionList legendreOpts;
    legendreOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::GaussLegendre);
    legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    legendreOpts.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    legendreOpts.addSwitch("overInt","Apply over-integration or not?",false);
    gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(thb2.tensorLevel(thb2.maxLevel()), legendreOpts);

    // Construct a patch rule for every level
    gsOptionList patchOpts;
    patchOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::PatchRule);
    patchOpts.addReal("quA", "Order of the target space", order  );
    patchOpts.addInt ("quB", "Regularity of the targed space", regularity    );
    patchOpts.addSwitch("overInt","Apply over-integration or not?",overInt);
    gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(thb2.tensorLevel(thb2.maxLevel()), patchOpts);

    std::vector<gsQuadRule<real_t>::uPtr> rules(thb2.maxLevel()+1);
    for (index_t l = 0; l!=thb2.maxLevel()+1; l++)
    {
        rules[l] = gsQuadrature::getPtr(thb2.tensorLevel(l),patchOpts);
    }

    gsVector<real_t> integral(thb2.size());
    integral.setZero();

    gsMatrix<index_t> acts;
    gsMatrix<real_t> vals;
    gsMatrix<real_t> qNodes;
    gsVector<real_t> qWeights;
    gsDebugVar(rules.size());
    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        
        gsInfo<<helem->getLevel()<<": "<<helem->lowerCorner().transpose()<<" , "<<helem->upperCorner().transpose()<<"\n";
        rules[helem->getLevel()]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        patchRule->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        thb2.active_into(qNodes,acts);
        thb2.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            integral[acts(r,0)]+= vals(r,0);
    }

    gsInfo<<integral.transpose();

    // // ======================================================================
    // // some properties
    // // ======================================================================
    // if (verbose)
    //     gsInfo << tbsb2 <<"\n";

    // // ======================================================================
    // // Define quadrature rules
    // // ======================================================================

    // // PatchRule
    // gsOptionList patchOpts;
    // patchOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::PatchRule);
    // patchOpts.addReal("quA", "Order of the target space", order  );
    // patchOpts.addInt ("quB", "Regularity of the targed space", regularity    );
    // patchOpts.addSwitch("overInt","Apply over-integration or not?",overInt);
    // gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(tbsb2, patchOpts);

    // gsMatrix<> points;
    // gsVector<> weights;

    // // --------------------------------------------------------------------------------------

    // typename gsBasis<>::domainIter domIt = tbsb2.domain()->beginAll();
    // typename gsBasis<>::domainIter domItEnd = tbsb2.domain()->endAll();

    // gsMatrix<> TensorPatch(tbsb2.dim(),0);
    // index_t start;
    // for (; domIt<domItEnd; ++domIt )
    // {
    //     if (verbose)
    //     {
    //         gsInfo<<"---------------------------------------------------------------------------\n";
    //         gsInfo  <<"Element with corners (lower) "
    //                 <<domIt.lowerCorner().transpose()<<" and (higher) "
    //                 <<domIt.upperCorner().transpose()<<" :\n";
    //     }
    //     //---------------------------------------------------------------------------
    //     //  Patch-rule
    //     patchRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
    //                     points, weights);
    //     if (verbose)
    //     {
    //         gsInfo  <<"* \t PatchRule ("<<( overInt ? "" : "not " )<<"overintegrated)\n"
    //                 <<"- points:\n"<<points<<"\n"
    //                 <<"- weights:\n"<<weights.transpose()<<"\n";
    //     }
    //     start = TensorPatch.cols();
    //     TensorPatch.conservativeResize(gsEigen::NoChange,TensorPatch.cols()+points.cols());
    //     TensorPatch.block(0,start,TensorPatch.rows(),points.cols()) = points;
    //     //---------------------------------------------------------------------------
    // }
    // if (verbose)
    //     gsInfo<<"---------------------------------------------------------------------------\n";

    // if (plot)
    // {
    //     gsWriteParaviewPoints(TensorPatch,"Points_Patch");
    // }
    // else
    // {
    //     gsInfo<<"No plot produced! Re-run with --plot to export points and basis to Paraview\n";
    // }

    // patchRule = gsQuadrature::getPtr(tbsb2, patchOpts,1);


    // boxSide side(4);
    // typename gsBasis<>::domainIter bIt = tbsb2.domain()->beginBdr(side);
    // typename gsBasis<>::domainIter bItEnd = tbsb2.domain()->endBdr(side);
    // // Start iteration over elements
    // TensorPatch.resize(tbsb2.dim(),0);
    // for (; bIt<bItEnd; ++bIt )
    // {
    //     if (verbose)
    //     {
    //         gsInfo<<"---------------------------------------------------------------------------\n";
    //         gsInfo  <<"Element with corners (lower) "
    //                 <<bIt.lowerCorner().transpose()<<" and (higher) "
    //                 <<bIt.upperCorner().transpose()<<" :\n";
    //     }

    //     //  Patch-rule
    //     patchRule->mapTo( bIt.lowerCorner(), bIt.upperCorner(),
    //                     points, weights);
    //     if (verbose)
    //     {
    //         gsInfo  <<"* \t PatchRule ("<<( overInt ? "" : "not " )<<"overintegrated)\n"
    //                 <<"- points:\n"<<points<<"\n"
    //                 <<"- weights:\n"<<weights.transpose()<<"\n";
    //     }
    //     start = TensorPatch.cols();
    //     TensorPatch.conservativeResize(gsEigen::NoChange,TensorPatch.cols()+points.cols());
    //     TensorPatch.block(0,start,TensorPatch.rows(),points.cols()) = points;

    // }

    // if (verbose)
    //     gsInfo<<"---------------------------------------------------------------------------\n";

    // if (plot)
    // {
    //     gsWriteParaviewPoints(TensorPatch,"BoundaryPoints_Patch");
    // }
    // else
    // {
    //     gsInfo<<"No plot produced! Re-run with --plot to export points and basis to Paraview\n";
    // }

    return 0;
}
