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
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementMarker.h>

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

    gsKnotVector<> kv(0, 1.0, 0, 3, 1);
    kv.uniformRefine();

    gsBSplineBasis<real_t> bsb(kv);

    gsTHBSplineBasis<1,real_t,true> thb(bsb);

    // Refine the basis randomly over some iterations
    {
        thb.uniformRefine();
        for (index_t it=0; it!=8; ++it)
        {
            gsHElementMarker<1, real_t> marker(thb);
            // take a random element index between 0 and the number of elements
            const index_t numElements = thb.domain()->numElements();
            index_t elemIndex = std::rand() % numElements;
            std::vector<real_t> errors(numElements, 0);
            errors[elemIndex] = 1.0; // Mark the randomly chosen element for refinement
            marker.setErrors(errors);
            auto marked = marker.markRef();
            std::vector<index_t> boxes = marker.helper().toRefBoxes(marked);
            thb.refineElements(boxes);
        }
    }

    // std::vector<index_t> box(6);
    // for (index_t i=0; i!=4; ++i)
    // {
    //     box[0] = i+1; // level
    //     box[1] = 0; // lower x
    //     box[2] = 2; // upper x
    //     box[3] = i+1; // level
    //     box[4] = 4*math::pow(2,i)-2; // lower x
    //     box[5] = 4*math::pow(2,i);   // upper x
    //     thb.refineElements(box);
    // }

    if (plot)
        gsWriteParaview(thb,"basis",1000,true);

    // Gauss Legendre
    gsOptionList legendreOpts;
    legendreOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::GaussLegendre);
    legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    legendreOpts.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    legendreOpts.addSwitch("overInt","Apply over-integration or not?",false);
    gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(thb.tensorLevel(thb.maxLevel()), legendreOpts);

    // Construct a patch rule for every level
    gsOptionList patchOpts;
    patchOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::PatchRule);
    patchOpts.addReal("quA", "Order of the target space", order  );
    patchOpts.addInt ("quB", "Regularity of the targed space", regularity    );
    patchOpts.addSwitch("overInt","Apply over-integration or not?",overInt);
    gsQuadRule<real_t>::uPtr patchRule = gsQuadrature::getPtr(thb.tensorLevel(thb.maxLevel()), patchOpts);

    std::vector<gsQuadRule<real_t>::uPtr> rules(thb.maxLevel()+1);
    for (index_t l = 0; l!=thb.maxLevel()+1; l++)
    {
        rules[l] = gsQuadrature::getPtr(thb.tensorLevel(l),patchOpts);
    }

    gsMatrix<index_t> acts;
    gsMatrix<real_t> vals;
    gsMatrix<real_t> qNodes;
    gsVector<real_t> qWeights;

    //////////////////////////////////////////////////////////////////////////
    // Base case: integrate every element with Gauss-Legendre quadrature
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> reference(thb.size());
    reference.setZero();
    for (auto & elem : thb.domain()->allElements())
    {
        legendre->mapTo(elem.lowerCorner(),elem.upperCorner(),qNodes,qWeights);
        thb.active_into(qNodes,acts);
        thb.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            reference[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 1: integrate every element with the patch rule of the finest level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case1(thb.size());
    case1.setZero();
    for (auto & elem : thb.domain()->allElements())
    {
        rules.back()->mapTo(elem.lowerCorner(),elem.upperCorner(),qNodes,qWeights);
        thb.active_into(qNodes,acts);
        thb.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            case1[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 2: integrate every element with the patch rule of its own level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case2(thb.size());
    case2.setZero();
    for (auto & elem : thb.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
        rules[helem->getLevel()]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        thb.active_into(qNodes,acts);
        thb.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            case2[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 3: integrate every basis function with the patch rule of its own level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case3(thb.size());
    case3.setZero();
    for (index_t i = 0; i!=thb.size(); i++)
    {
        gsMatrix<real_t> support = thb.support(i);
        // INEFFICIENT: we loop over all elements to find the support of the i-th basis function, but this is just for demonstration purposes
        for (auto & elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[helem->getLevel()]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb.active_into(qNodes,acts);
                thb.eval_into(qNodes,vals);
                vals = vals * qWeights;
                for (index_t r = 0; r!=acts.rows(); r++)
                    if (acts(r,0) == i)
                        case3[i]+= vals(r,0);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Case 4: Integrate every basis function on the level corresponding to 
    //         level of the finest level its support intersects with
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case4(thb.size());
    case4.setZero();
    for (index_t i = 0; i!=thb.size(); i++)
    {
        gsMatrix<real_t> support = thb.support(i);
        index_t level = 0;
        for (auto & elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                level = std::max(level, (index_t)helem->getLevel());
            }
        }

        for (auto & elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[level]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb.active_into(qNodes,acts);
                thb.eval_into(qNodes,vals);
                vals = vals * qWeights;
                for (index_t r = 0; r!=acts.rows(); r++)
                    if (acts(r,0) == i)
                        case4[i]+= vals(r,0);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Case 5: Integrate functions on their representation level
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case5(thb.size());
    case5.setZero();
    for (index_t i = 0; i!=thb.size(); i++)
    {
        // Representation level
        index_t repLevel = thb.repLevel(i);
        gsMatrix<real_t> support = thb.support(i);
        for (auto & elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[repLevel]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb.active_into(qNodes,acts);
                thb.eval_into(qNodes,vals);
                vals = vals * qWeights;
                for (index_t r = 0; r!=acts.rows(); r++)                    if (acts(r,0) == i)
                        case5[i]+= vals(r,0);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Case 6: Integrate only truncated functions in their own level
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case6(thb.size());
    case6.setZero();
    for (index_t i = 0; i!=thb.size(); i++)
    {
        if (thb.isTruncated(i))
        {
            // Representation level
            // index_t repLevel = thb.repLevel(i);
            index_t level = thb.levelOf(i);
            gsMatrix<real_t> support = thb.support(i);
            for (auto & elem : thb.domain()->allElements())
            {
                auto helem = dynamic_cast<gsHDomainIterator<real_t,1> *>(&elem);
                if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                    (support.col(1).array() >= helem->upperCorner().array()).all())
                {
                    rules[level]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                    thb.active_into(qNodes,acts);
                    thb.eval_into(qNodes,vals);
                    vals = vals * qWeights;
                    for (index_t r = 0; r!=acts.rows(); r++)
                        if (acts(r,0) == i)
                            case6[i]+= vals(r,0);
                }
            }
        }
    }   

    ////////////////////////////////////////////////////////////////////////
    // Print results
    ////////////////////////////////////////////////////////////////////////
    gsInfo<<"Case overview:\n";
    gsInfo<<"Reference: Uses Gauss-Legendre quadrature on every element\n";
    gsInfo<<"Case 1: Uses the patch rule of the finest level on every element\n";
    gsInfo<<"Case 2: Uses the patch rule of the own level on every element\n";
    gsInfo<<"Case 3: Uses the patch rule of the own level on every basis function\n";
    gsInfo<<"Case 4: Uses the patch rule of the finest level intersecting with the support of the basis function on every basis function\n";
    gsInfo<<"Case 5: Uses the patch rule of the representation level on every basis function\n";
    gsInfo<<"Case 6: Uses the patch rule of the own level on truncated basis functions and does not integrate non-truncated functions\n";


    gsInfo<<"Reference:\n"<<reference.transpose()<<"\n";
    gsInfo<<"Case 1:\n"<<case1.transpose()<<"\n";
    gsInfo<<"Case 2:\n"<<case2.transpose()<<"\n";
    gsInfo<<"Case 3:\n"<<case3.transpose()<<"\n";
    gsInfo<<"Case 4:\n"<<case4.transpose()<<"\n";
    gsInfo<<"Case 5:\n"<<case5.transpose()<<"\n";
    gsInfo<<"Case 6:\n"<<case6.transpose()<<"\n";

    // Differences
    gsInfo<<"Case 1 - Reference:\n"<<(case1-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 2 - Reference:\n"<<(case2-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 3 - Reference:\n"<<(case3-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 4 - Reference:\n"<<(case4-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 5 - Reference:\n"<<(case5-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 6 - Reference:\n"<<(case6-reference).cwiseAbs().maxCoeff()<<"\n";

    /////////////////////////////////////////////////////////////////////////
    // Construct a basis 
    /////////////////////////////////////////////////////////////////////////
    gsTHBSplineBasis<1,real_t,true> thb_mod(thb);
    gsHElementHelper<1, real_t> helper(thb_mod);

    // Loop over all truncated functions
    for (index_t i = 0; i!=thb_mod.size(); i++)
    {
        // If the function is truncated, we compute its representation level and refine the element corresponding to the support of the function on that level
        if (thb_mod.isTruncated(i))
        {
            index_t repLevel = thb_mod.repLevel(i);
            gsMatrix<real_t> support = thb_mod.support(i);
            // Convert support to element (box coordinates at repLevel)
            gsVector<real_t, 1> low = support.col(0);
            gsVector<real_t, 1> upp = support.col(1);
            
            auto element = helper.toElement(low, upp, repLevel);
            // Now you have element with box coordinates
            // Access them via:
            const auto & lowCorner = element.lowerCorner();  // Vector<index_t, 2>
            const auto & uppCorner = element.upperCorner();  // Vector<index_t, 2>
            gsInfo << "Element for truncated function " << i << " at representation level " << repLevel << ":\n";
            gsInfo << "  Low corner: " << lowCorner.transpose() << "\n";
            gsInfo << "  Upp corner: " << uppCorner.transpose() << "\n";
        }
    }



    return 0;
}