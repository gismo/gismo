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

    gsKnotVector<> kv1(0, 1.0, 0, 3, 1);
    kv1.uniformRefine();
    kv1.uniformRefine();
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2,real_t> tbsb2(kv1,kv2);

    gsTHBSplineBasis<2,real_t,true> thb2(tbsb2);

    // L-shaped hierarchical refinement: 6 levels with orthogonal transition bands.
    // Refines the full-width bottom arm and full-height left arm at each level,
    // creating two transition bands (horizontal and vertical) meeting at a refined corner.
    // This is more demanding than corner-only (more mixed elements across two bands).
    {
        thb2.uniformRefine();   // ensure level 1 exists
        for (index_t lvl = 1; lvl <= 5; ++lvl)
        {
            const index_t N = 4 * (1 << lvl);  // 4*2^lvl elements per direction at level lvl
            std::vector<index_t> boxes;
            // Horizontal arm: bottom 4 rows, full width at this level
            boxes.insert(boxes.end(), {lvl, 0, 0, N, 4});
            // Vertical arm: left 4 columns, full height at this level
            boxes.insert(boxes.end(), {lvl, 0, 0, 4, N});
            thb2.refineElements(boxes);
        }
    }

    // std::vector<index_t> box(10);
    // for (index_t i=0; i!=4; ++i)
    // {
    //     box[0] = i+1; // level
    //     box[1] = 0; // lower x
    //     box[2] = 0; // lower y
    //     box[3] = 2; // upper x
    //     box[4] = 2; // upper y
    //     box[5] = i+1; // level
    //     box[6] = 4*math::pow(2,i)-2; // lower x
    //     box[7] = 4*math::pow(2,i)-2; // lower y
    //     box[8] = 4*math::pow(2,i);   // upper x
    //     box[9] = 4*math::pow(2,i);   // upper y
    //     thb2.refineElements(box);
    // }

    if (plot)
    {
        gsMesh<> mesh(thb2);
        gsWriteParaview(mesh,"mesh");
        // gsWriteParaview(thb2,"basis",1000,true);
    }

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

    // Mass-aware patch rules (integrand B_i*B_j has degree 2*order).
    gsOptionList patchOptsM = patchOpts;
    patchOptsM.setReal("quA", 2*order);
    std::vector<gsQuadRule<real_t>::uPtr> rulesM(thb2.maxLevel()+1);
    for (index_t l = 0; l!=thb2.maxLevel()+1; l++)
    {
        rulesM[l] = gsQuadrature::getPtr(thb2.tensorLevel(l),patchOptsM);
    }

    gsMatrix<index_t> acts;
    gsMatrix<real_t> vals;
    gsMatrix<real_t> qNodes;
    gsVector<real_t> qWeights;

    //////////////////////////////////////////////////////////////////////////
    // Base case: integrate every element with Gauss-Legendre quadrature
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> reference(thb2.size());
    reference.setZero();
    for (auto & elem : thb2.domain()->allElements())
    {
        legendre->mapTo(elem.lowerCorner(),elem.upperCorner(),qNodes,qWeights);
        thb2.active_into(qNodes,acts);
        thb2.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            reference[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 1: integrate every element with the patch rule of the finest level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case1(thb2.size());
    case1.setZero();
    for (auto & elem : thb2.domain()->allElements())
    {
        rules.back()->mapTo(elem.lowerCorner(),elem.upperCorner(),qNodes,qWeights);
        thb2.active_into(qNodes,acts);
        thb2.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            case1[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 2: integrate every element with the patch rule of its own level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case2(thb2.size());
    case2.setZero();
    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        rules[helem->getLevel()]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        thb2.active_into(qNodes,acts);
        thb2.eval_into(qNodes,vals);
        vals = vals * qWeights;
        for (index_t r = 0; r!=acts.rows(); r++)
            case2[acts(r,0)]+= vals(r,0);
    }

    //////////////////////////////////////////////////////////////////////////
    // Case 3: integrate every basis function with the patch rule of its own level
    //////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case3(thb2.size());
    case3.setZero();
    for (index_t i = 0; i!=thb2.size(); i++)
    {
        gsMatrix<real_t> support = thb2.support(i);
        // INEFFICIENT: we loop over all elements to find the support of the i-th basis function, but this is just for demonstration purposes
        for (auto & elem : thb2.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[helem->getLevel()]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb2.active_into(qNodes,acts);
                thb2.eval_into(qNodes,vals);
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
    gsVector<real_t> case4(thb2.size());
    case4.setZero();
    for (index_t i = 0; i!=thb2.size(); i++)
    {
        gsMatrix<real_t> support = thb2.support(i);
        index_t level = 0;
        for (auto & elem : thb2.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                level = std::max(level, (index_t)helem->getLevel());
            }
        }

        for (auto & elem : thb2.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[level]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb2.active_into(qNodes,acts);
                thb2.eval_into(qNodes,vals);
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
    gsVector<real_t> case5(thb2.size());
    case5.setZero();
    for (index_t i = 0; i!=thb2.size(); i++)
    {
        // Representation level
        index_t repLevel = thb2.repLevel(i);
        gsMatrix<real_t> support = thb2.support(i);
        for (auto & elem : thb2.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                (support.col(1).array() >= helem->upperCorner().array()).all())
            {
                rules[repLevel]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                thb2.active_into(qNodes,acts);
                thb2.eval_into(qNodes,vals);
                vals = vals * qWeights;
                for (index_t r = 0; r!=acts.rows(); r++)                    if (acts(r,0) == i)
                        case5[i]+= vals(r,0);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Case 6: Integrate only truncated functions in their own level
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case6(thb2.size());
    case6.setZero();
    for (index_t i = 0; i!=thb2.size(); i++)
    {
        if (thb2.isTruncated(i))
        {
            // Representation level
            // index_t repLevel = thb2.repLevel(i);
            index_t level = thb2.levelOf(i);
            gsMatrix<real_t> support = thb2.support(i);
            for (auto & elem : thb2.domain()->allElements())
            {
                auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
                if ((support.col(0).array() <= helem->lowerCorner().array()).all() &&
                    (support.col(1).array() >= helem->upperCorner().array()).all())
                {
                    rules[level]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
                    thb2.active_into(qNodes,acts);
                    thb2.eval_into(qNodes,vals);
                    vals = vals * qWeights;
                    for (index_t r = 0; r!=acts.rows(); r++)
                        if (acts(r,0) == i)
                            case6[i]+= vals(r,0);
                }
            }
        }
    }   

    ////////////////////////////////////////////////////////////////////////
    // Case 8: per-pair patch-rule dispatch.
    //
    //   For each integrand pick a single patch rule by representation level:
    //     per-basis I_i  -> rules[ repLevel(i) ]
    //     per-pair  M_ij -> rules[ max(repLevel(i), repLevel(j)) ]
    //
    //   No Gauss in this variant; no element-level classification. To swap
    //   in a patch-safe classifier (Loibl analogue: Gauss fallback when the
    //   pair isn't level-pure / both untruncated), only the chooseLevel
    //   body needs to change.
    ////////////////////////////////////////////////////////////////////////

    // Reference mass matrix M_ref via the existing all-Gauss machinery.
    // Also accumulate Gauss points / weights for Paraview export.
    std::vector<real_t> allGX, allGY, allGW;
    gsMatrix<real_t> Mref = gsMatrix<real_t>::Zero(thb2.size(), thb2.size());
    for (auto & elem : thb2.domain()->allElements())
    {
        legendre->mapTo(elem.lowerCorner(),elem.upperCorner(),qNodes,qWeights);
        thb2.active_into(qNodes, acts);
        thb2.eval_into (qNodes, vals);
        for (index_t k = 0; k != qWeights.size(); ++k)
        {
            allGX.push_back(qNodes(0,k));
            allGY.push_back(qNodes(1,k));
            allGW.push_back(qWeights(k));
            for (index_t r = 0; r != acts.rows(); ++r)
                for (index_t c = 0; c != acts.rows(); ++c)
                    Mref(acts(r,0),acts(c,0)) += qWeights(k) * vals(r,k) * vals(c,k);
        }
    }

    gsVector<real_t> case8(thb2.size());                       case8.setZero();
    gsMatrix<real_t> M8(thb2.size(),thb2.size());              M8.setZero();
    const index_t Lmax = thb2.maxLevel();
    std::vector<index_t> ptsPerLevel(Lmax+1, 0);
    index_t totalGaussPts8 = 0;
    // Accumulate patch points / weights per level for Paraview export.
    std::vector<std::vector<real_t>> allPX(Lmax+1), allPY(Lmax+1), allPW(Lmax+1);

    // Local storage reused across elements.
    gsMatrix<real_t> gNodes8; gsVector<real_t> gWeights8;
    gsMatrix<real_t> gVals8;  gsMatrix<index_t> gActs8;
    gsMatrix<real_t> ctr2d(2,1);
    gsMatrix<index_t> actCtr;
    std::vector<gsMatrix<real_t>> pNodes  (Lmax+1);
    std::vector<gsVector<real_t>> pWeights(Lmax+1);
    std::vector<gsMatrix<real_t>> pVals   (Lmax+1);
    std::vector<gsMatrix<index_t>> pActs  (Lmax+1);

    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        const auto & low = helem->lowerCorner();
        const auto & upp = helem->upperCorner();

        // Classify via one-point active query (element center).
        ctr2d.col(0) = 0.5 * (low + upp);
        thb2.active_into(ctr2d, actCtr);

        bool isClean = (actCtr.rows() > 0);
        index_t L_elem = isClean ? thb2.repLevel(actCtr(0,0)) : -1;
        for (index_t r = 0; r != actCtr.rows() && isClean; ++r)
        {
            const index_t i = actCtr(r,0);
            if (thb2.isTruncated(i) || thb2.repLevel(i) != L_elem)
                isClean = false;
        }

        if (isClean && L_elem >= 0 && rulesM[L_elem])
        {
            // Clean element: patch rule ONLY — no Gauss.
            rulesM[L_elem]->mapTo(low, upp, gNodes8, gWeights8);
            thb2.active_into(gNodes8, gActs8);
            thb2.eval_into  (gNodes8, gVals8);
            ptsPerLevel[L_elem] += gWeights8.size();
            for (index_t k = 0; k != gWeights8.size(); ++k)
            {
                allPX[L_elem].push_back(gNodes8(0,k));
                allPY[L_elem].push_back(gNodes8(1,k));
                allPW[L_elem].push_back(gWeights8(k));
            }

            // All pairs on a clean element use the patch rule (uniform dispatch).
            for (index_t r = 0; r != gActs8.rows(); ++r)
            {
                const index_t i = gActs8(r,0);
                for (index_t k = 0; k != gWeights8.size(); ++k)
                    case8[i] += gWeights8(k) * gVals8(r,k);
            }
            for (index_t r = 0; r != gActs8.rows(); ++r)
            for (index_t c = r; c != gActs8.rows(); ++c)
            {
                real_t s = 0.0;
                for (index_t k = 0; k != gWeights8.size(); ++k)
                    s += gWeights8(k) * gVals8(r,k) * gVals8(c,k);
                const index_t i = gActs8(r,0), jj = gActs8(c,0);
                M8(i,jj) += s;
                if (jj != i) M8(jj,i) += s;
            }
        }
        else
        {
            // Mixed element: Gauss for all pairs + selective patch for
            // non-truncated same-level pairs (entry-classifier on this element only).
            legendre->mapTo(low, upp, gNodes8, gWeights8);
            thb2.active_into(gNodes8, gActs8);
            thb2.eval_into  (gNodes8, gVals8);
            totalGaussPts8 += gWeights8.size();

            // Determine levels needing patch on this mixed element.
            std::set<index_t> patchLevels;
            for (index_t r = 0; r != gActs8.rows(); ++r)
            {
                const index_t i = gActs8(r,0);
                if (!thb2.isTruncated(i))
                    patchLevels.insert(thb2.repLevel(i));
            }
            for (index_t L = 0; L <= Lmax; ++L) pNodes[L].resize(0,0);
            for (index_t L : patchLevels)
            {
                if (!rulesM[L]) continue;
                rulesM[L]->mapTo(low, upp, pNodes[L], pWeights[L]);
                if (pNodes[L].cols() > 0)
                {
                    thb2.active_into(pNodes[L], pActs[L]);
                    thb2.eval_into  (pNodes[L], pVals[L]);
                    ptsPerLevel[L] += pWeights[L].size();
                    for (index_t k = 0; k != pWeights[L].size(); ++k)
                    {
                        allPX[L].push_back(pNodes[L](0,k));
                        allPY[L].push_back(pNodes[L](1,k));
                        allPW[L].push_back(pWeights[L](k));
                    }
                }
            }

            // Per-basis: non-truncated → patch; truncated → Gauss.
            for (index_t r = 0; r != gActs8.rows(); ++r)
            {
                const index_t i = gActs8(r,0);
                const index_t L = thb2.repLevel(i);
                if (!thb2.isTruncated(i) && pNodes[L].cols() > 0)
                {
                    for (index_t rr = 0; rr != pActs[L].rows(); ++rr)
                        if (pActs[L](rr,0) == i)
                        {
                            for (index_t k = 0; k != pWeights[L].size(); ++k)
                                case8[i] += pWeights[L](k) * pVals[L](rr,k);
                            break;
                        }
                }
                else
                {
                    for (index_t k = 0; k != gWeights8.size(); ++k)
                        case8[i] += gWeights8(k) * gVals8(r,k);
                }
            }

            // Per-pair mass: non-trunc same-level → patch; else → Gauss.
            for (index_t r = 0; r != gActs8.rows(); ++r)
            for (index_t c = r; c != gActs8.rows(); ++c)
            {
                const index_t i  = gActs8(r,0);
                const index_t jj = gActs8(c,0);
                const index_t Li = thb2.repLevel(i);
                const index_t Lj = thb2.repLevel(jj);
                real_t s = 0.0;
                if (!thb2.isTruncated(i) && !thb2.isTruncated(jj) && Li == Lj
                    && pNodes[Li].cols() > 0)
                {
                    index_t ri = -1, rj = -1;
                    for (index_t rr = 0; rr != pActs[Li].rows(); ++rr)
                    {
                        if (pActs[Li](rr,0) == i)  ri = rr;
                        if (pActs[Li](rr,0) == jj) rj = rr;
                    }
                    if (ri >= 0 && rj >= 0)
                        for (index_t k = 0; k != pWeights[Li].size(); ++k)
                            s += pWeights[Li](k) * pVals[Li](ri,k) * pVals[Li](rj,k);
                }
                else
                {
                    for (index_t k = 0; k != gWeights8.size(); ++k)
                        s += gWeights8(k) * gVals8(r,k) * gVals8(c,k);
                }
                const index_t i2 = gActs8(r,0), jj2 = gActs8(c,0);
                M8(i2,jj2) += s;
                if (jj2 != i2) M8(jj2,i2) += s;
            }   // end for c
        }   // end else (mixed element)
    }   // end for elem

    ////////////////////////////////////////////////////////////////////////
    // Case 9: own-level patch rule on ALL elements — no Gauss, no classification.
    // Uses rulesM[L_elem] for every leaf element regardless of clean/mixed.
    // Point count matches the "Total (naive own-level)" in the Case 8 report.
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case9(thb2.size());
    case9.setZero();
    gsMatrix<real_t> M9(thb2.size(), thb2.size());
    M9.setZero();

    gsMatrix<real_t> gNodes9;
    gsVector<real_t> gWeights9;
    gsMatrix<index_t> gActs9;
    gsMatrix<real_t>  gVals9;

    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        const index_t L_e = helem->getLevel();
        if (!rulesM[L_e]) continue;
        rulesM[L_e]->mapTo(helem->lowerCorner(), helem->upperCorner(),
                           gNodes9, gWeights9);
        if (gWeights9.size() == 0) continue;
        thb2.active_into(gNodes9, gActs9);
        thb2.eval_into  (gNodes9, gVals9);
        for (index_t r = 0; r != gActs9.rows(); ++r)
        {
            const index_t i = gActs9(r,0);
            for (index_t k = 0; k != gWeights9.size(); ++k)
                case9[i] += gWeights9(k) * gVals9(r,k);
        }
        for (index_t r = 0; r != gActs9.rows(); ++r)
        for (index_t c = r; c != gActs9.rows(); ++c)
        {
            real_t s = 0.0;
            for (index_t k = 0; k != gWeights9.size(); ++k)
                s += gWeights9(k) * gVals9(r,k) * gVals9(c,k);
            const index_t i = gActs9(r,0), jj = gActs9(c,0);
            M9(i,jj) += s;
            if (jj != i) M9(jj,i) += s;
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Case 10: rep-level dispatch — use rulesM[max(repLevel(i), repLevel(j))]
    // Groups pairs by their max repLevel and evaluates the appropriate rule.
    // Hypothesis: T-T pairs need their own representation-level rule, not element level.
    ////////////////////////////////////////////////////////////////////////
    gsVector<real_t> case10(thb2.size());
    case10.setZero();
    gsMatrix<real_t> M10(thb2.size(), thb2.size());
    M10.setZero();

    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        const auto & low = helem->lowerCorner();
        const auto & upp = helem->upperCorner();

        // Get active functions via center point
        gsMatrix<real_t> ctr(2,1);
        ctr.col(0) = 0.5 * (low + upp);
        gsMatrix<index_t> actCtr;
        thb2.active_into(ctr, actCtr);

        // Build set of distinct max(repLevel) values
        std::set<index_t> Lp_set;
        for (index_t r = 0; r != actCtr.rows(); ++r)
        for (index_t c = r; c != actCtr.rows(); ++c)
        {
            const index_t Li = thb2.repLevel(actCtr(r,0));
            const index_t Lj = thb2.repLevel(actCtr(c,0));
            Lp_set.insert(std::max(Li, Lj));
        }

        // Evaluate each rule once and accumulate matching pairs
        for (index_t Lp : Lp_set)
        {
            if (!rulesM[Lp]) continue;
            gsMatrix<real_t> nodes;
            gsVector<real_t> weights;
            rulesM[Lp]->mapTo(low, upp, nodes, weights);
            if (weights.size() == 0) continue;

            gsMatrix<index_t> acts;
            gsMatrix<real_t> vals;
            thb2.active_into(nodes, acts);
            thb2.eval_into(nodes, vals);

            // Per-basis: accumulate only for functions whose repLevel == Lp
            for (index_t r = 0; r != acts.rows(); ++r)
            {
                const index_t i = acts(r,0);
                if (thb2.repLevel(i) != Lp) continue;
                for (index_t k = 0; k != weights.size(); ++k)
                    case10[i] += weights(k) * vals(r,k);
            }

            // Per-pair: accumulate pairs whose max(repLevel) == Lp
            for (index_t r = 0; r != acts.rows(); ++r)
            for (index_t c = r; c != acts.rows(); ++c)
            {
                const index_t i = acts(r,0);
                const index_t j = acts(c,0);
                const index_t Li = thb2.repLevel(i);
                const index_t Lj = thb2.repLevel(j);
                if (std::max(Li, Lj) != Lp) continue;

                real_t s = 0.0;
                for (index_t k = 0; k != weights.size(); ++k)
                    s += weights(k) * vals(r,k) * vals(c,k);
                M10(i,j) += s;
                if (j != i) M10(j,i) += s;
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
    gsInfo<<"Case 8: per-pair patch-rule dispatch (I_i: rules[repLevel(i)]; M_ij: rules[max(repLevel(i),repLevel(j))])\n";
    gsInfo<<"Case 9: own-level patch rule on ALL elements, no Gauss, no classification\n";
    gsInfo<<"Case 10: rep-level dispatch — use rules[max(repLevel(i),repLevel(j))] for all pairs\n";
    gsInfo<<"\n";

    // Differences
    gsInfo<<"Case 1 - Reference:\n"<<(case1-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 2 - Reference:\n"<<(case2-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 3 - Reference:\n"<<(case3-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 4 - Reference:\n"<<(case4-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 5 - Reference:\n"<<(case5-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 6 - Reference:\n"<<(case6-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 8 - Reference:\n"<<(case8-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 9 - Reference:\n"<<(case9-reference).cwiseAbs().maxCoeff()<<"\n";
    gsInfo<<"Case 10 - Reference:\n"<<(case10-reference).cwiseAbs().maxCoeff()<<"\n";
    // Basis composition
    {
        std::map<index_t,index_t> lvlCount;
        for (index_t i = 0; i != thb2.size(); ++i)
            ++lvlCount[thb2.levelOf(i)];
        gsInfo<<"Basis composition (levelOf -> count): ";
        for (auto & p : lvlCount) gsInfo<<p.first<<":"<<p.second<<" ";
        gsInfo<<"  maxLevel="<<thb2.maxLevel()<<"\n";

        std::map<index_t,index_t> leafLvlCount;
        for (auto & elem : thb2.domain()->allElements())
        {
            auto h = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            ++leafLvlCount[h->getLevel()];
        }
        gsInfo<<"Leaf level distribution: ";
        for (auto & p : leafLvlCount) gsInfo<<p.first<<":"<<p.second<<" ";
        gsInfo<<"\n";
    }

    // Compute all-Gauss baseline for point counts
    std::map<index_t,index_t> gaussPtsPerLeafLvl, patchPtsPerLeafLvl;
    for (auto & elem : thb2.domain()->allElements())
    {
        auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
        const index_t lv = helem->getLevel();
        legendre->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        gaussPtsPerLeafLvl[lv] += qWeights.size();
        rulesM[lv]->mapTo(helem->lowerCorner(),helem->upperCorner(),qNodes,qWeights);
        patchPtsPerLeafLvl[lv] += qWeights.size();
    }
    index_t totalG = 0, totalP = 0;
    for (auto & kv : gaussPtsPerLeafLvl)
    {
        totalG += kv.second;
        totalP += patchPtsPerLeafLvl.count(kv.first) ? patchPtsPerLeafLvl.at(kv.first) : 0;
    }

    // Case 8 detailed report.
    {
        const real_t Ierr = (case8 - reference).cwiseAbs().maxCoeff();
        const real_t Merr = (M8   - Mref     ).cwiseAbs().maxCoeff();
        gsInfo<<"Case 8 per-basis ||I - Iref||_inf = "<<Ierr<<"\n";
        gsInfo<<"Case 8 mass     ||M - Mref||_inf = "<<Merr<<"\n";

        real_t IerrU = 0, IerrT = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        {
            const real_t e = std::abs(case8[i] - reference[i]);
            if (thb2.isTruncated(i)) IerrT = std::max(IerrT, e);
            else                     IerrU = std::max(IerrU, e);
        }
        gsInfo<<"  per-basis max err  untruncated="<<IerrU<<"  truncated="<<IerrT<<"\n";

        real_t MerrUU = 0, MerrUT = 0, MerrTT = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        for (index_t j = i; j != thb2.size(); ++j)
        {
            const real_t e = std::abs(M8(i,j) - Mref(i,j));
            const bool ti = thb2.isTruncated(i), tj = thb2.isTruncated(j);
            if (!ti && !tj)      MerrUU = std::max(MerrUU, e);
            else if ( ti &&  tj) MerrTT = std::max(MerrTT, e);
            else                 MerrUT = std::max(MerrUT, e);
        }
        gsInfo<<"  mass max err  U-U="<<MerrUU
              <<"  U-T="<<MerrUT
              <<"  T-T="<<MerrTT<<"\n";

        // Print top-5 absolute mass errors with classification.
        std::vector<std::tuple<real_t,index_t,index_t>> tops;
        for (index_t i = 0; i != thb2.size(); ++i)
            for (index_t j = i; j != thb2.size(); ++j)
                tops.emplace_back(std::abs(M8(i,j)-Mref(i,j)), i, j);
        std::partial_sort(tops.begin(), tops.begin()+std::min<size_t>(8,tops.size()),
                          tops.end(),
                          [](const auto&a, const auto&b){return std::get<0>(a)>std::get<0>(b);});
        gsInfo<<"  top mass errors (|err|, i, j, repL_i, repL_j, truncI, truncJ):\n";
        for (size_t k = 0; k != std::min<size_t>(8,tops.size()); ++k)
        {
            const auto & t = tops[k];
            const index_t i = std::get<1>(t), j = std::get<2>(t);
            gsInfo<<"    "<<std::get<0>(t)<<"  ("<<i<<","<<j<<")"
                  <<" repL=("<<thb2.repLevel(i)<<","<<thb2.repLevel(j)<<")"
                  <<" trunc=("<<thb2.isTruncated(i)<<","<<thb2.isTruncated(j)<<")"
                  <<" M="<<M8(i,j)<<" Mref="<<Mref(i,j)<<"\n";
        }

        // Quadrature point count comparison (using pre-computed values).
        gsInfo<<"\n  Quadrature point comparison (per leaf level):\n";
        gsInfo<<"  leafLvl | Gauss pts | Patch pts (own level)\n";
        for (auto & kv : gaussPtsPerLeafLvl)
        {
            const index_t lv = kv.first;
            const index_t gp = kv.second;
            const index_t pp = patchPtsPerLeafLvl.count(lv) ? patchPtsPerLeafLvl.at(lv) : 0;
            gsInfo<<"  L"<<lv<<"      | "<<gp<<" | "<<pp<<"\n";
        }
        gsInfo<<"  Total (naive own-level)  | "<<totalG<<" | "<<totalP<<"\n";

        const index_t totalPatch8 = std::accumulate(ptsPerLevel.begin(),ptsPerLevel.end(),0);
        gsInfo<<"\n  Case 8 entry-classifier actual point counts:\n";
        gsInfo<<"  Gauss pts (fallback, all leaves) = "<<totalGaussPts8<<"\n";
        gsInfo<<"  Patch pts per level (only where needed):";
        for (index_t L = 0; L <= Lmax; ++L)
            gsInfo<<" L"<<L<<":"<<ptsPerLevel[L];
        gsInfo<<"\n";
        gsInfo<<"  Patch pts total = "<<totalPatch8<<"\n";
        gsInfo<<"  Grand total (Gauss+Patch, with overlap) = "
              <<totalGaussPts8 + totalPatch8<<"\n";
        gsInfo<<"  vs all-Gauss = "<<totalG
              <<"  (ratio "<<(real_t)(totalGaussPts8+totalPatch8)/totalG<<"x)\n";
        gsInfo<<"  Case 9 own-level patch total pts = "<<totalP
              <<"  (ratio "<<(real_t)totalP/totalG<<"x vs all-Gauss)\n";
    }

    // Case 9 detailed report.
    {
        const real_t Ierr9 = (case9 - reference).cwiseAbs().maxCoeff();
        const real_t Merr9 = (M9   - Mref     ).cwiseAbs().maxCoeff();
        gsInfo<<"Case 9 per-basis ||I - Iref||_inf = "<<Ierr9<<"\n";
        gsInfo<<"Case 9 mass     ||M - Mref||_inf = "<<Merr9<<"\n";

        real_t IerrU9 = 0, IerrT9 = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        {
            const real_t e = std::abs(case9[i] - reference[i]);
            if (thb2.isTruncated(i)) IerrT9 = std::max(IerrT9, e);
            else                     IerrU9 = std::max(IerrU9, e);
        }
        gsInfo<<"  per-basis max err  untruncated="<<IerrU9<<"  truncated="<<IerrT9<<"\n";

        real_t MerrUU9 = 0, MerrUT9 = 0, MerrTT9 = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        for (index_t j = i; j != thb2.size(); ++j)
        {
            const real_t e = std::abs(M9(i,j) - Mref(i,j));
            const bool ti = thb2.isTruncated(i), tj = thb2.isTruncated(j);
            if (!ti && !tj)      MerrUU9 = std::max(MerrUU9, e);
            else if ( ti &&  tj) MerrTT9 = std::max(MerrTT9, e);
            else                 MerrUT9 = std::max(MerrUT9, e);
        }
        gsInfo<<"  mass max err  U-U="<<MerrUU9
              <<"  U-T="<<MerrUT9
              <<"  T-T="<<MerrTT9<<"\n";
    }

    // Case 10 detailed report (rep-level dispatch).
    {
        const real_t Ierr10 = (case10 - reference).cwiseAbs().maxCoeff();
        const real_t Merr10 = (M10   - Mref     ).cwiseAbs().maxCoeff();
        gsInfo<<"Case 10 per-basis ||I - Iref||_inf = "<<Ierr10<<"\n";
        gsInfo<<"Case 10 mass     ||M - Mref||_inf = "<<Merr10<<"\n";

        real_t IerrU10 = 0, IerrT10 = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        {
            const real_t e = std::abs(case10[i] - reference[i]);
            if (thb2.isTruncated(i)) IerrT10 = std::max(IerrT10, e);
            else                     IerrU10 = std::max(IerrU10, e);
        }
        gsInfo<<"  per-basis max err  untruncated="<<IerrU10<<"  truncated="<<IerrT10<<"\n";

        real_t MerrUU10 = 0, MerrUT10 = 0, MerrTT10 = 0;
        for (index_t i = 0; i != thb2.size(); ++i)
        for (index_t j = i; j != thb2.size(); ++j)
        {
            const real_t e = std::abs(M10(i,j) - Mref(i,j));
            const bool ti = thb2.isTruncated(i), tj = thb2.isTruncated(j);
            if (!ti && !tj)      MerrUU10 = std::max(MerrUU10, e);
            else if ( ti &&  tj) MerrTT10 = std::max(MerrTT10, e);
            else                 MerrUT10 = std::max(MerrUT10, e);
        }
        gsInfo<<"  mass max err  U-U="<<MerrUU10
              <<"  U-T="<<MerrUT10
              <<"  T-T="<<MerrTT10<<"\n";

        // Case 10 point count using rep-level dispatch
        std::map<index_t,index_t> patchPtsPerLevelC10;
        for (auto & elem : thb2.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            const auto & low = helem->lowerCorner();
            const auto & upp = helem->upperCorner();

            gsMatrix<real_t> ctr(2,1); ctr.col(0) = 0.5 * (low + upp);
            gsMatrix<index_t> actCtr; thb2.active_into(ctr, actCtr);

            std::set<index_t> Lp_set;
            for (index_t r = 0; r != actCtr.rows(); ++r)
            for (index_t c = r; c != actCtr.rows(); ++c)
                Lp_set.insert(std::max(thb2.repLevel(actCtr(r,0)), thb2.repLevel(actCtr(c,0))));

            for (index_t Lp : Lp_set)
            {
                if (!rulesM[Lp]) continue;
                gsMatrix<real_t> nodes; gsVector<real_t> weights;
                rulesM[Lp]->mapTo(low, upp, nodes, weights);
                patchPtsPerLevelC10[Lp] += weights.size();
            }
        }
        index_t totalPatch10 = 0;
        for (auto & kv : patchPtsPerLevelC10) totalPatch10 += kv.second;
        gsInfo<<"  Case 10 total patch pts = "<<totalPatch10
              <<"  (ratio "<<(real_t)totalPatch10/totalG<<"x vs all-Gauss)\n";
    }

    ////////////////////////////////////////////////////////////////////////
    // Assembly timing tests: 5-level L-shaped mesh, Gauss vs Case 10
    ////////////////////////////////////////////////////////////////////////
    gsInfo<<"\n========== ASSEMBLY TIMING TESTS ==========\n";

    // Create 5-level mesh (one level fewer than the 6-level study mesh)
    gsKnotVector<> kv1_asm(0, 1.0, 0, 3, 1);
    kv1_asm.uniformRefine();
    kv1_asm.uniformRefine();
    gsKnotVector<> kv2_asm = kv1_asm;
    gsTensorBSplineBasis<2,real_t> tbsb2_asm(kv1_asm, kv2_asm);
    gsTHBSplineBasis<2,real_t,true> thb2_asm(tbsb2_asm);
    {
        thb2_asm.uniformRefine();
        for (index_t lvl = 1; lvl <= 4; ++lvl)  // 5 levels total (0-4)
        {
            const index_t N = 4 * (1 << lvl);
            std::vector<index_t> boxes;
            boxes.insert(boxes.end(), {lvl, 0, 0, N, 4});
            boxes.insert(boxes.end(), {lvl, 0, 0, 4, N});
            thb2_asm.refineElements(boxes);
        }
    }

    // Create assembly rules for 5-level mesh
    std::vector<gsQuadRule<real_t>::uPtr> rulesM_asm(thb2_asm.maxLevel()+1);
    gsOptionList patchOptsAsmM;
    patchOptsAsmM.addInt("quRule", "Quadrature rule", gsQuadrature::PatchRule);
    patchOptsAsmM.addReal("quA", "Order of target space", 2*order);
    patchOptsAsmM.addInt("quB", "Regularity of target space", regularity);
    patchOptsAsmM.addSwitch("overInt", "Over-integration", false);
    for (index_t l = 0; l != thb2_asm.maxLevel()+1; ++l)
        rulesM_asm[l] = gsQuadrature::getPtr(thb2_asm.tensorLevel(l), patchOptsAsmM);

    // Stiffness rules. The integrand ∇B_i·∇B_j is degree 2(p-1), continuity
    // C^{r-1}, so it needs reduced regularity (quB=regularity-1) to capture the
    // C^0 kink that the C^r mass rule cannot represent. The degree, however,
    // must stay at 2*order: under rep-level dispatch the degree-2(p-1) rule is
    // NOT exact (the region-telescoping argument that makes dispatch exact
    // requires the degree-2p product moments). Target space S(2p, C^{r-1}).
    std::vector<gsQuadRule<real_t>::uPtr> rulesK_asm(thb2_asm.maxLevel()+1);
    gsOptionList patchOptsAsmK = patchOptsAsmM;
    patchOptsAsmK.setInt("quB", std::max(0, regularity - 1));
    for (index_t l = 0; l != thb2_asm.maxLevel()+1; ++l)
        rulesK_asm[l] = gsQuadrature::getPtr(thb2_asm.tensorLevel(l), patchOptsAsmK);

    gsOptionList legendreOptsAsm;
    legendreOptsAsm.addInt("quRule", "Quadrature rule", gsQuadrature::GaussLegendre);
    legendreOptsAsm.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
    legendreOptsAsm.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
    legendreOptsAsm.addSwitch("overInt", "Over-integration", false);
    gsQuadRule<real_t>::uPtr legendre_asm = gsQuadrature::getPtr(thb2_asm.tensorLevel(thb2_asm.maxLevel()), legendreOptsAsm);

    index_t n_elem_asm = 0;
    for (auto & elem : thb2_asm.domain()->allElements()) { (void)elem; ++n_elem_asm; }
    gsInfo<<"Assembly basis: "<<thb2_asm.size()<<" functions, "<<n_elem_asm<<" elements\n";

    // 3a: Gauss assembly with timing (baseline)
    gsMatrix<real_t> M_gauss, K_gauss;
    double time_gauss = 0.0;
    {
        gsStopwatch sw;
        M_gauss = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());
        K_gauss = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());

        for (auto & elem : thb2_asm.domain()->allElements())
        {
            gsMatrix<real_t> nodes; gsVector<real_t> weights;
            legendre_asm->mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
            gsMatrix<index_t> acts; gsMatrix<real_t> vals; gsMatrix<real_t> derivs;
            thb2_asm.active_into(nodes, acts);
            thb2_asm.eval_into(nodes, vals);
            thb2_asm.deriv_into(nodes, derivs);
            for (index_t k = 0; k != weights.size(); ++k)
            {
                for (index_t r = 0; r != acts.rows(); ++r)
                for (index_t c = r; c != acts.rows(); ++c)
                {
                    const index_t i = acts(r,0), j = acts(c,0);
                    real_t sm = weights(k) * vals(r,k) * vals(c,k);
                    M_gauss(i,j) += sm;
                    if (j != i) M_gauss(j,i) += sm;
                    // deriv_into layout: block of dim rows per function (function r,
                    // direction d at row 2*r+d). See gsFunctionSet.h.
                    real_t sk = weights(k) * (derivs(2*r+0,k)*derivs(2*c+0,k) +
                                              derivs(2*r+1,k)*derivs(2*c+1,k));
                    K_gauss(i,j) += sk;
                    if (j != i) K_gauss(j,i) += sk;
                }
            }
        }
        time_gauss = sw.stop();
        gsInfo<<"Gauss assembly time: "<<time_gauss<<" s\n";
    }

    // 3b: Case 10 (rep-level dispatch) assembly with profiling
    {
        gsStopwatch sw_total;
        gsStopwatch sw_op;
        double time_ctr = 0, time_set = 0, time_mapto = 0, time_eval = 0, time_accum = 0;

        gsMatrix<real_t> M_c10 = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());
        gsMatrix<real_t> K_c10 = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());

        for (auto & elem : thb2_asm.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
            const auto & low = helem->lowerCorner();
            const auto & upp = helem->upperCorner();

            sw_op.restart();
            gsMatrix<real_t> ctr(2,1); ctr.col(0) = 0.5 * (low + upp);
            gsMatrix<index_t> actCtr; thb2_asm.active_into(ctr, actCtr);
            time_ctr += sw_op.stop();

            sw_op.restart();
            std::set<index_t> Lp_set;
            for (index_t r = 0; r != actCtr.rows(); ++r)
            for (index_t c = r; c != actCtr.rows(); ++c)
                Lp_set.insert(std::max(thb2_asm.repLevel(actCtr(r,0)), thb2_asm.repLevel(actCtr(c,0))));
            time_set += sw_op.stop();

            for (index_t Lp : Lp_set)
            {
                if (!rulesM_asm[Lp] || !rulesK_asm[Lp]) continue;

                sw_op.restart();
                gsMatrix<real_t> nodes; gsVector<real_t> weights;
                rulesM_asm[Lp]->mapTo(low, upp, nodes, weights);
                if (weights.size() == 0) continue;
                time_mapto += sw_op.stop();

                sw_op.restart();
                gsMatrix<index_t> acts; gsMatrix<real_t> vals; gsMatrix<real_t> derivs;
                thb2_asm.active_into(nodes, acts);
                thb2_asm.eval_into(nodes, vals);

                // Use reduced-regularity rule for derivatives (stiffness)
                gsMatrix<real_t> nodes_K; gsVector<real_t> weights_K;
                rulesK_asm[Lp]->mapTo(low, upp, nodes_K, weights_K);
                gsMatrix<index_t> acts_K; gsMatrix<real_t> derivs_K;
                if (weights_K.size() > 0)
                {
                    thb2_asm.active_into(nodes_K, acts_K);
                    thb2_asm.deriv_into(nodes_K, derivs_K);
                }
                time_eval += sw_op.stop();

                sw_op.restart();
                // Accumulate M with standard rules
                for (index_t r = 0; r != acts.rows(); ++r)
                for (index_t c = r; c != acts.rows(); ++c)
                {
                    const index_t i = acts(r,0), j = acts(c,0);
                    const index_t Li = thb2_asm.repLevel(i), Lj = thb2_asm.repLevel(j);
                    if (std::max(Li, Lj) != Lp) continue;

                    real_t sm = 0.0;
                    for (index_t k = 0; k != weights.size(); ++k)
                        sm += weights(k) * vals(r,k) * vals(c,k);
                    M_c10(i,j) += sm;
                    if (j != i) M_c10(j,i) += sm;
                }

                // Accumulate K with reduced-regularity rules
                if (weights_K.size() > 0)
                {
                    for (index_t r = 0; r != acts_K.rows(); ++r)
                    for (index_t c = r; c != acts_K.rows(); ++c)
                    {
                        const index_t i = acts_K(r,0), j = acts_K(c,0);
                        const index_t Li = thb2_asm.repLevel(i), Lj = thb2_asm.repLevel(j);
                        if (std::max(Li, Lj) != Lp) continue;

                        real_t sk = 0.0;
                        for (index_t k = 0; k != weights_K.size(); ++k)
                            sk += weights_K(k) * (derivs_K(2*r+0,k)*derivs_K(2*c+0,k) +
                                                  derivs_K(2*r+1,k)*derivs_K(2*c+1,k));
                        K_c10(i,j) += sk;
                        if (j != i) K_c10(j,i) += sk;
                    }
                }
                time_accum += sw_op.stop();
            }
        }
        double time_c10 = sw_total.stop();
        real_t err_M = (M_c10 - M_gauss).cwiseAbs().maxCoeff();
        real_t err_K = (K_c10 - K_gauss).cwiseAbs().maxCoeff();
        gsInfo<<"Case 10 assembly time: "<<time_c10<<" s (speedup: "<<time_gauss/time_c10<<"x)\n";
        gsInfo<<"  ||M_c10 - M_gauss||_inf = "<<err_M<<", ||K_c10 - K_gauss||_inf = "<<err_K<<"\n";
        gsInfo<<"  Breakdown: ctr="<<time_ctr<<" set="<<time_set<<" mapTo="<<time_mapto<<" eval="<<time_eval<<" accum="<<time_accum<<" s\n";
    }

    ////////////////////////////////////////////////////////////////////////
    // Priority-2 diagnostic: can one rule integrate the OTHER matrix?
    // Assemble BOTH M and K with a single given rule set (rep-level dispatch).
    //   - With the M rule (S(2p,C^r)): M exact, K should FAIL (C^r cannot
    //     represent the C^0 kink of the derivative product).
    //   - With the K rule (S(2p,C^{r-1})): K exact, and M should ALSO be exact
    //     (S(2p,C^{r-1}) contains the C^r mass integrand).
    ////////////////////////////////////////////////////////////////////////
    {
        auto assembleBoth = [&](const std::vector<gsQuadRule<real_t>::uPtr> & rset,
                                gsMatrix<real_t> & M, gsMatrix<real_t> & K)
        {
            M = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());
            K = gsMatrix<real_t>::Zero(thb2_asm.size(), thb2_asm.size());
            for (auto & elem : thb2_asm.domain()->allElements())
            {
                auto helem = dynamic_cast<gsHDomainIterator<real_t,2> *>(&elem);
                const auto & low = helem->lowerCorner();
                const auto & upp = helem->upperCorner();
                gsMatrix<real_t> ctr(2,1); ctr.col(0) = 0.5 * (low + upp);
                gsMatrix<index_t> actCtr; thb2_asm.active_into(ctr, actCtr);
                std::set<index_t> Lp_set;
                for (index_t r = 0; r != actCtr.rows(); ++r)
                for (index_t c = r; c != actCtr.rows(); ++c)
                    Lp_set.insert(std::max(thb2_asm.repLevel(actCtr(r,0)), thb2_asm.repLevel(actCtr(c,0))));
                for (index_t Lp : Lp_set)
                {
                    if (!rset[Lp]) continue;
                    gsMatrix<real_t> nodes; gsVector<real_t> weights;
                    rset[Lp]->mapTo(low, upp, nodes, weights);
                    if (weights.size() == 0) continue;
                    gsMatrix<index_t> acts; gsMatrix<real_t> vals, derivs;
                    thb2_asm.active_into(nodes, acts);
                    thb2_asm.eval_into(nodes, vals);
                    thb2_asm.deriv_into(nodes, derivs);
                    for (index_t r = 0; r != acts.rows(); ++r)
                    for (index_t c = r; c != acts.rows(); ++c)
                    {
                        const index_t i = acts(r,0), j = acts(c,0);
                        if (std::max(thb2_asm.repLevel(i), thb2_asm.repLevel(j)) != Lp) continue;
                        real_t sm = 0.0, sk = 0.0;
                        for (index_t k = 0; k != weights.size(); ++k)
                        {
                            sm += weights(k) * vals(r,k) * vals(c,k);
                            sk += weights(k) * (derivs(2*r+0,k)*derivs(2*c+0,k) +
                                                derivs(2*r+1,k)*derivs(2*c+1,k));
                        }
                        M(i,j) += sm; K(i,j) += sk;
                        if (j != i) { M(j,i) += sm; K(j,i) += sk; }
                    }
                }
            }
        };

        gsMatrix<real_t> M_fromM, K_fromM, M_fromK, K_fromK;
        assembleBoth(rulesM_asm, M_fromM, K_fromM);
        assembleBoth(rulesK_asm, M_fromK, K_fromK);
        gsInfo<<"Cross-rule diagnostic (rep-level dispatch):\n";
        gsInfo<<"  M with M-rule ||.-M_gauss|| = "<<(M_fromM-M_gauss).cwiseAbs().maxCoeff()
              <<"  (expect roundoff)\n";
        gsInfo<<"  K with M-rule ||.-K_gauss|| = "<<(K_fromM-K_gauss).cwiseAbs().maxCoeff()
              <<"  (expect FAIL: C^r rule misses C^0 kink)\n";
        gsInfo<<"  M with K-rule ||.-M_gauss|| = "<<(M_fromK-M_gauss).cwiseAbs().maxCoeff()
              <<"  (expect roundoff: K-rule degree 2p contains M integrand)\n";
        gsInfo<<"  K with K-rule ||.-K_gauss|| = "<<(K_fromK-K_gauss).cwiseAbs().maxCoeff()
              <<"  (expect roundoff)\n";
    }

    if (plot)
    {
        // Export Gauss points (value = weight).
        {
            const index_t n = static_cast<index_t>(allGX.size());
            gsMatrix<real_t> X(1,n), Y(1,n), Z(1,n), W(1,n);
            Z.setZero();
            for (index_t k = 0; k != n; ++k)
            { X(0,k)=allGX[k]; Y(0,k)=allGY[k]; W(0,k)=allGW[k]; }
            gsWriteParaviewPoints(X, Y, Z, W, "pts_gauss");
            gsInfo<<"Exported "<<n<<" Gauss quad points to pts_gauss.vts\n";
        }
        // Export patch-rule points per level (value = weight).
        for (index_t L = 0; L <= Lmax; ++L)
        {
            if (allPX[L].empty()) continue;
            const index_t n = static_cast<index_t>(allPX[L].size());
            gsMatrix<real_t> X(1,n), Y(1,n), Z(1,n), W(1,n);
            Z.setZero();
            for (index_t k = 0; k != n; ++k)
            { X(0,k)=allPX[L][k]; Y(0,k)=allPY[L][k]; W(0,k)=allPW[L][k]; }
            std::string fn = "pts_patch_L" + std::to_string(L);
            gsWriteParaviewPoints(X, Y, Z, W, fn);
            gsInfo<<"Exported "<<n<<" patch-rule points (L"<<L<<") to "<<fn<<".vts\n";
        }
    }

    /////////////////////////////////////////////////////////////////////////
    // Construct a basis 
    /////////////////////////////////////////////////////////////////////////
    gsTHBSplineBasis<2,real_t,true> thb2_mod(thb2);
    gsHElementHelper<2, real_t> helper(thb2_mod);

    // Loop over all truncated functions
    std::vector<index_t> refboxes;
    for (index_t i = 0; i!=thb2_mod.size(); i++)
    {
        // If the function is truncated, we compute its representation level and refine the element corresponding to the support of the function on that level
        if (thb2_mod.isTruncated(i))
        {
            index_t repLevel = thb2_mod.repLevel(i);
            gsMatrix<real_t> support = thb2_mod.support(i);
            // Convert support to element (box coordinates at repLevel)
            gsVector<real_t, 2> low = support.col(0);
            gsVector<real_t, 2> upp = support.col(1);
            
            auto element = helper.toElement(low, upp, repLevel);
            // Now you have element with box coordinates
            // Access them via:
            const auto & lowCorner = element.lowerCorner();  // Vector<index_t, 2>
            const auto & uppCorner = element.upperCorner();  // Vector<index_t, 2>
            gsInfo << "Element for truncated function " << i << " at representation level " << repLevel << ":\n";
            gsInfo << "  Low corner: " << lowCorner.transpose() << "\n";
            gsInfo << "  Upp corner: " << uppCorner.transpose() << "\n";

            std::vector<index_t> box(5);
            box[0] = repLevel; // level
            box[1] = lowCorner(0); // lower x
            box[2] = lowCorner(1); // lower y
            box[3] = uppCorner(0); // upper x
            box[4] = uppCorner(1); // upper y
            refboxes.push_back(box[0]);
            refboxes.push_back(box[1]);
            refboxes.push_back(box[2]);
            refboxes.push_back(box[3]);
            refboxes.push_back(box[4]);
        }
    }

    thb2_mod.refineElements(refboxes);

    if (plot)
    {
        gsMesh<> mesh(thb2_mod);
        gsWriteParaview(mesh,"mesh_refined");
        // gsWriteParaview(thb2_mod,"basis_refined",1000,true);
    }

     // for (auto & box : refboxes)
     // {
     //     gsInfo << "Refine box: level " << box[0] << ", low (" << box[1] << ", " << box[2] << "), upp (" << box[3] << ", " << box[4] << ")\n";


    return 0;
}