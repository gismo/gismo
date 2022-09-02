/** @file thbRefinement_example.cpp

    @brief Demonstates THB refinement and provides info on the resulting basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsHBoxUtils.h>
#include <gsAssembler/gsAdaptiveMeshing.h>
#include <gsAssembler/gsAdaptiveMeshingCompare.h>

using namespace gismo;

void writeCells(const std::string name, const gsMatrix<real_t> corners)
{
    std::ofstream file;
    file.open(name,std::ofstream::out);
    for (index_t c=0; c!=corners.cols()/2; c++)
    {
        for (index_t d=0; d!=2; d++)
            for (index_t r=0; r!=corners.rows(); r++)
            {
                file<<corners(r,2*c+d);
                if (!(r==corners.rows()-1 && d==1))
                    file<<",";
                else
                    file<<"\n";
            }
    }
    file.close();
}

int main(int argc, char *argv[])
{
    index_t degree    = 1;
    index_t m         = 2;

    index_t numHref   = 2;

    index_t testCase  = 0;
    index_t verbose   = 0;

    index_t rule = 3;

    bool plot     = false;
    bool Hneigh   = false;

    gsCmdLine cmd("Create standard refined THB meshes.");
    cmd.addInt("m","jump",
               "parameter m", m);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addInt("r","numHref",
               "Number of uniform refinements to be performed", numHref);
    cmd.addInt("t","testCase",
               "Test configuration", testCase);

    cmd.addInt("R","rule",
               "Rule for refinement/coarsening", rule);


    cmd.addInt("v","verbose",
               "Verbose output", verbose);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("Hneigh", "H-neighborhood if true, T-neighborhood if false (default)", Hneigh);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mpBspline, mp;

    gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    if (degree>1) bspline.degreeElevate(degree-1);

    mpBspline.addPatch(bspline);

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }

    std::vector<index_t> boxes(5);

    // Initial refinement
    boxes[0] = 1;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 2;
    boxes[4] = 2;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 4;
    boxes[4] = 2;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 2;
    boxes[2] = 2;
    boxes[3] = 4;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 3;
    boxes[1] = 2;
    boxes[2] = 0;
    boxes[3] = 6;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 3;
    boxes[1] = 4;
    boxes[2] = 4;
    boxes[3] = 8;
    boxes[4] = 8;
    mp.patch(0).refineElements(boxes);

    // boxes[0] = 4;
    // boxes[1] = 6;
    // boxes[2] = 2;
    // boxes[3] = 8;
    // boxes[4] = 4;
    // mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 10;
    boxes[2] = 10;
    boxes[3] = 12;
    boxes[4] = 12;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1000,true);


    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    gsHBoxContainer<2> markedRef, markedCrs;
    gsHBox<2> cell;
    basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    gsVector<index_t,2> low,upp;
    index_t lvl;

    low <<7,3;
    upp <<8,4;
    lvl = 4;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);

    low <<6,6;
    upp <<7,7;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);

    low <<4,3;
    upp <<5,4;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);

    markedRef.markTadmissible(m);
    markedRef.makeUnitBoxes();

    gsDebugVar(markedRef);
    gsWriteParaview(mp,"mp",1000,true);

    gsHBox<2>::Container c = markedRef.toContainer();
    for (typename gsHBox<2>::cIterator it=c.begin(); it!=c.end(); it++)
    {
        gsDebugVar(*it);
        gsHBoxContainer<2> siblings(it->getSiblings());
        gsDebugVar(siblings);
    }

    for (typename gsHBox<2>::Iterator it=c.begin(); it!=c.end(); it++)
    {
        gsDebugVar(*it);
        gsHBoxContainer<2> siblings(gsHBoxUnique<2,real_t>(it->getTneighborhood(2)));
        // siblings.makeUnique();
        gsDebugVar(siblings);
    }

    gsHBoxContainer<2> resContainer2(cell.getParent().getCneighborhood(m));
    gsDebugVar(resContainer2);



    low <<0,0;
    upp <<1,1;
    lvl = 2;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);

    low <<3,3;
    upp <<4,4;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);

    low <<5,1;
    upp <<6,2;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);

    low <<4,7;
    upp <<5,8;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);

    gsOverlapCompare<2,real_t> overlapCompare(markedRef.toContainer(),2);
    gsHBoxContainer<2> coarsened;
    for (typename gsHBoxContainer<2>::HIterator Hit=markedCrs.begin(); Hit!=markedCrs.end(); Hit++)
        for (typename gsHBoxContainer<2>::cIterator Cit=Hit->begin(); Cit!=Hit->end(); Cit++)
        {
            if (overlapCompare.check(*Cit))
                coarsened.add(*Cit);
        }

    gsDebugVar(coarsened);
    // mp.patch(0).refineElements(markedRef.toRefBoxes());
    // gsWriteParaview(mp,"mp_ref",1000,true);
    // mp.patch(0).unrefineElements(coarsened.toCrsBoxes());
    // gsWriteParaview(mp,"mp_crs",1000,true);
    return 0;
    //
    for (typename gsHBoxContainer<2>::HIterator Hit=markedCrs.begin(); Hit!=markedCrs.end(); Hit++)
        for (typename gsHBoxContainer<2>::cIterator Cit=Hit->begin(); Cit!=Hit->end(); Cit++)
        {
            for (typename gsHBoxContainer<2>::HIterator Hit2=markedCrs.begin(); Hit2!=markedCrs.end(); Hit2++)
            {
                for (typename gsHBoxContainer<2>::cIterator Cit2=Hit2->begin(); Cit2!=Hit2->end(); Cit2++)
                {
                    gsDebug<<"Checking "<<*Cit2<<" against "<<*Cit<<"\n";
                }
            }
        }


    typename gsBasis<>::domainIter domIt = mp.basis(0).makeDomainIterator();
    gsHDomainIterator<real_t,2> * domHIt = nullptr;
    domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
    std::vector<real_t> errors(gsMultiBasis<>(mp).totalElements());
    index_t i=0;
    for (; domHIt->good(); domHIt->next(), i++ )
    {
        errors[i] = i;
        gsHBox<2> box(domHIt);
        gsDebugVar(box);
        gsDebugVar(i);
    }

    gsDebugVar(errors.size());
    for (std::vector<real_t>::const_iterator it = errors.begin(); it!=errors.end(); it++)
        gsDebugVar(*it);

    errors[4]  = 1004;
    errors[28] = 1028;
    errors[35] = 1035;

    errors[21] = 1./1021;
    errors[31] = 1./1031;
    errors[41] = 1./1041;

    gsAdaptiveMeshing<real_t> mesher(mp);
    mesher.options().setInt("RefineRule",rule);
    mesher.options().setInt("CoarsenRule",rule);
    mesher.options().setReal("RefineParam",0.2);
    mesher.options().setReal("CoarsenParam",0.9999);
    mesher.getOptions();
    std::vector<bool> refine, coarsen;
    mesher.mark_into(errors,refine,false);
    mesher.mark_into(errors,coarsen,true);

    index_t j=0;
    for (std::vector<real_t>::const_iterator it = errors.begin(); it!=errors.end(); it++, j++)
    {
        std::string what;
        if (refine[j])
            what = "refined";
        else if (coarsen[j])
            what = "coarsened";
        else
            what = "nothing";
        gsDebug<<"Element "<<j<<": error = "<<errors[j]<<" is "<<what<<"\n";
    }

    std::vector<bool> filter(markedCrs.totalSize());
    // std::fill(filter.begin(), filter.end(), true);
    index_t k=0;
    for (typename gsHBoxContainer<2>::HIterator Hit=markedCrs.begin(); Hit!=markedCrs.end(); Hit++)
    {
        for (typename gsHBoxContainer<2>::cIterator Cit=Hit->begin(); Cit!=Hit->end(); Cit++, k++)
        {
            bool clean = true;

            // Rule #1: Support extension of the coarsening cell contains a refinement cell
            index_t l=0;
            for (typename gsHBoxContainer<2>::HIterator Hit2=markedCrs.begin(); Hit2!=markedCrs.end(); Hit2++)
            {
                for (typename gsHBoxContainer<2>::cIterator Cit2=Hit2->begin(); Cit2!=Hit2->end(); Cit2++, l++)
                {
                    if (!clean || Cit->isSame(*Cit2))
                        continue;

                    gsHBox<2> parent1(Cit->getParent());
                    gsHBox<2> parent2(Cit2->getParent());

                    clean &= !(parent1.contains(*Cit2));
                    if (!clean)
                    {
                        gsDebug<<"Rule #1: Box "<<*Cit<<" eliminated because it overlaps with coarsening box "<<*Cit2<<"\n";
                        break;
                    }
                }
            }

            if (!clean)
                continue;

            // Rule #2: Support extension of the coarsening cell contains a refinement cell ------>>>>> What about the T/H-Neighborhood???
            for (typename gsHBox<2>::HIterator HRit=markedRef.begin(); HRit!=markedRef.end(); HRit++)
            {
                for (typename gsHBoxContainer<2>::cIterator CRit=HRit->begin(); CRit!=HRit->end(); CRit++)
                {
                    if (!clean || Cit->isSame(*CRit))
                        continue;

                    gsHBox<2> box(*CRit);
                    gsHBox<2> parent(Cit->getParent());

                    clean &= !(parent.contains(box));
                    if (!clean)
                    {
                        gsDebug<<"Rule #2: Box "<<*Cit<<" eliminated because contains refinement box "<<*CRit<<"\n";
                        break;
                    }
                }
            }
            filter.at(k) = clean;
        }
    }


    gsWriteParaview(mp,"mp",1000,true);
    return 0;

    // Create a gsHBox for the marked cell
    if (verbose>0)
    {
        gsInfo<<"Cells marked for refinement:\n";
        gsInfo<<markedRef<<"\n";
    }

    // gsHBoxContainer<2> tmp;
    // gsHBoxContainer<2>::HContainer container;
    // if (verbose>1 &&   Hneigh  )
    // {
    //     gsInfo<<"------------------------H-Neighborhood m\n";
    //     tmp = gsHBoxContainer<2>(cell.getHneighborhood(m));
    //     container = tmp.boxes();
    //     for (typename gsHBoxContainer<2>::HIterator hit = container.begin(); hit!=container.end(); hit++)
    //         for (typename gsHBoxContainer<2>::Iterator it = hit->begin(); it!=hit->end(); it++)
    //             gsInfo<<it->getCoordinates()<<"\n";

    // }

    // if (verbose>1 && !(Hneigh) )
    // {
    //     gsInfo<<"------------------------T-Neighborhood m\n";
    //     tmp = gsHBoxContainer<2>(cell.getTneighborhood(m));
    //     container = tmp.boxes();
    //     for (typename gsHBoxContainer<2>::HIterator hit = container.begin(); hit!=container.end(); hit++)
    //         for (typename gsHBoxContainer<2>::Iterator it = hit->begin(); it!=hit->end(); it++)
    //             gsInfo<<it->getCoordinates()<<"\n";
    // }

    if (Hneigh)
        markedRef.markHadmissible(m);
    else
        markedRef.markTadmissible(m);

    if (verbose>1)
    {
        if (Hneigh)
        {
            gsInfo<<"------------------------Marked H-neighborhood\n";
            gsInfo<<markedRef;
        }
        else
        {
            gsInfo<<"------------------------Marked T-neighborhood\n";
            gsInfo<<markedRef;
        }
    }


    gsHBoxContainer<2> container(markedRef.toUnitBoxes());
    // gsDebugVar(container);

    std::vector<index_t> refboxes = container.toRefBoxes();
    for (index_t k=0; k!=refboxes.size()/5; k++)
    {
        for (index_t b=0; b!=5; b++)
            gsDebug<<refboxes[5*k+b]<<",";
        gsDebug<<"\n";
    }

    mp.patch(0).refineElements(container.toRefBoxes());
    gsInfo<<"Mesh has "<<mp.patch(0).coefs().rows()<<" DoFs\n";


    if (plot)
    {
        gsInfo<<"Plotting geometry...";
        gsWriteParaview(mp,"mp",1000,true);
        gsInfo<<"done\n";
        // gsInfo<<"Plotting hierarchical basis...";
        // gsWriteParaview(*basis,"basis",1000,true);
        // gsInfo<<"done\n";
        // gsInfo<<"Plotting individual levels...";
        // for (index_t k=0; k!=numHref+1; k++)
        // {
        //     gsInfo<<" "<<k;
        //     gsWriteParaview((basis->tensorLevel(k)),"basis_" + std::to_string(k),1000);
        // }
        // gsInfo<<"done\n";

    }

    return 0;
}

