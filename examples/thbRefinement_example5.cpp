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

    boxes[0] = 4;
    boxes[1] = 6;
    boxes[2] = 2;
    boxes[3] = 8;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 10;
    boxes[2] = 10;
    boxes[3] = 12;
    boxes[4] = 12;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1000,true);


    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    gsHBoxContainer<2,real_t> markedRef, markedCrs;
    gsHBox<2,real_t> cell;
    basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    gsVector<index_t,2> low,upp;
    index_t lvl;

    low <<7,3;
    upp <<8,4;
    lvl = 4;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedRef.add(cell);

    low <<6,6;
    upp <<7,7;
    lvl = 3;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedRef.add(cell);

    low <<4,3;
    upp <<5,4;
    lvl = 3;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedRef.add(cell);

    markedRef.makeUnitBoxes();
    // markedRef.markTadmissible(m);

    gsDebugVar(markedRef);
    mp.patch(0).refineElements(markedRef.toRefBoxes());

    low <<0,0;
    upp <<1,1;
    lvl = 2;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedCrs.add(cell);

    low <<3,3;
    upp <<4,4;
    lvl = 3;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedCrs.add(cell);

    low <<5,1;
    upp <<6,2;
    lvl = 3;
    cell = gsHBox<2,real_t>(low,upp,lvl,basis);
    markedCrs.add(cell);


    for (typename gsHBoxContainer<2,real_t>::HIterator Hit=markedCrs.begin(); Hit!=markedCrs.end(); Hit++)
    {
        for (typename gsHBoxContainer<2,real_t>::cIterator Cit=Hit->begin(); Cit!=Hit->end(); Cit++)
        {
            bool clean = true;
            for (typename gsHBox<2,real_t>::HIterator HRit=markedRef.begin(); HRit!=markedRef.end(); HRit++)
            {
                for (typename gsHBoxContainer<2,real_t>::cIterator CRit=HRit->begin(); CRit!=HRit->end(); CRit++)
                {
                    gsHBox<2,real_t> box(*Cit);
                    gsHBox<2,real_t> parent(Cit->getParent());
                    gsDebugVar(box);
                    gsDebugVar(parent);
                    gsDebugVar(*CRit);

                    gsDebugVar(parent.contains(box));
                    gsDebugVar(box.isContained(parent));

                    gsDebugVar(parent.contains(*CRit));
                    gsDebugVar(CRit->isContained(box));
                    gsDebugVar(CRit->isContained(parent));
                }
            }
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

    // gsHBoxContainer<2,real_t> tmp;
    // gsHBoxContainer<2,real_t>::HContainer container;
    // if (verbose>1 &&   Hneigh  )
    // {
    //     gsInfo<<"------------------------H-Neighborhood m\n";
    //     tmp = gsHBoxContainer<2,real_t>(cell.getHneighborhood(m));
    //     container = tmp.boxes();
    //     for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
    //         for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
    //             gsInfo<<it->getCoordinates()<<"\n";

    // }

    // if (verbose>1 && !(Hneigh) )
    // {
    //     gsInfo<<"------------------------T-Neighborhood m\n";
    //     tmp = gsHBoxContainer<2,real_t>(cell.getTneighborhood(m));
    //     container = tmp.boxes();
    //     for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
    //         for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
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


    gsHBoxContainer<2,real_t> container(markedRef.toUnitBoxes());
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

