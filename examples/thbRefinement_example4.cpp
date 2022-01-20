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
    index_t refmode   = 0;
    index_t degree    = 2;

    index_t numHref   = 2;

    index_t testCase  = 0;

    index_t steps     = 1;

    bool plot     = false;

    gsCmdLine cmd("Create standard refined THB meshes.");
    cmd.addInt("m","mode",
               "Refinement mode (0, 1, 2 3 4)", refmode);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addInt("r","numHref",
               "Number of uniform refinements to be performed", numHref);
    cmd.addInt("s","steps",
               "Number of steps to be performed", steps);

    cmd.addInt("t","testCase",
               "Test configuration", testCase);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mpBspline, mp;

    gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    bspline.degreeElevate(degree);

    gsVector<> coords;
    if (testCase==0)
    {
        coords.resize(1);
        coords<<0.5;
    }

    for (index_t k=0; k!=coords.size(); k++)
    {
        bspline.insertKnot(coords.at(k),0);
        bspline.insertKnot(coords.at(k),1);
    }

    mpBspline.addPatch(bspline);

    gsInfo<<"mpBspline.patch(0) = \n";
    gsInfo<<mpBspline.patch(0)<<"\n";

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }


    index_t m;
    std::vector<index_t> boxes(5), marked(5);

    // Initial refinement
    if (testCase==0)
    {
        boxes[0] = 1;
        boxes[1] = 0;
        boxes[2] = 0;
        boxes[3] = 1;
        boxes[4] = 1;
    }
    mp.patch(0).refineElements(boxes);

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    for (index_t s = 0; s!=steps; s++)
    {
        gsInfo<<"-----------------step "<<s<<"-----------------\n";
        gsVector<index_t,2> low,upp;
        index_t lvl = -1;
        if (testCase==0)
        {
            m = 2;
            lvl = s+1;
            low<<0,0;
            upp<<1,1;
        }

        basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

        // /// Write the cells

        // for (index_t l=0; l!=basis->maxLevel()+1; l++)
        // {
        //     typename gsBasis<real_t>::domainIter domIt = basis->tensorLevel(l).makeDomainIterator();
        //     gsMatrix<real_t> cells(2,2*domIt->numElements());
        //     index_t c = 0;
        //     for (; domIt->good(); domIt->next() )
        //     {
        //         cells.col(2*c)   = domIt->lowerCorner();
        //         cells.col(2*c+1) = domIt->upperCorner();
        //         c++;
        //     }
        //     writeCells("level" + std::to_string(l) + ".csv",cells);
        // }

        // typename gsBasis<real_t>::domainIter domIt = basis->makeDomainIterator();
        // gsHDomainIterator<real_t,2> * domHIt = nullptr;
        // domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
        // GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

        // gsMatrix<real_t> hcells(2,2*domHIt->numElements());
        // index_t c = 0;

        // gsVector<index_t> idx;
        // std::vector<index_t> parent;
        // std::vector<std::vector<index_t>> suppext;
        // for (; domHIt->good(); domHIt->next() )
        // {
        //     hcells.col(2*c)   = domHIt->lowerCorner();
        //     hcells.col(2*c+1) = domHIt->upperCorner();
        //     c++;
        // }
        // writeCells("hcells.csv",hcells);

        // Create a gsHBox for the marked cell
        gsHBox<2,real_t> cell(low,upp,lvl,basis);
        gsMatrix<> cellCoords = cell.getCoordinates();
        // writeCells("selected.csv",cellCoords);
        gsDebugVar(cellCoords);
        gsDebugVar(cell);

        gsHBoxContainer<2,real_t> tmp;
        gsHBoxContainer<2,real_t>::HContainer container;

        gsDebugVar("------------------------Support Extension");
        tmp = gsHBoxContainer<2,real_t>(cell.getSupportExtension());
        container = tmp.boxes();
        for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
            for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                gsDebugVar(it->getCoordinates());

        gsDebugVar("------------------------Multi-Level Support Extension 0");
        tmp = gsHBoxContainer<2,real_t>(cell.getMultiLevelSupportExtension(0));
        container = tmp.boxes();
        for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
            for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                gsDebugVar(it->getCoordinates());

        gsDebugVar("------------------------Multi-Level Support Extension 1");
        tmp = gsHBoxContainer<2,real_t>(cell.getMultiLevelSupportExtension(1));
        container = tmp.boxes();
        for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
            for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                gsDebugVar(it->getCoordinates());

        gsDebugVar("------------------------H-Neighborhood 2");
        tmp = gsHBoxContainer<2,real_t>(cell.getHneighborhood(m));
        container = tmp.boxes();
        for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
            for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                gsDebugVar(it->getCoordinates());

        gsDebugVar("------------------------T-Neighborhood 2");
        tmp = gsHBoxContainer<2,real_t>(cell.getTneighborhood(m));
        container = tmp.boxes();
        for (typename gsHBoxContainer<2,real_t>::HIterator hit = container.begin(); hit!=container.end(); hit++)
            for (typename gsHBoxContainer<2,real_t>::Iterator it = hit->begin(); it!=hit->end(); it++)
                gsDebugVar(it->getCoordinates());

        gsHBoxContainer<2,real_t> marked(cell);
        for (index_t l = 0; l!=basis->maxLevel()+1; l++)
        {
            gsDebugVar(l);
            marked.markHrecursive(l,m);
        }

        gsDebugVar(basis->maxLevel());
        gsDebugVar(marked);

        mp.patch(0).refineElements(cell.toRefBox());

    }

    if (plot)
    {
        gsInfo<<"Plotting geometry...";
        gsWriteParaview(mp,"mp",1000,true);
        gsInfo<<"done\n";
        gsInfo<<"Plotting hierarchical basis...";
        gsWriteParaview(*basis,"basis",1000,true);
        gsInfo<<"done\n";
        gsInfo<<"Plotting individual levels...";
        for (index_t k=0; k!=numHref+1; k++)
        {
            gsInfo<<" "<<k;
            gsWriteParaview((basis->tensorLevel(k)),"basis_" + std::to_string(k),1000);
        }
        gsInfo<<"done\n";

    }

    return 0;
}

