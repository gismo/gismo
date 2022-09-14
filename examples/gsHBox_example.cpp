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

typename gsHBox<2,real_t>::HContainer markRecursive(const typename gsHBox<2,real_t>::HContainer & marked, index_t lvl, index_t m)
{
    typename gsHBox<2,real_t>::HContainer marked_copy = marked;
    typename gsHBox<2,real_t>::Container marked_l = marked[lvl];
    typename gsHBox<2,real_t>::Container marked_k;

    gsHBoxContainer<2,real_t> neighbors;
    for (typename gsHBox<2,real_t>::Iterator it = marked_l.begin(); it!=marked_l.end(); it++)
    {
        neighbors.add(it->getTneighborhood(m));
    }
    gsDebugVar(neighbors.totalSize());
    index_t k = lvl - m + 1;
    if (neighbors.boxes().size()!=0)
    {
        marked_k = marked_copy[k];
        gsHBoxContainer<2,real_t> boxUnion = gsHBoxUtils<2,real_t>::Union(neighbors,gsHBoxContainer<2,real_t>(marked_k));
        marked_copy[k] = boxUnion.getActivesOnLevel(k);
        marked_copy = markRecursive(marked_copy,k,m);
    }
    return marked_copy;
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
    boxes[1] = 0;
    boxes[2] = 2;
    boxes[3] = 2;
    boxes[4] = 4;
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
    boxes[2] = 4;
    boxes[3] = 8;
    boxes[4] = 6;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 10;
    boxes[2] = 12;
    boxes[3] = 12;
    boxes[4] = 14;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1,true);


    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    gsHBoxContainer<2> markedRef, markedRef2, markedCrs;
    gsHBox<2> cell;
    basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    gsVector<index_t,2> low,upp;
    index_t lvl;
    // low <<10,13;
    // upp <<11,14;
    // lvl = 4;
    // cell = gsHBox<2>(low,upp,lvl,basis);
    // markedRef.add(cell);
    //
    low <<7,5;
    upp <<8,6;
    lvl = 4;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);
    gsDebugVar(markedRef);
    gsHBoxContainer<2> tmp(gsHBoxUtils<2,real_t>::markAdmissible(cell,m));
    gsDebugVar(tmp);

    // low <<6,6;
    // upp <<7,7;
    // lvl = 3;
    // cell = gsHBox<2>(low,upp,lvl,basis);
    // markedRef.add(cell);

    // low <<4,3;
    // upp <<5,4;
    // lvl = 3;
    // cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);
    gsWriteParaview(markedRef,"markedRef");
    // markedRef2 = markedRef;

    // markedRef.markTadmissible(m);
    // gsWriteParaview(markedRef,"markedRef_T-admissible");

    // markedRef2.markHadmissible(m);
    // gsWriteParaview(markedRef2,"markedRef_H-admissible");

    // markedRef.makeUnitBoxes();

    // gsDebugVar(markedRef);

    // low <<0,0;
    // upp <<1,1;
    // lvl = 2;
    // cell = gsHBox<2>(low,upp,lvl,basis);
    // markedCrs.add(cell);

    low <<4,2;
    upp <<5,3;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);

    gsWriteParaview(cell,"crsCell");
    gsHBoxContainer<2> Cextension(cell.getParent().getCextension(m));
    gsWriteParaview(Cextension,"crsCell_Cextension");

    gsHBoxContainer<2> Cneighborhood(cell.getParent().getCneighborhood(m));
    gsWriteParaview(Cneighborhood,"crsCell_Cneighborhood");

    gsHBox<2,real_t>::Container neighborhood, markedRedContainer = markedRef.toContainer();
    for (typename gsHBox<2,real_t>::Iterator mit=markedRedContainer.begin(); mit!=markedRedContainer.end(); mit++)
    {
        index_t lvl = mit->level();
        index_t k = lvl - m + 2;
        if (k-1>=0)
        {
            // Get multi level support extension on level k
            gsHBox<2,real_t>::Container extension = mit->getMultiLevelSupportExtension(k);
            gsHBoxContainer<2,real_t> extensionContainer(extension);
            gsWriteParaview(extensionContainer,"extension");
            // Eliminate elements which are too low
            gsHBox<2,real_t>::Container parents = extensionContainer.getParents();

            for (gsHBox<2,real_t>::Iterator it = parents.begin(); it!=parents.end(); it++)
            {
                it->computeCenter(); // needed to check active
                if (it->isActive())
                    neighborhood.push_back(*it);
            }
            // neighborhood = parents;
        }
    }
    gsHBoxContainer<2,real_t> neighborhoodContainer(neighborhood);
    gsWriteParaview(neighborhoodContainer,"neighborhood");


    typename gsHBoxContainer<2,real_t>::HContainer unitBoxes = markedRef.toUnitHBoxes();
    gsHBoxContainer<2,real_t> plt(unitBoxes);
    gsWriteParaview(plt,"unitBoxes_" + std::to_string(0));
    for (size_t l = 0; l!=unitBoxes.size(); l++)
    {
        unitBoxes = markRecursive(unitBoxes,l,m);
        plt = gsHBoxContainer<2,real_t> (unitBoxes);
        gsWriteParaview(plt,"unitBoxes_" + std::to_string(l+1));
    }

    unitBoxes = gsHBoxUtils<2,real_t>::Unique(unitBoxes);

    return 0;
}

