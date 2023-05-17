/** @file gsHBox_example

    @brief Demonstrates functionality of the gsHBox

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (TU Delft 2019-...)
*/

#include <iostream>

#include <gismo.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsHBoxUtils.h>
#include <gsAssembler/gsAdaptiveMeshing.h>
#include <gsAssembler/gsAdaptiveMeshingCompare.h>

using namespace gismo;

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

    gsHBoxContainer<2> markedRef, markedRef2, markedCrs;
    gsHBox<2> cell;

    gsVector<index_t,2> low,upp;
    index_t lvl;
    low <<7,5;
    upp <<8,6;
    lvl = 4;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);
    gsWriteParaview(markedRef,"markedRef");

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

    markedRef.markAdmissible(2);
    gsWriteParaview(markedRef,"refCell_Admissible");

    return 0;
}

