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

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t degree    = 1;
    index_t m         = 2;
    gsCmdLine cmd("Example of gsHBox.");
    cmd.addInt("m","jump",
               "parameter m", m);
    cmd.addInt("p","degree",
               "Spline degree", degree);
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
    gsInfo<<"Initial basis constructed:\n"<<mp.basis(0)<<"\n";

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
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
    gsInfo<<"Added one refinement box (markedRef.pvd):\n"<<markedRef<<"\n";

    markedRef.markAdmissible(m);
    gsWriteParaview(markedRef,"refCell_Admissible");
    gsInfo<<"Admissibility region of the refined cell plotted in refCell_Admissible.pvd\n";

    low <<4,2;
    upp <<5,3;
    lvl = 3;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedCrs.add(cell);
    gsWriteParaview(cell,"crsCell");
    gsInfo<<"Added one coarsening box (markedCrs.pvd):\n"<<markedCrs<<"\n";

    gsHBoxContainer<2> Cextension(cell.getParent().getCextension(m));
    gsWriteParaview(Cextension,"crsCell_Cextension");
    gsInfo<<"Coarsening extension plotted in crsCell_Cextension.pvd\n";

    gsHBoxContainer<2> Cneighborhood(cell.getParent().getCneighborhood(m));
    gsWriteParaview(Cneighborhood,"crsCell_Cneighborhood");
    gsInfo<<"Coarsening neighborhood plotted in crsCell_Cneighborhood.pvd\n";
    return 0;
}

