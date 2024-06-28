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
    boxes[1] = 10;
    boxes[2] = 12;
    boxes[3] = 12;
    boxes[4] = 14;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1,true);


    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    // Assigns errors to the elements, based on their level: e = 10^{l-1}
    typename gsBasis<>::domainIter domIt = mp.basis(0).makeDomainIterator();
    gsHDomainIterator<real_t,2> * domHIt = nullptr;
    domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
    std::vector<real_t> errors(gsMultiBasis<>(mp).totalElements());
    index_t i=0;
    for (; domHIt->good(); domHIt->next(), i++ )
        errors[i] = math::pow(10,domHIt->getLevel()-1);// + i;

    index_t offset = 0;
    for (size_t p = 0; p!=mp.nPatches(); p++)
    {
        auto first = errors.begin() + offset;
        offset += mp.basis(p).numElements();
        auto last = errors.begin() + offset;
        std::vector<real_t> tmpErrors(first,last);
        gsElementErrorPlotter<real_t> err_eh(mp.basis(p),tmpErrors);
        const gsField<> elemError_eh( mp.patch(p), err_eh, true );
        gsWriteParaview<>( elemError_eh, "error_elem_ref" + std::to_string(p), 10000, false);
    }

    gsAdaptiveMeshing<real_t> mesher(mp);
    mesher.options().setInt("RefineRule",rule);
    mesher.options().setInt("CoarsenRule",rule);
    mesher.options().setSwitch("Admissible",true);
    if (rule==1)
    {
        mesher.options().setReal("RefineParam",0.3);
        mesher.options().setReal("CoarsenParam",0.1);
    }
    else if (rule==2)
    {
        mesher.options().setReal("RefineParam",0.3);
        mesher.options().setReal("CoarsenParam",0.1);
    }
    else if (rule==3)
    {
        mesher.options().setReal("RefineParam",0.1);
        mesher.options().setReal("CoarsenParam",0.01);
    }

    mesher.getOptions();
    gsHBoxContainer<2,real_t> refine,coarsen;
    mesher.markRef_into(errors,refine);

    gsVector<index_t,2> low,upp;
    index_t lvl;
    low <<4,5;
    upp <<5,6;
    lvl = 3;
    gsHBox<2>cell(low,upp,lvl,basis);

    refine.add(cell);
    mesher.markCrs_into(errors,refine,coarsen);

    gsInfo<<"Cells marked for refinement:\n";
    gsInfo<<refine<<"\n";
    gsWriteParaview(refine,"marked4ref");

    gsInfo<<"Cells marked for coarsening:\n";
    gsInfo<<coarsen<<"\n";
    gsWriteParaview(coarsen,"marked4crs");

    mesher.refine(refine);
    mesher.unrefine(coarsen);
    gsWriteParaview(mp,"end",1,true);

    return 0;
}

