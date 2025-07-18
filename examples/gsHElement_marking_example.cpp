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
#include <gsHSplines/gsHElementMarker.h>

// #include <gsUtils/gsCombinatorics.h>

using namespace gismo;


template<short_t d>
void run(gsTensorBSpline<d,real_t> & spline, index_t degree, index_t m, index_t numRef);

int main(int argc, char *argv[])
{

    index_t degree    = 1;
    index_t m         = 2;
    index_t numRef    = 5;
    short_t d         = 2;
    gsCmdLine cmd("Example of gsHBox.");
    cmd.addInt("m","jump",
               "parameter m", m);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addInt("n","numRef",
               "Number of refinement iterations", numRef);
    cmd.addInt("d","dim",
               "Dimension of the problem", d);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (d==2)
    {
        gsTensorBSpline<2,real_t> spline = *gsNurbsCreator<>::BSplineSquare(1,0,0);
        run(spline, degree, m, numRef);
    }
    else if (d==3)
    {
        gsTensorBSpline<3,real_t> spline = *gsNurbsCreator<>::BSplineCube(1,0,0,0);
        run(spline, degree, m, numRef);
    }
    return EXIT_SUCCESS;
}

template<short_t d>
void run(gsTensorBSpline<d,real_t> & spline, index_t degree, index_t m, index_t numRef)
{
    typedef typename gsHElementHelper<d,real_t>::element_t element_t;
    typedef typename gsHElementHelper<d,real_t>::HElementContainer HElementContainer;
    gsMultiPatch<> mpBspline, mp;

    if (degree>1) spline.degreeElevate(degree-1);

    // for (index_t k = 0; k < 6; ++k)
    //     bspline.uniformRefine();

    mpBspline.addPatch(spline);

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<d,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<d,real_t> *geo = dynamic_cast< gsTensorBSpline<d,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<d,real_t>(*geo);
        mp.addPatch(thb);
    }

    gsHElementMarker<d,real_t> marker(mp.basis(0));
    marker.options().setInt("MaxLevel",10);
    marker.options().setInt("RefineRule",1);
    marker.options().setReal("RefineParam",0.99);
    marker.options().setInt("CoarsenRule",1);
    marker.options().setReal("CoarsenParam",0.1);

    gsParaviewCollection meshes("mesh");
    gsParaviewCollection refined("markedRef");
    gsParaviewCollection coarsened("markedCrs");
    gsStopwatch timer;
    for (index_t i = 0; i!=numRef; i++)
    {
        // Plot stuff
        gsMesh<real_t> mesh(mp.basis(0));
        gsWriteParaview(mesh,"mesh_"+util::to_string(i),false);
        meshes.addPart("mesh_"+util::to_string(i)+".vtp",i,"Mesh");

        gsInfo<<"Refinement iteration "<<i<<":\n";
        gsInfo<<"Number of elements: "<<mp.basis(0).numElements()<<"\n";

        // Make a random error vector, for testing purposes
        gsVector<real_t> errorVec(mp.basis(0).numElements());
        errorVec.setRandom();
        errorVec.array() += 1.0; // Make sure the errors are positive
        errorVec.array() *= 0.5; // Scale the errors down
        std::vector<real_t> errs(errorVec.data(), errorVec.data() + errorVec.size());

        marker.setErrors(errs);
        // gsInfo<<marker<<"\n";

        timer.restart();
        HElementContainer markedRef = marker.markRef();
        gsInfo<<"Marked "<<markedRef.size()<<" elements for refinement in "<<timer.stop()<<" seconds.\n";
        // gsInfo<<"Marked elements:\n";
        // for (const auto & elem : markedRef)
        //     gsInfo<<elem<<"\n";

        timer.restart();
        HElementContainer markedCrs = marker.markCrs(markedRef);
        gsInfo<<"Marked "<<markedCrs.size()<<" elements for coarsening in "<<timer.stop()<<" seconds.\n";
        // gsInfo<<"Marked elements:\n";
        // for (const auto & elem : markedCrs)
        //     gsInfo<<elem<<"\n";

        gsMatrix<> boxes;
        gsVector<size_t> levels;
        std::tie(boxes,levels) = marker.helper().toBoxesAndLevels(markedRef);
        gsWriteParaview(boxes,"markedRef_"+util::to_string(i),gsVector<real_t>(levels.cast<real_t>()));
        refined.addPart("markedRef_"+util::to_string(i)+".vtu",i,"Solution");
        std::tie(boxes,levels) = marker.helper().toBoxesAndLevels(markedCrs);
        gsWriteParaview(boxes,"markedCrs_"+util::to_string(i),gsVector<real_t>(levels.cast<real_t>()));
        coarsened.addPart("markedCrs_"+util::to_string(i)+".vtu",i,"Solution");

        std::vector<index_t> refBox = marker.toRefBoxes(markedRef);
        mp.basis(0).refineElements(refBox);
        std::vector<index_t> crsBox = marker.toCrsBoxes(markedCrs);
        mp.basis(0).unrefineElements(crsBox);
    }

    meshes.save();
    refined.save();
    coarsened.save();
}

