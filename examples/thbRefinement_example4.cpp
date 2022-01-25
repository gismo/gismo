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

    index_t steps     = 1;

    bool plot     = false;
    bool Hneigh   = false;

    gsCmdLine cmd("Create standard refined THB meshes.");
    cmd.addInt("m","jump",
               "parameter m", m);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addInt("r","numHref",
               "Number of uniform refinements to be performed", numHref);
    cmd.addInt("s","steps",
               "Number of steps to be performed", steps);

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

    if (testCase==1)
    {
        mpBspline.uniformRefine();
        mpBspline.uniformRefine();
    }

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


    std::vector<index_t> boxes(5), marked(5);

    // Initial refinement
    if (testCase==0)
    {
        boxes[0] = 1;
        boxes[1] = 0;
        boxes[2] = 0;
        boxes[3] = 1;
        boxes[4] = 1;
        mp.patch(0).refineElements(boxes);
    }

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    for (index_t s = 0; s!=steps; s++)
    {
        gsHBoxContainer<2,real_t> marked;
        gsHBox<2,real_t> cell;
        basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

        if (verbose>0) gsInfo<<"-----------------step "<<s<<"-----------------\n";
        gsVector<index_t,2> low,upp;
        index_t lvl = -1;
        if (testCase==0)
        {
            lvl = s+1;
            low<<0,0;
            upp<<1,1;
            cell = gsHBox<2,real_t>(low,upp,lvl,basis);
            marked.add(cell);
        }
        else if (testCase==1)
        {
            lvl = s;
            index_t nb = basis->tensorLevel(basis->maxLevel()).knots(0).uSize()-1;
            index_t extent = 2*math::ceil((degree+1)/2.)-1;
            for ( index_t i = 0; i!= nb; ++i)
            {
                low<<std::max(i-extent/2,0),   std::max(i,0);
                upp<<std::min(i+1+extent/2,nb),std::min(i+1,nb);

                cell = gsHBox<2,real_t>(low,upp,lvl,basis);
                marked.add(cell);
            }
            marked.makeUnitBoxes();
        }

        // Create a gsHBox for the marked cell
        if (verbose>0)
        {
            gsInfo<<"Cells marked for refinement:\n";
            gsInfo<<marked<<"\n";
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
            marked.markHadmissible(m);
        else
            marked.markTadmissible(m);

        if (verbose>1)
        {
            if (Hneigh)
            {
                gsInfo<<"------------------------Marked H-neighborhood\n";
                gsInfo<<marked;
            }
            else
            {
                gsInfo<<"------------------------Marked T-neighborhood\n";
                gsInfo<<marked;
            }
        }


        gsHBoxContainer<2,real_t> container(marked.toUnitBoxes());
        // gsDebugVar(container);

        mp.patch(0).refineElements(container.toRefBoxes());
        gsInfo<<"Mesh has "<<mp.patch(0).coefs().rows()<<" DoFs\n";

    }

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

