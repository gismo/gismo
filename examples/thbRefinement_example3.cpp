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
    index_t refLevels = 5;
    index_t refmode   = 0;
    index_t degree    = 2;
    index_t CExt      = 0;

    index_t numHref   = 2;


    real_t xmin = 0.0;
    real_t xmax = 0.0625;
    real_t ymin = 0.0;
    real_t ymax = 0.0625;

    bool plot     = false;

    gsCmdLine cmd("Create standard refined THB meshes.");
    cmd.addInt("E","extension",
               "Coarsening extension", CExt);
    cmd.addInt("l","levels",
               "Number of refinement levels", refLevels);
    cmd.addInt("m","mode",
               "Refinement mode (0, 1, 2 3 4)", refmode);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    cmd.addInt("r","numHref",
               "Number of uniform refinements to be performed", numHref);

    cmd.addReal("x","xmin","xmin", xmin);
    cmd.addReal("y","ymin","ymin", ymin);
    cmd.addReal("X","xmax","xmax", xmax);
    cmd.addReal("Y","ymax","ymax", ymax);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mpBspline, mp;
    mpBspline.addPatch(gsNurbsCreator<>::BSplineSquare(1,0,0));
    mpBspline.degreeElevate(degree);



    // // Create a tensor-producte basis
    // gsKnotVector<> KV (0, 1, numknots, degree+1, 1);
    // gsTensorBSplineBasis<2> tp(KV,KV);
    // gsInfo<< "Coarsest level: "<< tp <<"\n";

    gsTensorBSplineBasis<2> & tp = dynamic_cast<gsTensorBSplineBasis<2> & >(mpBspline.basis(0));
    gsKnotVector<> KV = tp.knots(0);
    index_t numknots = KV.uSize() - 2;
    gsDebugVar(KV);
    gsDebugVar(numknots);

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }

    std::vector<index_t> boxes(5);
    boxes[0] = numHref;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = math::pow(2,numHref);
    boxes[4] = math::pow(2,numHref);

    // mp.uniformRefine();
    // mp.uniformRefine();

    // std::vector<index_t> boxes(5);
    // boxes[0] = 2;
    // boxes[1] = 2;
    // boxes[2] = 2;
    // boxes[3] = 8;
    // boxes[4] = 8;



    for (index_t k=0; k!=5; k++)
        gsDebugVar(boxes[k]);

    mp.patch(0).refineElements(boxes);

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    for (index_t l=0; l!=basis->maxLevel()+1; l++)
    {
        typename gsBasis<real_t>::domainIter domIt = basis->tensorLevel(l).makeDomainIterator();
        gsMatrix<real_t> cells(2,2*domIt->numElements());
        index_t c = 0;
        for (; domIt->good(); domIt->next() )
        {
            cells.col(2*c)   = domIt->lowerCorner();
            cells.col(2*c+1) = domIt->upperCorner();
            c++;
        }
        writeCells("level" + std::to_string(l) + ".csv",cells);
    }

    typename gsBasis<real_t>::domainIter domIt = basis->makeDomainIterator();
    gsHDomainIterator<real_t,2> * domHIt = nullptr;
    domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
    GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

    gsMatrix<real_t> hcells(2,2*domHIt->numElements());
    index_t c = 0;

    gsVector<index_t> idx;
    std::vector<index_t> parent;
    std::vector<std::vector<index_t>> suppext;
    for (; domHIt->good(); domHIt->next() )
    {
        hcells.col(2*c)   = domHIt->lowerCorner();
        hcells.col(2*c+1) = domHIt->upperCorner();
        c++;
    }
    writeCells("hcells.csv",hcells);

    domHIt->reset();
    index_t cellIdx = 5;
    for (index_t k=0; k!=cellIdx; k++)
        domHIt->next();

    gsHBox<2,real_t> cell(domHIt);
    gsMatrix<> cellCoords = cell.getCoordinates();

    writeCells("selected.csv",cellCoords);

    gsHBoxContainer<2,real_t> tmp;
    gsHBox<2,real_t> box1(domHIt);
    gsDebugVar(box1);
    gsDebugVar(box1.getCoordinates());
    gsDebugVar(box1.getParent());
    gsDebugVar(box1.getParent().getParent());
    gsDebugVar(box1.getAncestor(1));
    gsHBox<2,real_t> box2(box1.lowerIndex(),box1.upperIndex(),box1.level());
    gsDebugVar(box2);
    gsHBox<2,real_t> box3(box1.lowerIndex(),box1.upperIndex(),box1.level(),basis);
    gsDebugVar(box3);

    gsDebugVar(box1.contains(box1));

    tmp = gsHBoxContainer<2,real_t>(box1.getSupportExtension());
    gsDebugVar(tmp);

    tmp = gsHBoxContainer<2,real_t>(box1.getMultiLevelSupportExtension(1));
    gsDebugVar(tmp);
    tmp = gsHBoxContainer<2,real_t>(box1.getMultiLevelSupportExtension(2));
    gsDebugVar(tmp);
    tmp = gsHBoxContainer<2,real_t>(box1.getMultiLevelSupportExtension(3));
    gsDebugVar(tmp);

    tmp = gsHBoxContainer<2,real_t>(box1.getHneighborhood(2));
    gsDebugVar(tmp);
    tmp = gsHBoxContainer<2,real_t>(box1.getTneighborhood(2));
    gsDebugVar(tmp);
    // tmp = gsHBoxContainer<2,real_t>(box1.toHBoxContainer().markHrecursive(3,2));
    // gsDebugVar(tmp);


    gsHBoxContainer<2,real_t> HBContainer1;
    for (index_t k=0; k!=5 && domHIt->good(); k++)
    {
        HBContainer1.add(gsHBox<2,real_t>(domHIt));
        domHIt->next();
    }
    gsDebugVar(HBContainer1);


    gsHBoxContainer<2,real_t> HBContainer2;
    for (index_t k=0; k!=5 && domHIt->good(); k++)
    {
        HBContainer2.add(gsHBox<2,real_t>(domHIt));
        domHIt->next();
    }
    gsDebugVar(HBContainer2);

    gsHBoxContainer<2,real_t> HBContainer3 = HBContainer1.boxUnion(HBContainer2);
    gsDebugVar(HBContainer3);

    gsHBoxContainer<2,real_t> HBContainer4(HBContainer3.getParents());
    gsDebugVar(HBContainer4);

    gsHBoxContainer<2,real_t> marked(cell);
    gsHBoxContainer<2,real_t> HBContainer5(marked.markHrecursive(marked.boxes(),cell.level(),1));
    gsDebugVar(HBContainer5);


    // CellContainerType H_cells = get_H_neighborhood<2>(basis,domHIt,1);
    // gsMatrix<> Hneighborhood;
    // std::tie(Hneighborhood,std::ignore) = getCoordinates<2>(basis,H_cells);
    // writeCells("Hneighborhood.csv",Hneighborhood);

    // CellContainerType T_cells = get_T_neighborhood<2>(basis,domHIt,1);
    // gsMatrix<> Tneighborhood;
    // std::tie(Tneighborhood,std::ignore) = getCoordinates<2>(basis,T_cells);
    // writeCells("Tneighborhood.csv",Tneighborhood);



    // typename gsBasis<real_t>::domainIter domIt2 = basis->makeDomainIterator();
    // gsHDomainIterator<real_t,2> * domHIt2 = nullptr;
    // domHIt2 = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt2.get());
    // GISMO_ENSURE(domHIt2!=nullptr,"Domain should be 2 dimensional for flattening");
    // index_t m = 1;
    // // get_multilevel_support_extension
    // for (; domHIt2->good(); domHIt2->next() )
    // {
    //     const gsVector<real_t>& cen = domHIt2->centerPoint();

    //     gsMatrix<index_t> acts;
    //     basis->tensorLevel(m).active_into(cen,acts);

    //     gsMatrix<real_t> boxes(2,2*acts.rows());
    //     for (index_t act = 0; act!=acts.rows(); act++)
    //         boxes.block(0,2*act,2,2) = basis->tensorLevel(m).support(acts(act,0));

    //     gsDebugVar(boxes);
    //     gsDebugVar(cen.transpose());
    //     gsDebugVar(acts.transpose());
    // }


    gsWriteParaview(mp,"mp",1000,true);
    gsWriteParaview(*basis,"basis",1000,true);

    for (index_t k=0; k!=numHref+1; k++)
        gsWriteParaview((basis->tensorLevel(k)),"basis_" + std::to_string(k),1000);

    gsWrite(mp,"mp");


    return 0;
}

