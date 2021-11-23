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


template<int d>
std::pair<gsMatrix<real_t>,std::vector<index_t>> getCoordinates(const gsHTensorBasis<d,real_t> * basis, const std::vector<index_t> boxes)
{
    GISMO_ASSERT(boxes.size() % (2*d+1) == 0,"Size of box vector is incorrect");
    index_t nboxes = boxes.size() / (2*d+1);
    gsMatrix<real_t> result(d,2*nboxes);
    std::vector<index_t> levels(nboxes);
    index_t lvl;
    index_t offset = 0;
    gsVector<real_t> low(d), upp(d);

    gsInfo<<"ancestor k = 1: ";
    for (index_t k=0; k!=5; k++)
        gsInfo<<boxes.at(k)<<"\t";
    gsInfo<<"\n";

    for (index_t k = 0; k!=nboxes; k++)
    {
        offset = k*(2*d+1);
        lvl = boxes[offset];
        levels[k] = lvl;
        for (index_t i=0; i!=d; i++)
        {
            const gsKnotVector<real_t> & kv = basis->tensorLevel(lvl).knots(i);
            low[i] = kv.uValue(boxes[1+i+offset]);
            upp[i] = kv.uValue(boxes[1+i+d+offset]);
        }
        result.col(2*k)   = low;
        result.col(2*k+1) = upp;
    }
    return std::make_pair(result,levels);
}

// Based on gsHTensorBasis.asElements
template<int d>
std::vector<index_t> getIndices(const gsHTensorBasis<d,real_t> * basis, const gsMatrix<real_t> & boxes, const std::vector<index_t> & levels)
{
    GISMO_ASSERT(boxes.cols() / 2 == levels.size(), "Number of boxes should be equal to the number of levels");

    // Each box will be represented by 2*d+1 entries specifying
    // <level to be refined to>,<lower corner>,<upper corner>
    const int offset = 2*d+1;

    // Initialize vector of size
    // "entries per box" times "number of boxes":
    std::vector<index_t> refVector( offset * levels.size() );
    gsMatrix<real_t> ctr(d,1);

    // Loop over all boxes:
    for(index_t i = 0; i < levels.size(); i++)
    {
        ctr = ( boxes.col( 2*i ) + boxes.col( 2*i+1) )/2;

        for(index_t j = 0; j < boxes.rows();j++)
        {
            // Convert the parameter coordinates to (unique) knot indices
            const gsKnotVector<real_t> & kv = basis->tensorLevel(levels[i]).knots(j);
            int k1 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),
                                       boxes(j,2*i  ) ) - 1).uIndex();
            int k2 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,
                                       boxes(j,2*i+1) ) - 1).uIndex();

            // Trivial boxes trigger some refinement
            if ( k1 == k2)
            {
                if (0!=k1) {--k1;}
                ++k2;
            }

            // Store the data...
            refVector[i*offset]       = levels.at(i);
            refVector[i*offset+1+j]   = k1;
            refVector[i*offset+1+j+d] = k2;
        }
    }
    // gsDebug<<"begin\n";
    // for (std::vector<unsigned>::const_iterator i = refVector.begin(); i != refVector.end(); ++i)
    //     std::cout << *i << ' ';
    // gsDebug<<"end\n";

    return refVector;
}

template<int d>
std::vector<index_t> limitBoxesToLevel(const gsHTensorBasis<d,real_t> * basis, const std::vector<index_t> & boxes, index_t lvl)
{
    GISMO_ASSERT(boxes.size() % (2*d+1) == 0,"Size of box vector is incorrect");
    index_t nboxes = boxes.size() / (2*d+1);
    std::vector<index_t> result;
    gsMatrix<real_t> corners;

    std::tie(corners,std::ignore) = getCoordinates<d>(basis,boxes);

    index_t offset1 = 0;
    gsInfo<<"boxes: \n";
    for (index_t b=0; b!=boxes.size() / (2*d+1); b++)
    {
        for (index_t k=0; k!=2*d+1; k++)
            gsInfo<<boxes.at(offset1 + k)<<"\t";
        gsInfo<<"\n";
        offset1 += 2*d+1;
    }

    gsVector<real_t> cen(d);
    index_t offset = 0;
    std::vector<bool> keep(nboxes);
    for (index_t k=0; k!=nboxes; k++)
    {
        cen = (corners.col(2*k) + corners.col(2*k+1))/2;
        offset = k*(2*d+1);
        gsInfo<<"Handling box ";
        for (index_t k=0; k!=2*d+1; k++)
            gsInfo<<boxes.at(offset + k)<<"\t";
        gsInfo<<" ...";

        gsInfo<<" target level = "<<lvl<<" current level "<<basis->getLevelAtPoint(cen)<<"\n";

        if (basis->getLevelAtPoint(cen) != lvl)
            continue;

        for (index_t i=0; i!=2*d+1; i++)
            result.push_back(boxes.at(offset + i));
    }
    return result;
}


template<int d>
std::vector<index_t> get_parent_of_cell(const gsHDomainIterator<real_t,d> * domHIt)
{
    std::vector<index_t> parent(2*d+1);
    gsVector<index_t> idx;
    index_t lvl = domHIt->getLevel();

    // get_parent_of_cell
    idx = domHIt->elementMultiIndex();

    parent[0] = lvl-1;
    for (index_t i = 0; i!=d; ++i)
    {
        parent[i+1]    = idx[i]/2;
        parent[i+1+d]    = idx[i]/2+1;
    }

    return parent;
}

template<int d>
std::vector<index_t> get_parent_of_cell(const std::vector<index_t> boxes)
{
    std::vector<index_t> parents(boxes.size());
    GISMO_ASSERT(boxes.size() % (2*d+1) == 0,"Size of box vector is incorrect");
    index_t nboxes = boxes.size() / (2*d+1);
    index_t offset = 0;
    for (index_t k = 0; k!=nboxes; k++)
    {
        parents[offset] = boxes[offset]-1;
        for (index_t i = 0; i!=d; i++)
        {
            parents[1 + i + offset]     = boxes[1 + i + offset] / 2;
            parents[1 + i + d + offset] = boxes[1 + i + d + offset] / 2 +  (boxes[1 + i + d + offset] % 2 != 0);
        }
        offset += 2*d+1;
    }
    return parents;
}


template<int d>
std::vector<index_t> get_ancestor_of_cell(const gsHDomainIterator<real_t,d> * domHIt, const index_t k)
{
    index_t lvl = domHIt->getLevel();
    GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k > 0,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);

    std::vector<index_t> ancestors;

    ancestors = get_parent_of_cell<d>(domHIt);
    while (k < lvl-1)
    {
        ancestors = get_parent_of_cell<d>(ancestors);
        lvl--;
    }


    return ancestors;
}

template<int d>
std::vector<index_t> get_ancestor_of_cell(const std::vector<index_t> boxes, const index_t k)
{
    std::vector<index_t> parents(boxes.size());
    GISMO_ASSERT(boxes.size() == (2*d+1),"Size of box vector is incorrect, should be "<<2*d+1<<" (single box)");
    std::vector<index_t> ancestors;

    GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k > 0,"Current level should be larger than 0!");

    index_t lvl = boxes[0];
    ancestors = get_parent_of_cell<d>(boxes);
    while (k < lvl-1)
    {
        ancestors = get_parent_of_cell<d>(ancestors);
        lvl--;
    }
    return ancestors;
}

/**
 * @brief      For a given element with center \a cen and level \a lvl, compute the support extension
 *
 * For a given element \f$Q\f$ with center \a cen and level \a lvl, compute the support extension, defined as
 *
 * \f$ \tilde{Q} = \left\{ Q'   \in G^l : \exists\beta^l\in\mathcal{B}^l, \text{supp } \beta^l \cap Q' \neq \empty \wedge \text{supp }\beta^l\cap Q \neq \empty \right\} \f$
 *
 *  With \f$G^l\f$ the grid at level \f$l\f$, \f$\mathcal{B}^l\f$ the basis at level \f$l\f$.
 *  In other words, the elements that have basis functions with support in elements \f$Q'\f$ and \f$Q\f$
 *
 * @param[in]  basis  The basis
 * @param[in]  cen    The center of the element
 * @param[in]  lvl    The level of the element
 *
 * @tparam     d      Domain dimension
 *
 * @return     The support extension in knot-value format
 */
template<int d>
gsMatrix<real_t> _get_support_extension(const gsHTensorBasis<d,real_t> * basis, const gsVector<real_t> & cen, index_t lvl)
{
    // Compute actives
    gsMatrix<index_t> acts;
    basis->tensorLevel(lvl).active_into(cen,acts);

    // Support extension
    gsMatrix<real_t> boxes(d,2*acts.rows());
    for (index_t act = 0; act!=acts.rows(); act++)
        boxes.block(0,2*act,d,2) = basis->tensorLevel(lvl).support(acts(act,0));

    return boxes;
}

/**
 * @brief      For a given element with center \a cen and level \a lvl, compute the support extension
 *
 * For a given element \f$Q\f$ with center \a cen and level \a lvl, compute the support extension, defined as
 *
 * \f$ \tilde{Q} = \left\{ Q'   \in G^l : \exists\beta^l\in\mathcal{B}^l, \text{supp } \beta^l \cap Q' \neq \empty \wedge \text{supp }\beta^l\cap Q \neq \empty \right\} \f$
 *
 *  With \f$G^l\f$ the grid at level \f$l\f$, \f$\mathcal{B}^l\f$ the basis at level \f$l\f$.
 *  In other words, the elements that have basis functions with support in elements \f$Q'\f$ and \f$Q\f$
 *
 * @param[in]  basis  The basis
 * @param[in]  cen    The center of the element
 * @param[in]  lvl    The level of the element
 *
 * @tparam     d      Domain dimension
 *
 * @return     The support extension in box format
 */
template<int d>
std::vector<index_t> get_support_extension(const gsHTensorBasis<d,real_t> * basis, const gsVector<real_t> & cen, index_t lvl)
{
    gsMatrix<real_t> boxes = _get_support_extension<d>(basis,cen,lvl);
    std::vector<index_t> lvls(boxes.cols() / 2);
    std::fill(lvls.begin(),lvls.end(),lvl);
    std::vector<index_t> Boxes = getIndices<d>(basis,boxes,lvls);

    index_t offset = 0;
    gsInfo<<"support extension for box with center: "<<cen.transpose()<<" on level "<<lvl<<"\n";
    for (index_t b=0; b!=Boxes.size() / (2*d+1); b++)
    {
        for (index_t k=0; k!=2*d+1; k++)
            gsInfo<<Boxes.at(offset + k)<<"\t";
        gsInfo<<"\n";
        offset += 2*d+1;
    }
    return Boxes;
}

template<int d>
std::vector<index_t> get_support_extension(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt)
{
    // Get center of the element
    const gsVector<real_t>& cen = domHIt->centerPoint();

    index_t lvl = domHIt->getLevel();

    return get_support_extension<d>(basis,cen,lvl);
}

// Is this one correct?
template<int d>
std::vector<index_t> get_support_extension(const gsHTensorBasis<d,real_t> * basis, const std::vector<index_t> boxes)
{
    // GISMO_NO_IMPLEMENTATION;

    // Obtain the center of each box in knot-coordinates
    gsMatrix<real_t> coords;
    std::vector<index_t> levels;
    std::tie(coords,levels) = getCoordinates<2>(basis,boxes);

    gsMatrix<real_t> centers(coords.rows(),coords.cols()/2);
    for (index_t k=0; k!=centers.cols(); k++)
        centers.col(k) = ( coords.col(2*k) + coords.col(2*k+1) ) / 2;

    std::vector<index_t> suppext, result;
    for (index_t k=0; k!=coords.cols() / 2; k++)
    {
        suppext = get_support_extension<d>(basis,centers.col(k),levels.at(k));
        result.insert(result.end(), suppext.begin(), suppext.end());
    }
    return result;
}

template<int d>
std::vector<index_t> get_multilevel_support_extension(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t k)
{
    index_t lvl = domHIt->getLevel();
    if (k==lvl)
    {
        return get_support_extension<d>(basis,domHIt);
    }
    else
    {
        std::vector<index_t> ancestors = get_ancestor_of_cell<d>(domHIt,k);
        return get_support_extension<d>(basis,ancestors);
    }
}

template<int d>
std::vector<index_t> get_multilevel_support_extension(const gsHTensorBasis<d,real_t> * basis, const std::vector<index_t> boxes, const index_t k)
{
    GISMO_ASSERT(boxes.size() == (2*d+1),"Size of box vector is incorrect, should be "<<2*d+1<<" (single box)");
    index_t lvl = boxes[0];
    if (k==lvl)
    {
        return get_support_extension<d>(basis,boxes);
    }
    else
    {
        std::vector<index_t> ancestors = get_ancestor_of_cell<d>(boxes,k);
        return get_support_extension<d>(basis,ancestors);
    }
}

template<int d>
std::vector<index_t> get_H_neighborhood(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t m)
{
    std::vector<index_t> neighborhood, extension;
    index_t lvl = domHIt->getLevel();
    index_t k = lvl - m - 1;
    if (k < 0)
        return neighborhood;
    else
    {
        // Get multi level support extension on level k
        extension = get_multilevel_support_extension<d>(basis,domHIt,k);
        // Eliminate elements which are too low
        neighborhood = limitBoxesToLevel<d>(basis,extension,k);
    }
    return neighborhood;
}

template<int d>
std::vector<index_t> get_T_neighborhood(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t m)
{
    std::vector<index_t> neighborhood, extension;
    index_t lvl = domHIt->getLevel();
    index_t k = lvl - m - 1;
    if (k - 1< 0)
        return neighborhood;
    else
    {
        // Get multi level support extension on level k
        extension = get_multilevel_support_extension<d>(basis,domHIt,k);
        // Get the parents
        extension = get_parent_of_cell<d>(extension);
        // Eliminate elements which are too low
        neighborhood = limitBoxesToLevel<d>(basis,extension,k-1);
    }
    return neighborhood;
}

template<int d>
std::vector<index_t> mark_recursive()
{

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
    for (index_t k=0; k!=mpBspline.nPatches(); ++k)
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

    for (index_t k=0; k!=5; k++)
        gsDebugVar(boxes[k]);

        mp.patch(0).refineElements(boxes);



    // std::vector<index_t> boxes(5);
    // boxes[0] = 2;
    // boxes[1] = 0;
    // boxes[2] = 0;
    // boxes[3] = 4;
    // boxes[4] = 4;

    // mp.patch(0).refineElements(boxes);

    // boxes[0] = 3;
    // boxes[1] = 0;
    // boxes[2] = 3;
    // boxes[3] = 3;
    // boxes[4] = 6;

    // mp.patch(0).refineElements(boxes);

    // gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));

    // gsDebugVar(basis->maxLevel());

    // boxes[0] = 1;
    // boxes[1] = 1;
    // boxes[2] = 1;
    // boxes[3] = 2;
    // boxes[4] = 2;

    // mp.patch(0).unrefineElements(boxes);

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));


    gsDebugVar(basis->maxLevel());
    return 0;

    // get_support_extension
    typename gsBasis<real_t>::domainIter domIt = basis->makeDomainIterator();
    gsHDomainIterator<real_t,2> * domHIt = nullptr;
    domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
    GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

    gsVector<index_t> idx;
    index_t lvl;
    std::vector<index_t> parent;
    std::vector<index_t> suppext;
    for (; domHIt->good(); domHIt->next() )
    {
        // support extension
        suppext = get_support_extension<2>(basis,domHIt);

        // original_cell
        gsVector<index_t> idx = domHIt->elementMultiIndex();
        parent.resize(5);
        parent[0] = domHIt->getLevel();
        parent[1] = idx[0];
        parent[2] = idx[1];
        parent[3] = idx[0]+1;
        parent[4] = idx[1]+1;


        gsInfo<<"lvl\tbox\n";
        for (index_t k = domHIt->getLevel(); k >= 0; k--)
        {
            gsInfo<<parent.at(0)<<"\t\t";
            for (index_t k=1; k!=5; k++)
                gsInfo<<parent.at(k)<<"\t";
            gsInfo<<"\n";
            parent = get_parent_of_cell<2>(parent);
        }


        // std::vector<index_t> tmp;
        // // get_parent_of_cell
        // tmp = get_parent_of_cell<2>(domHIt);
        // gsInfo<<"box: ";
        // for (index_t k=0; k!=5; k++)
        //     gsInfo<<tmp.at(k)<<"\t";
        // gsInfo<<"\n";

        // // get_parent_of_cell
        // parent = get_parent_of_cell<2>(parent);
        // gsInfo<<"box: ";
        // for (index_t k=0; k!=5; k++)
        //     gsInfo<<parent.at(k)<<"\t";
        // gsInfo<<"\n";

        // parent = get_parent_of_cell<2>(parent);
        // gsInfo<<"box: ";
        // for (index_t k=0; k!=5; k++)
        //     gsInfo<<parent.at(k)<<"\t";
        // gsInfo<<"\n";

        // parent = get_parent_of_cell<2>(parent);
        // gsInfo<<"box: ";
        // for (index_t k=0; k!=5; k++)
        //     gsInfo<<parent.at(k)<<"\t";
        // gsInfo<<"\n";


        std::vector<index_t> ancestor = get_ancestor_of_cell<2>(domHIt,1);
        gsInfo<<"ancestor k = 1: ";
        for (index_t k=0; k!=5; k++)
            gsInfo<<ancestor.at(k)<<"\t";
        gsInfo<<"\n";

        suppext = get_multilevel_support_extension<2>(basis,domHIt,1);
        // gsInfo<<"supp ext k = 1: ";
        // for (index_t k=0; k!=5; k++)
        //     gsInfo<<suppext.at(k)<<"\t";
        // gsInfo<<"\n";
        //
        suppext = get_H_neighborhood<2>(basis,domHIt,1);
        suppext = get_T_neighborhood<2>(basis,domHIt,1);

    }

    mp.patch(0).unrefineElements(parent);



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

