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

typedef typename std::vector<index_t>                           CellType;
typedef typename std::vector<std::vector<index_t>>              CellContainerType;


template<short_t d>
void printCells(const CellContainerType & cells)
{
    index_t offset = 0;
    gsInfo<<cells.size()<<" cells: \n";
    for (size_t b=0; b!=cells.size(); b++)
    {
        if (cells[b].size()==0)
            gsInfo<<"empty!\n";
        else
        {
            for (size_t k=0; k!=cells[b].size(); k++)
                gsInfo<<cells[b].at(k)<<"\t";
            gsInfo<<"\n";
            offset += cells[b].size();
        }
    }
}

// template<short_t d>
// CellContainerType removeEmptyCells(CellContainerType & cells)
// {
//     CellContainerType newcells;
//     for (index_t k=0; k!=cells.size(); k++)
//     {
//         if (cells[k].size()==0)
//             continue;

//         GISMO_ASSERT(cells[k].size() == (2*d+1),"Size of cell vector is incorrect, must be "<<2*d+1<<" but "<<cells[k].size()<<" is provided.");
//         newcells.push_back(cells[k]);
//     }
//     return newcells;
// }


// CALL STD SORT
template<short_t d>
CellContainerType cellUnion(CellContainerType & region1, CellContainerType & region2)
{
    CellContainerType result;

    gsInfo<<"region1:\n";
    printCells<d>(region1);
    gsInfo<<"region2:\n";
    printCells<d>(region2);
    if (region1.size()!=0 && region2.size()!=0)
    {
        auto comp = [](CellType & a, CellType & b)
                        {
                            gsInfo<<"a = :";
                            for (index_t k =0; k!=a.size(); k++)
                                gsInfo<<a.at(k)<<"\t";
                            gsInfo<<"\n";

                            gsInfo<<"b = :";
                            for (index_t k =0; k!=b.size(); k++)
                                gsInfo<<b.at(k)<<"\t";
                            gsInfo<<"\n";

                            GISMO_ASSERT(a.size()==b.size(),"Cells must have the same size. a.size()="<<a.size()<<", b.size()="<<b.size());
                            bool same = true;
                            for (index_t k=0; k!=a.size(); k++)
                                same &= a[k]==b[k];
                            return same;
                        };

        std::sort(region1.begin(),region1.end(),comp);
        std::sort(region2.begin(),region2.end(),comp);

        /// SORT FIRST
        std::set_union(region1.begin(),region1.end(),region2.begin(),region2.end(),result.begin(),comp
                       );
    }
    else if (region1.size()!=0 && region2.size()==0)
        result.insert(result.end(),region1.begin(),region1.end());
    else if (region1.size()==0 && region2.size()!=0)
        result.insert(result.end(),region2.begin(),region2.end());
    else
        // Do nothing

    return result;
}

// // Strict weak ordering
// // Should be contained but not same!
// template<short_t d>
// CellContainerType removeDuplicates(CellContainerType & region)
// {

//     CellContainerType result;

//     gsInfo<<"region:\n";
//     printCells<d>(region);
//     std::sort(region.begin(),region.end(),
//                     [](CellType & a, CellType & b)
//                     {
//                         return a.contains(b) && !(a.isSame(b));
//                     }
//                    );

//     bool !finished = false;
//     auto it = region.begin();
//     auto it1 = it+1;
//     CellContainerType cells;
//     while(!finished && it1!=region.end())
//     {
//         if (it1->contains(*it+1))
//             it++;
//         else
//             cells.push_back(*it);
//     }

//     std::sort(region.begin(),region.end(),
//                     [](CellType & a, CellType & b)
//                     {
//                         return a.contains(b);
//                     }
//                    );

//     return result;
// }

template<short_t d>
std::pair<gsMatrix<real_t>,index_t> getCoordinates(const gsHTensorBasis<d,real_t> * basis, const CellType cell)
{
    gsMatrix<real_t> result(d,2);
    index_t level;
    gsVector<real_t> low(d), upp(d);

    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, must be "<<2*d+1<<" but "<<cell.size()<<" is provided.");
    level = cell[0];
    for (index_t i=0; i!=d; i++)
    {
        const gsKnotVector<real_t> & kv = basis->tensorLevel(level).knots(i);
        low[i] = kv.uValue(cell[1+i]);
        upp[i] = kv.uValue(cell[1+i+d]);
    }
    result.col(0) = low;
    result.col(1) = upp;
    return std::make_pair(result,level);
}

template<short_t d>
std::pair<gsMatrix<real_t>,std::vector<index_t>> getCoordinates(const gsHTensorBasis<d,real_t> * basis, const CellContainerType cells)
{
    gsInfo<<"cells for coordinates:\n";
    printCells<d>(cells);
    index_t ncells = cells.size();
    gsMatrix<real_t> result(d,2*ncells);
    std::vector<index_t> levels(ncells);
    gsVector<real_t> low(d), upp(d);

    for (index_t k = 0; k!=ncells; k++)
    {
        if (cells[k].size()==0)
        {
            result.col(2*k).setZero();
            result.col(2*k+1).setZero();
            levels[k] = -1;
            continue;
        }

        GISMO_ASSERT(cells[k].size() == (2*d+1),"Size of cell vector is incorrect, must be "<<2*d+1<<" but "<<cells[k].size()<<" is provided.");
        levels[k] = cells[k][0];
        for (index_t i=0; i!=d; i++)
        {
            const gsKnotVector<real_t> & kv = basis->tensorLevel(levels[k]).knots(i);
            low[i] = kv.uValue(cells[k][1+i]);
            upp[i] = kv.uValue(cells[k][1+i+d]);
        }
        result.col(2*k)   = low;
        result.col(2*k+1) = upp;
    }
    return std::make_pair(result,levels);
}

// Based on gsHTensorBasis.asElements
template<short_t d>
CellContainerType getIndices(const gsHTensorBasis<d,real_t> * basis, const gsMatrix<real_t> & cells, const std::vector<index_t> & levels)
{
    GISMO_ASSERT(cells.cols() / 2 == static_cast<index_t>(levels.size()), "Number of cells should be equal to the number of levels");

    // Each cell will be represented by 2*d+1 entries specifying
    // <level to be refined to>,<lower corner>,<upper corner>

    // Initialize vector of size
    // "entries per cell" times "number of cells":
    CellContainerType container( levels.size() );
    gsMatrix<real_t> ctr(d,1);

    // Loop over all cells:
    for(size_t i = 0; i < levels.size(); i++)
    {
        ctr = ( cells.col( 2*i ) + cells.col( 2*i+1) )/2;

        container[i].resize(2*d+1);
        container[i][0] = levels.at(i);

        for(index_t j = 0; j < cells.rows();j++)
        {
            // Convert the parameter coordinates to (unique) knot indices
            const gsKnotVector<real_t> & kv = basis->tensorLevel(levels[i]).knots(j);
            int k1 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),
                                       cells(j,2*i  ) ) - 1).uIndex();
            int k2 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,
                                       cells(j,2*i+1) ) - 1).uIndex();

            // Trivial cells trigger some refinement
            if ( k1 == k2)
            {
                if (0!=k1) {--k1;}
                ++k2;
            }

            // Store the data...
            container[i][1+j]   = k1;
            container[i][1+j+d] = k2;
        }
    }
    // gsDebug<<"begin\n";
    // for (std::vector<unsigned>::const_iterator i = refVector.begin(); i != refVector.end(); ++i)
    //     std::cout << *i << ' ';
    // gsDebug<<"end\n";

    return container;
}

// Based on gsHTensorBasis.asElements
template<short_t d>
CellType getIndices(const gsHTensorBasis<d,real_t> * basis, const gsMatrix<real_t> & cell, const index_t level)
{
    GISMO_ASSERT(cell.cols()==2,"A single cell must be defined by a d x 2 matrix, but received "<<cell.rows()<<" x "<<cell.cols());
    std::vector<index_t> levels = {level};
    CellContainerType cells = getIndices<d>(basis,cell,levels);
    return cells.at(0);
}

/**
 * @brief      Get the cell indices of a cell represented by a H-domain iterator
 *
 * @param[in]  basis   The H-basis
 * @param[in]  domHIt  The H-domain iterator
 *
 * @tparam     d       domain dimension
 *
 * @return     The cell in index coordinates
 */
template<short_t d>
CellType getIndices(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt)
{
    gsMatrix<real_t> cell(d,2);
    std::vector<index_t> levels(1);
    levels[0] =  domHIt->getLevel();
    cell.col(0) =domHIt->lowerCorner() ;
    cell.col(1) =domHIt->upperCorner() ;

    CellContainerType cells = getIndices<d>(basis,cell,levels);
    return cells[0];
}

// template<short_t d>
// CellType getIndices(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<d,real_t> * domHIt)
// {
//     gsVector<index_t> idx;
//     idx = domHIt->elementMultiIndex();

//     CellType cell(2*d+1);
//     cell[0] = domHIt->getLevel()
//     for (index_t i = 0; i != d; i++)
//     {
//         cell[i+1]    = idx[i];
//         cell[i+1+d]  = idx[i]+1;
//     }
//     return cell;
// }

template<short_t d>
CellContainerType getActivesOnLevel(const gsHTensorBasis<d,real_t> * basis, const CellContainerType & cells, index_t lvl)
{
    CellContainerType actives;
    CellType active;
    gsMatrix<real_t> corners;
    gsVector<real_t> cen;

    for (index_t c=0; c!=cells.size(); c++)
    {
        if (cells[c].size()==0)
            continue;


        std::tie(corners,std::ignore) = getCoordinates<d>(basis,cells[c]);
        cen = (corners.col(0) + corners.col(1))/2;
        if (basis->getLevelAtPoint(cen) != lvl)
            continue;

        for (index_t i=0; i!=2*d+1; i++)
            active.push_back(cells[c].at(i));

        actives.push_back(active);
    }
    return actives;
}

// Works for multiple elements
template<short_t d>
CellContainerType limitBoxesToLevel(const gsHTensorBasis<d,real_t> * basis, const CellContainerType & cells, index_t lvl)
{
    index_t ncells = cells.size();
    CellContainerType result(ncells);

    gsMatrix<real_t> corners;

    gsInfo<<"limitBoxesToLevel:\n";
    printCells<d>(cells);

    std::tie(corners,std::ignore) = getCoordinates<d>(basis,cells);

    printCells<d>(cells);

    gsVector<real_t> cen(d);
    std::vector<bool> keep(ncells);
    for (index_t b=0; b!=ncells; b++)
    {
        GISMO_ASSERT(cells[b].size() == (2*d+1),"Size of cell vector is incorrect, must be "<<2*d+1<<" but "<<cells[b].size()<<" is provided.");

        cen = (corners.col(2*b) + corners.col(2*b+1))/2;
        gsInfo<<"Handling cell ";
        for (size_t k=0; k!=cells[b].size(); k++)
            gsInfo<<cells[b][k]<<"\t";
        gsInfo<<" ...";

        gsInfo<<" target level = "<<lvl<<" current level "<<basis->getLevelAtPoint(cen)<<"\n";

        if (basis->getLevelAtPoint(cen) != lvl)
            continue;

        for (index_t i=0; i!=2*d+1; i++)
            result[b].push_back(cells[b].at(i));
    }
    return result;
}

// Works for a single element
template<short_t d>
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

// Works for multiple elements
template<short_t d>
CellType get_parent_of_cell(const CellType cell)
{
    CellType parent(cell.size());
    GISMO_ASSERT(cell.size() % (2*d+1) == 0,"Size of cell vector is incorrect, must be "<<2*d+1<<" but "<<cell.size()<<" is provided.");
    parent[0] = cell[0]-1;
    for (index_t i = 0; i!=d; i++)
    {
        parent[1 + i]     = cell[1 + i] / 2;
        parent[1 + i + d] = cell[1 + i + d] / 2 +  (cell[1 + i + d] % 2 != 0);
    }
    return parent;
}

// Works for multiple elements
template<short_t d>
CellContainerType get_parent_of_cells(const CellContainerType cells)
{
    CellContainerType parents(cells.size());
    for (index_t k = 0; k!=cells.size(); k++)
    {
        GISMO_ASSERT(cells[k].size() == (2*d+1),"Size of cell "<<k<<" is incorrect, should be "<<2*d+1);
        parents[k] = cells[k];
        parents[k][0] = cells[k][0]-1;
        for (index_t i = 0; i!=d; i++)
        {
            parents[k][1 + i]     = cells[k][1 + i] / 2;
            parents[k][1 + i + d] = cells[k][1 + i + d] / 2 +  (cells[k][1 + i + d] % 2 != 0);
        }
    }
    return parents;
}


// Works for a single element
template<short_t d>
CellType get_ancestor_of_cell(const gsHDomainIterator<real_t,d> * domHIt, const index_t k)
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

// Works for a single element
template<short_t d>
CellType get_ancestor_of_cell(const CellType cell, const index_t k)
{
    CellType parent(2*d+1);
    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, should be "<<2*d+1<<" (single cell)");
    CellType ancestor;

    index_t lvl = cell[0];
    GISMO_ASSERT(lvl > k,"Current level should be larger than requested level, l = "<<lvl<<", k = "<<k);
    GISMO_ASSERT(k >= 0,"Current level should be larger than or equal to 0!");

    ancestor = get_parent_of_cell<d>(cell);
    while (k < lvl-1)
    {
        ancestor = get_parent_of_cell<d>(ancestor);
        lvl--;
    }
    return ancestor;
}

/**
 * @brief      For a single given element with center \a cen and level \a lvl, compute the support extension
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
template<short_t d>
gsMatrix<real_t> _get_support_extension_of_cell(const gsHTensorBasis<d,real_t> * basis, const gsVector<real_t> & cen, index_t lvl)
{
    // Compute actives
    gsMatrix<index_t> acts;
    basis->tensorLevel(lvl).active_into(cen,acts);

    // Support extension
    gsMatrix<real_t> cells(d,2*acts.rows());
    for (index_t act = 0; act!=acts.rows(); act++)
        cells.block(0,2*act,d,2) = basis->tensorLevel(lvl).support(acts(act,0));

    return cells;
}

/**
 * @brief      For a single given element with center \a cen and level \a lvl, compute the support extension
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
 * @return     The support extension in index format
 */
template<short_t d>
CellContainerType get_support_extension_of_cell(const gsHTensorBasis<d,real_t> * basis, const gsVector<real_t> & cen, index_t lvl)
{
    gsMatrix<real_t> cells = _get_support_extension_of_cell<d>(basis,cen,lvl);
    std::vector<index_t> lvls(cells.cols() / 2);
    std::fill(lvls.begin(),lvls.end(),lvl);
    CellContainerType Boxes = getIndices<d>(basis,cells,lvls);

    // gsInfo<<"support extension for cell with center: "<<cen.transpose()<<" on level "<<lvl<<"\n";
    // printCells<d>(Boxes);
    return Boxes;
}

// Based on a single cell (given by the domHIt)
template<short_t d>
CellContainerType get_support_extension_of_cell(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt)
{
    // Get center of the element
    const gsVector<real_t>& cen = domHIt->centerPoint();

    index_t lvl = domHIt->getLevel();
    return get_support_extension_of_cell<d>(basis,cen,lvl);
}

// Is this one correct?
// Works for multiple cells
template<short_t d>
CellContainerType get_support_extension_of_cell(const gsHTensorBasis<d,real_t> * basis, const CellType cell)
{
    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, should be a single cell");
    // Obtain the center of each cell in knot-coordinates
    gsMatrix<real_t> coords;
    index_t level;
    std::tie(coords,level) = getCoordinates<2>(basis,cell);

    gsMatrix<real_t> centers(coords.rows(),coords.cols()/2);
    GISMO_ASSERT(coords.cols()/2 == 1,"Only a single coordinate should be returned");

    centers.col(0) = ( coords.col(0) + coords.col(1) ) / 2;

    return get_support_extension_of_cell<d>(basis,centers.col(0),level);
}

template<short_t d>
CellContainerType get_multilevel_support_extension(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t k)
{
    index_t lvl = domHIt->getLevel();
    GISMO_ASSERT(k <= lvl,"Ancestor must be requested from a level lower than the current level (k < l), but k="<<k<<" and lvl="<<lvl);
    if (k==lvl)
    {
        return get_support_extension_of_cell<d>(basis,domHIt);
    }
    else
    {
        CellType ancestors = get_ancestor_of_cell<d>(domHIt,k);
        return get_support_extension_of_cell<d>(basis,ancestors);
    }
}

template<short_t d>
CellContainerType get_multilevel_support_extension(const gsHTensorBasis<d,real_t> * basis, const CellType cell, const index_t k)
{
    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, should be "<<2*d+1<<" (single cell)");
    index_t lvl = cell[0];
    GISMO_ASSERT(k <= lvl,"Ancestor must be requested from a level lower than the current level (k < l), but k="<<k<<" and lvl="<<lvl);
    if (k==lvl)
    {
        return get_support_extension_of_cell<d>(basis,cell);
    }
    else
    {
        CellType ancestor = get_ancestor_of_cell<d>(cell,k);
        return get_support_extension_of_cell<d>(basis,ancestor);
    }
}

template<short_t d>
CellContainerType get_H_neighborhood(const gsHTensorBasis<d,real_t> * basis, const CellType cell, const index_t m)
{
    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, should be "<<2*d+1<<" (single cell)");
    CellContainerType neighborhood, extension;
    index_t lvl = cell[0];
    index_t k = lvl - m + 1;
    if (k < 0)
        return neighborhood;
    else
    {
        // Get multi level support extension on level k
        extension = get_multilevel_support_extension<d>(basis,cell,k);
        // Eliminate elements which are too low
        neighborhood = getActivesOnLevel<d>(basis,extension,k);
    }
    return neighborhood;
}


template<short_t d>
CellContainerType get_H_neighborhood(const gsHTensorBasis<d,real_t> * basis, const CellContainerType cells, const index_t m)
{
    CellContainerType neighborhood, result;
    result = cells;
    for (index_t c=0; c!=cells.size(); c++)
    {
        neighborhood = get_H_neighborhood(basis,cells[c],m);
        if (neighborhood.size()!=0)
            neighborhood = cellUnion<d>(neighborhood,result);
    }
    gsInfo<<"H-neigh\n";
    printCells<d>(neighborhood);
    return neighborhood;
}

template<short_t d>
CellContainerType get_H_neighborhood(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t m)
{
    CellType cell = getIndices<d>(basis,domHIt);
    return get_H_neighborhood(basis,cell,m);
}

template<short_t d>
CellContainerType get_T_neighborhood(const gsHTensorBasis<d,real_t> * basis, const CellType cell, const index_t m)
{
    GISMO_ASSERT(cell.size() == (2*d+1),"Size of cell vector is incorrect, should be "<<2*d+1<<" (single cell)");
    CellContainerType neighborhood, extension;
    index_t lvl = cell[0];
    index_t k = lvl - m + 2;
    if (k - 1 < 0)
        return neighborhood;
    else
    {
        // Get multi level support extension on level k
        extension = get_multilevel_support_extension<d>(basis,cell,k);
        gsInfo<<"T neighborhood MLS extension:\n";
        printCells<d>(extension);
        // Get the parents
        extension = get_parent_of_cells<d>(extension);
        gsInfo<<"T neighborhood extension:\n";
        printCells<d>(extension);
        // Eliminate elements which are too low
        neighborhood = limitBoxesToLevel<d>(basis,extension,k-1);
    }
    return neighborhood;
}

template<short_t d>
CellContainerType get_T_neighborhood(const gsHTensorBasis<d,real_t> * basis, const gsHDomainIterator<real_t,d> * domHIt, const index_t m)
{
    CellType cell = getIndices<d>(basis,domHIt);
    return get_T_neighborhood(basis,cell,m);
}

template<short_t d>
CellContainerType mark_recursive(const gsHTensorBasis<d,real_t> * basis, const CellContainerType & marked, index_t l, index_t m)
{
    // make the set marked_l
    CellContainerType result;
    CellContainerType marked_l, marked_k, neighbors;
    marked_l = getActivesOnLevel<2>(basis,marked,l);

        gsDebugVar(l);
    gsInfo<<"marked_l:\n";
    printCells<d>(marked_l);
    // if T-admissibility
    //
    //

    index_t k;
    neighbors = get_H_neighborhood<d>(basis,marked_l,m);
    printCells<d>(marked_l);
    if (neighbors.size()!=0)
    {
        k = l-m+1;
        gsDebugVar(k);
        printCells<d>(marked);
        marked_k = getActivesOnLevel<2>(basis,marked,k);
        printCells<d>(marked_k);
        marked_k.insert(marked_k.end(),neighbors.begin(),neighbors.end());
        // marked_k = cellUnion<d>(marked_k,neighbors)
        result = mark_recursive<d>(basis,marked_k,k,m);
    }
    return result;

}

// template<d>
// CellContainerType admissible_refinement(const gsHTensorBasis<d,real_t> * basis, const CellContainerType cells, index_t m)
// {
//     // if T-admissibility
//     //
//     //

//     CellContainerType marked(cells.size());
//     for (index_t k=0; k!=cells.size(); k++)
//         marked[k] = {cells[k]};

//     for (index_t l = 0; l!=basis->maxLevel(); l++)
//     {
//         marked = mark_recursive<d>(basis,marked,m);
//     }
//     CellContainerType neighbors = get_H_neighborhood<d>(basis,cells,m);

// }

// template<short_t d>
// std::vector<index_t> mark_recursive()
// {

// }

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

    gsMatrix<> cell;
    CellType selectedIdx = getIndices<2>(basis,domHIt);
    CellContainerType temp(1);
    temp[0] = selectedIdx;
    printCells<2>(temp);
    std::tie(cell,std::ignore) = getCoordinates<2>(basis,selectedIdx);

    writeCells("selected.csv",cell);

    gsHBoxContainer<2,real_t> tmp;
    gsHBox<2,real_t> box0(selectedIdx);
    gsDebugVar(box0);
    gsHBox<2,real_t> box1(domHIt);
    gsDebugVar(box1);
    gsDebugVar(box1.getCoordinates());
    // gsDebugVar(box1.getParent());
    // gsDebugVar(box1.getParent().getParent());
    gsDebugVar(box1.getAncestor(1));
    gsHBox<2,real_t> box2(box0.lowerIndex(),box0.upperIndex(),box0.level());
    gsDebugVar(box2);
    gsHBox<2,real_t> box3(box0.lowerIndex(),box0.upperIndex(),box0.level(),basis);
    gsDebugVar(box3);

    gsDebugVar(box0.contains(box1));

    tmp = gsHBoxContainer<2,real_t>(box1.getSupportExtension();
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
    tmp = gsHBoxContainer<2,real_t>(box1.getContainer().markHrecursive(3,2));
    gsDebugVar(tmp);


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

    gsDebugVar(HBContainer3.getParents());




    CellContainerType marked(1);
    marked[0] = selectedIdx;
    gsInfo<<"Marked:\n";
    printCells<2>(marked);

    CellContainerType region, neighbors, result;

    index_t m =2;
    for (index_t l=0; l!=basis->maxLevel()+1; l++)
    {
        result = mark_recursive(basis,marked,l,m);
        marked.insert(marked.end(),result.begin(),result.end());

        gsInfo<<"Level "<<l<<"\n";
        printCells<2>(marked);
    }




return 0;

    CellContainerType H_cells = get_H_neighborhood<2>(basis,domHIt,1);
    gsMatrix<> Hneighborhood;
    std::tie(Hneighborhood,std::ignore) = getCoordinates<2>(basis,H_cells);
    writeCells("Hneighborhood.csv",Hneighborhood);

    CellContainerType T_cells = get_T_neighborhood<2>(basis,domHIt,1);
    gsMatrix<> Tneighborhood;
    std::tie(Tneighborhood,std::ignore) = getCoordinates<2>(basis,T_cells);
    writeCells("Tneighborhood.csv",Tneighborhood);



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

