/** @file gsHFitting.h

    @brief Data-driven fitting using hierarchical splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#pragma once

#include <gsHSplines/gsHFitting.h>

namespace gismo {

void print_video_message(){
    gsInfo<<"----------------------------------------------------------------\n";
    gsInfo<<"Hello, so far everything is fine.\n";
    gsInfo<<"----------------------------------------------------------------\n";
}

template <class T>
std::vector<T> vector_linspace(T a,
                               T b,
                               index_t N){
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}



template<class T>
gsMatrix<T> myPointGrid( gsVector<T> const & a,
                         gsVector<T> const & b,
                         gsVector<index_t> const & np )
{
    gsInfo<<"grid generation of dimension:\n";
    gsInfo<< a.size() << "x" << np.prod() <<"\n";

    gsInfo<<"Lower corener:"<<"\n";
    gsInfo<<a<<"\n";
    gsInfo<<"Upper corener:"<<"\n";
    gsInfo<<b<<"\n";

    gsInfo<<"np:"<<"\n";
    gsInfo<<np<<"\n";


    gsInfo<<"Can we generate the following matrix?\n";
    gsMatrix<T> pg(a.size(), np.prod());
    gsInfo<<"done.\n";

    gsInfo<<"grid iterator generation:\n";
    gsGridIterator<T,CUBE> pt(a, b, np);

    for(index_t c = 0; pt; ++pt, ++c)
        pg.col(c) = *pt;
    return pg;
}

template <short_t d, class T>
std::vector<index_t> gsHFitting<d, T>::getAllBoxes()
{
    // cells contains lower corners of elements marked for refinment from maxLevel
    std::vector<index_t> cells;

    // boxes contains elements marked for refinement from differnet levels,
    // format: { level lower-corners  upper-corners ... }
    std::vector<index_t> boxes;

    for (size_t index = 0; index != this->m_param_values.cols(); index++)
    {
        appendBox(boxes, cells, this->m_param_values.col(index));
    }

    return boxes;
}

/*
template <short_t d, class T>
std::vector<gsVector<T>> gsHFitting<d, T>::getAllCells(){
    std::vector<gsVector<T>> cells;
    //gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    typename gsTHBSplineBasis<d, T>::domainIter domIt = basis->makeDomainIterator();
        for (; domIt->good(); domIt->next() )    // loop over all elements
        {
            cells.push_back(domIt->lowerCorner()); // gsVector
            cells.push_back(domIt->upperCorner()); // geVector
        }
    return cells;
}
*/

template <class T>
gsMatrix<T> gsVectorASgsMatrix(const gsVector<T>& col0,
                               const gsVector<T>& col1)
{
    gsMatrix<T> elMatrix(2,2);
    elMatrix.col(0) << col0;
    elMatrix.col(1) << col1;
    return elMatrix;
}


/// returns the cells covering the whole domain.
template <short_t d, class T>
std::vector<gsMatrix<T>> gsHFitting<d, T>::getAllCells(){
    std::vector<gsMatrix<T>> cells;
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    typename gsTHBSplineBasis<d, T>::domainIter domIt = basis->makeDomainIterator();
        for (; domIt->good(); domIt->next() )    // loop over all elements
        {
            gsMatrix<T> elMatrix(2,2);
            elMatrix.col(0)<< domIt->lowerCorner(); // first column  = lower corner
            elMatrix.col(1)<< domIt->upperCorner(); // second column = upper corner
            cells.push_back(elMatrix);
        }
    return cells;
}

/// returns the cells covering the whole domain.
template <short_t d, class T>
std::vector<gsMatrix<T>> getAllCellsFromBasis(gsTHBSplineBasis<d, T>& basis){
    std::vector<gsMatrix<T>> cells;
    //gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    typename gsTHBSplineBasis<d, T>::domainIter domIt = basis.makeDomainIterator();
        for (; domIt->good(); domIt->next() )    // loop over all elements
        {
            gsMatrix<T> elMatrix(2,2);
            elMatrix.col(0)<< domIt->lowerCorner(); // first column  = lower corner
            elMatrix.col(1)<< domIt->upperCorner(); // second column = upper corner
            cells.push_back(elMatrix);
        }
    return cells;
}

template <class T>
bool is_point_inside_cell(const gsMatrix<T>& parameter,
                          const gsMatrix<T>& element)
{
    const real_t x = parameter(0, 0);
    const real_t y = parameter(1, 0);

    return element(0, 0) <= x && x < element(0, 1) &&
           element(1, 0) <= y && y < element(1, 1);
}

template <class T>
bool is_point_inside_cell(const T x,
                          const T y,
                          const gsMatrix<T>& element)
{
    bool condition = (element(0, 0) <= x && x <= element(0, 1) && element(1, 0) <= y && y <= element(1, 1));
    return condition;
}


// difference with is_point_inside_cell in the inclusion of the left and right interval extremes.
template <class T>
bool is_point_within_cell(const gsMatrix<T>& parameter,
                          const gsMatrix<T>& element)
{
    const real_t x = parameter(0, 0);
    const real_t y = parameter(1, 0);

    return element(0, 0) < x && x < element(0, 1) &&
           element(1, 0) < y && y < element(1, 1);
}

template <class T>
bool is_point_within_cell(const T x,
                          const T y,
                          const gsMatrix<T>& element)
{
    bool condition = (element(0, 0) < x && x < element(0, 1) && element(1, 0) < y && y < element(1, 1));
    return condition;
}

template <short_t d, class T>
std::vector<T> gsHFitting<d, T>::cellAvgErr(const std::vector<T>& errors,
                                            const gsMatrix<T>& parameters){
    std::vector<T> cellAvgErrs;
    std::vector<gsMatrix<T>> elements = this->getAllCells();

    for(index_t el = 0; el < elements.size(); el++){
        gsMatrix<T> a_cell = elements[el];
        std::vector<T> a_cellErrAvg;
        for(index_t it=0; it < parameters.cols(); it++){
            const T xx = parameters.col(it)(0);
            const T yy = parameters.col(it)(1);
                if (is_point_inside_cell(xx, yy, a_cell))
                {
                    a_cellErrAvg.push_back(errors[it]);
                }
            }
            T sum_of_elems = 0;
            for(typename std::vector<T>::iterator errIt = a_cellErrAvg.begin(); errIt != a_cellErrAvg.end(); ++errIt)
            {
               sum_of_elems += *errIt;
            }

            if (a_cellErrAvg.size() > 0){
                T tempErr = sum_of_elems/a_cellErrAvg.size();
                cellAvgErrs.push_back(tempErr);
            }
            else{
                cellAvgErrs.push_back(0);
            }
            a_cellErrAvg.clear();
    }

    return cellAvgErrs;
}

template <short_t d, class T>
std::vector<T> gsHFitting<d, T>::cellMaxErr(const std::vector<T>& errors,
                                            const gsMatrix<T>& parameters){
    std::vector<T> cellMaxErrs;
    std::vector<gsMatrix<T>> elements = this->getAllCells();

    for(index_t el = 0; el < elements.size(); el++){
        gsMatrix<T> a_cell = elements[el];
        std::vector<T> a_cellErrs;
        for(index_t it=0; it < parameters.cols(); it++){
            const T xx = parameters.col(it)(0);
            const T yy = parameters.col(it)(1);
                if (is_point_inside_cell(xx, yy, a_cell))
                {
                    a_cellErrs.push_back(errors[it]);
                }
            }
            T cell_max_err = 0;
            for(typename std::vector<T>::iterator errIt = a_cellErrs.begin(); errIt != a_cellErrs.end(); ++errIt)
            {
                if (*errIt > cell_max_err){
                    cell_max_err = *errIt;
                }
            }

            if (a_cellErrs.size() > 0){
                cellMaxErrs.push_back(cell_max_err);
            }
            else{
                cellMaxErrs.push_back(0);
            }
            a_cellErrs.clear();
    }

    return cellMaxErrs;
}









template<short_t d, class T>
void gsHFitting<d, T>::computeApproximation(index_t maxPcIter)
{
    // We run one fitting step and compute the errors
    this->compute(m_lambda);
    //parameter correction
    this->parameterCorrection(1e-7, maxPcIter, 1e-4);//closestPoint accuracy, orthogonality tolerance
    //Approximation error
    this->computeErrors();
}

template <short_t d, class T>
std::vector<T> gsHFitting<d, T>::BoxAvgErr(const std::vector<T>& errors,
                                            const gsMatrix<T>& parameters){
    std::vector<T> boxAvgErr;
    return boxAvgErr;
}

template <short_t d, class T>
std::vector<T> gsHFitting<d, T>::BoxMaxErr(const std::vector<T>& errors,
                                            const gsMatrix<T>& parameters){
    std::vector<T> boxMaxErr;
    return boxMaxErr;
}




template <short_t d, class T>
void gsHFitting<d, T>::appendBox(std::vector<gsMatrix<T>>& cells,
                                 std::vector<index_t>& boxes)
                                 // Last argument-> the one which is changed by the function
{
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    const int maxLvl = basis->maxLevel();
    const tensorBasis & tBasis = *(basis->getBases()[maxLvl]);

    // loop over the "already selected" cells
    for(auto el=cells.begin(); el!=cells.end(); el++){
        gsMatrix<T> matrix = *el;
        gsVector<index_t, d> a_cell;
        for (index_t dim = 0; dim != d; dim++)
        {
            const gsKnotVector<T> & kv = tBasis.component(dim).knots();
            a_cell(dim) = kv.uFind(matrix(dim,0)).uIndex();
        }

        //gsInfo<<"cell:\n"<< *el <<"\n";
        //gsInfo<<"a_cell:\n"<< a_cell <<"\n";

        // get level of a cell
        gsVector<index_t, d> a_cell_upp = a_cell + gsVector<index_t, d>::Ones();
        const int cell_lvl = basis->tree().query3(a_cell, a_cell_upp, maxLvl) + 1;

        // get the box
        gsVector<index_t> box(2 * d + 1);
        box[0] = cell_lvl;
        for (short_t dim = 0; dim != d; dim++)
        {
            const unsigned numBreaks = basis->numBreaks(cell_lvl, dim) - 1 ;

            unsigned lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const unsigned shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const unsigned shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            // apply extensions
            index_t low = ( (lowIndex > m_ext[dim]) ? (lowIndex - m_ext[dim]) : 0 );
            index_t upp = ( (lowIndex + m_ext[dim] + 1 < numBreaks) ?
                             (lowIndex + m_ext[dim] + 1) : numBreaks );

            box[1 + dim    ] = low;
            box[1 + d + dim] = upp;
        }

        append(boxes, box);
    }
}







} // namespace gismo