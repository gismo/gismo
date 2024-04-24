/** @file gsPreCICE.h

    @brief Header file for using gsPreCICE extension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-...)
*/

#pragma once

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo {

template <short_t DIM, class T, template<short_t _DIM, class _T> class basis_type>
inline void getKnots(const gsBasis<T> & source, std::vector<gsKnotVector<T>> & tensorKnots)
{
    if ( const basis_type<DIM,T> * basis = dynamic_cast<const basis_type<DIM,T>*>(&source) )
        for (index_t d=0; d!=DIM; d++)
            tensorKnots[d] = basis->knots(d);
}

/**
 * @brief      { function_description }
 *
 *  Make a matrix with the knot vectors
 *  [[x1_1, x1_2, ..., nan , nan ,..., nan , nan , ...]
 *   [nan , nan , ..., x2_1, x2_2,..., nan , nan , ...]
 *   [nan , nan , ..., nan , nan ,..., x2_1, x2_2, ...]]
 *
 * @param[in]  basis  The basis
 *
 * @tparam     T      { description }
 *
 * @return     { description_of_the_return_value }
 */
template<class T>
gsVector<T> knotsToVector(const gsBasis<T> & basis)
{
    const size_t DIM = 1;
    std::vector<gsKnotVector<T>> tensorKnots(DIM);
    switch (DIM)
    {
        case 1:
            getKnots<1,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<1,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        default:
            GISMO_ERROR("Basis type not understood");
            break;
    }
    gsDebugVar(tensorKnots[0]);
    std::vector<size_t> sizes(DIM);
    sizes[0] = tensorKnots[0].size();

    gsVector<T> knots(std::accumulate(sizes.begin(),sizes.end(),0));

    for (index_t i = 0; i < knots.size(); ++i)
        knots[i] = tensorKnots[0][i];
    return knots;
}

template<class T>
gsMatrix<T> knotsToMatrix(const gsBasis<T> & basis)
{
    const size_t DIM = basis.domainDim();
    std::vector<gsKnotVector<T>> tensorKnots(DIM);
    switch (DIM)
    {
        case 1:
            getKnots<1,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<1,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        case 2:
            getKnots<2,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<2,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        case 3:
            getKnots<3,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<3,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        default:
            GISMO_ERROR("Basis type not understood");
            break;
    }

    std::vector<size_t> sizes(DIM);
    std::vector<size_t> strides(DIM);
    for (size_t d=0; d!=DIM; d++)
        sizes[d] = tensorKnots[d].size();

    strides[0]=0;
    for (size_t d=1; d!=DIM; d++)
        strides[d]= strides[d-1]+sizes[d-1];

    gsMatrix<T> knots(DIM,std::accumulate(sizes.begin(),sizes.end(),0));
    knots.setConstant(std::nan("1"));
    for (size_t d=0; d!=DIM; d++)
        knots.block(d,strides[d],1,sizes[d]) = tensorKnots[d].asMatrix();

    return knots;
}

template<class T>
gsMatrix<T> knotVectorUnpack(const gsMatrix<T> & knots, index_t numBoundaries)
{
    gsMatrix<> kv_unpacked;
    kv_unpacked = knots.row(0);
    kv_unpacked.resize(knots.cols()/numBoundaries,numBoundaries);

    return kv_unpacked;
}

template<class T>
gsMatrix<T> unPackControlPoints(const gsMatrix<T> & controlPoints, const gsMatrix<T> & kv_unpacked, index_t knot_index, index_t cp_index)
{
    // number of cps n = N_knot - p - 1
    int counter = cp_index;
    gsVector<> temp = kv_unpacked.row(knot_index);

    // Calculate the amount of control points based on the knot vector
    for (index_t i = 0; i < temp.size(); ++i)
    {
        if(temp[i] == 0)
            counter ++;
    }
    counter++;

    gsMatrix<T> unpackedCps(controlPoints.rows(), counter - cp_index);


    // Simplified copying of control points (just a placeholder)
    for (index_t i = 0; i < controlPoints.rows(); ++i) 
    {
        index_t diff = counter - cp_index;
        index_t startIndex = cp_index - diff + 1;  // Calculate the start index
        for (size_t j = 0; j < diff; ++j) 
        {
            unpackedCps(i,j) = controlPoints(i,j + startIndex);
        }
    }
    return unpackedCps;
}



template<class T>
typename gsBasis<T>::Ptr knotMatrixToBasis(const gsMatrix<T> & knots)
{
    gsBasis<> * basis;
    const short_t DIM = knots.rows();

    std::vector<gsKnotVector<T>> KVs(DIM);
    index_t k=0;
    for (size_t d=0; d!=DIM; d++)
    {
        std::vector<T> tmp;
        std::copy_if(knots.row(d).begin(),knots.row(d).end(),
                     std::back_inserter(tmp),
                     [](T a){return !math::isnan(a);});
        KVs[d] = gsKnotVector<T>(tmp);
    }

    switch(DIM)
    {
        case 1:
            basis = new  gsBSplineBasis<T>(KVs[0]);
            break;
        case 2:
            basis = new  gsTensorBSplineBasis<2,T>(KVs);
            break;
        case 3:
            basis = new  gsTensorBSplineBasis<3,T>(KVs);
            break;
    }
    return memory::make_shared_not_owned(basis);
}

} //namespace gismo
