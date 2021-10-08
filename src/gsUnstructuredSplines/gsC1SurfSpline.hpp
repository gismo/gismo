/** @file gsC1SurfSpline.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include<gsUnstructuredSplines/gsC1SurfEdge.h>
#include<gsUnstructuredSplines/gsC1SurfVertex.h>

namespace gismo
{

template<short_t d,class T>
void gsC1SurfSpline<d,T>::defaultOptions()
{

}
template<short_t d,class T>
void gsC1SurfSpline<d,T>::init()
{

}   


template<short_t d,class T>
void gsC1SurfSpline<d,T>::compute()
{
    // Compute Inner Basis functions
    index_t shift_row = 0, shift_col = 0;
    for(size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        index_t dim_u = m_bases[np].getInnerBasis().component(0).size();
        index_t dim_v = m_bases[np].getInnerBasis().component(1).size();

        index_t row_i = 0;
        for (index_t j = 2; j < dim_v-2; ++j)
            for (index_t i = 2; i < dim_u-2; ++i)
            {
                // TODO
                m_matrix.insert(shift_row + row_i,shift_col + j*dim_u+i) = 1.0;
                ++row_i;
            }

        shift_row += m_bases[np].size_rows();
        shift_col += m_bases[np].size_cols();
    }

    // Compute Interface Basis functions
    /*
     *  (Side-1) * 2 + 0/1 = Index
     *  0 == for lower vertex index, 1 == higher vertex index
     *
     *  Side 1, Vertex 1 == 0
     *  Side 1, Vertex 3 == 1
     *  Side 2, Vertex 2 == 2
     *  Side 2, Vertex 4 == 3
     *  Side 3, Vertex 1 == 4
     *  ...
     */
    std::vector<std::vector<gsMultiPatch<T>>> vertex_bf(m_patches.nPatches(), std::vector<gsMultiPatch<T>>(8));
    for (size_t numInt = 0; numInt < m_patches.interfaces().size(); numInt++)
    {
        const boundaryInterface & item = m_patches.interfaces()[numInt];

        gsC1SurfEdge c1SurfEdge(m_patches, item);
        //TODO
        //gsApproxC1Edge<d, T> approxC1Edge(m_patches, m_bases, item, numInt, m_options);
        //approxC1Edge.saveBasisInterface(m_matrix);
    }
    // Compute Edge Basis functions
    for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
    {
        const patchSide & bit = m_patches.boundaries()[numBdy];

        gsC1SurfEdge c1SurfEdge(m_patches, bit);

        //TODO
        //gsApproxC1Edge<d, T> approxC1Edge(m_patches, m_bases, bit, numBdy, m_options);
        //approxC1Edge.saveBasisBoundary(m_matrix);
    }
    // Compute Vertex Basis functions
    for (size_t numVer = 0; numVer < m_patches.vertices().size(); numVer++)
    {
        std::vector<patchCorner> allcornerLists = m_patches.vertices()[numVer];
        std::vector<size_t> patchIndex;
        std::vector<size_t> vertIndex;
        for (size_t j = 0; j < allcornerLists.size(); j++)
        {
            patchIndex.push_back(allcornerLists[j].patch);
            vertIndex.push_back(allcornerLists[j].m_index);
        }

        gsC1SurfVertex c1SurVertex(m_patches, patchIndex, vertIndex);
        c1SurVertex.computeG1InternalVertexBasis(m_options);
        //TODO
        //gsApproxC1Vertex<d, T> approxC1Vertex(m_patches, m_bases, patchIndex, vertIndex, numVer, m_options);
        //approxC1Vertex.saveBasisVertex(m_matrix);

    }

    m_matrix.makeCompressed();

    if (m_options.getSwitch("info"))
    {
        gsInfo << "Dim for Patches: \n";
        for(size_t np = 0; np < m_patches.nPatches(); ++np)
        {
            gsInfo << "(" << m_bases[np].size_rows() << "," << m_bases[np].size_cols() << "), ";
        }
        gsInfo << "\n";
    }
}

template<short_t d,class T>
void gsC1SurfSpline<d,T>::writeParaviewSinglePatch(int patchID, std::string type)
{

}

template<short_t d,class T>
void gsC1SurfSpline<d,T>::plotParaview(std::string fn, int npts) 
{

}



template<short_t d,class T>
void gsC1SurfSpline<d,T>::createPlusMinusSpace(gsKnotVector<T> & kv1, gsKnotVector<T> & kv2,
                                               gsKnotVector<T> & kv1_patch, gsKnotVector<T> & kv2_patch,
                                               gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result)
{
    std::vector<real_t> knots_unique_1 = kv1.unique();
    std::vector<real_t> knots_unique_2 = kv2.unique();

    std::vector<index_t> knots_mult_1 = kv1.multiplicities();
    std::vector<index_t> knots_mult_2 = kv2.multiplicities();

    std::vector<real_t> patch_kv_unique_1 = kv1_patch.unique();
    std::vector<index_t> patch_kv_mult_1 = kv1_patch.multiplicities();

    std::vector<real_t> patch_kv_unique_2 = kv2_patch.unique();
    std::vector<index_t> patch_kv_mult_2 = kv2_patch.multiplicities();

    std::vector<real_t> knot_vector_plus, knot_vector_minus;

    if (knots_unique_1 != knots_unique_2)
        gsInfo << "NOT IMPLEMENTED YET 1: Plus, Minus space \n";

    if (kv1.degree() != kv2.degree())
        gsInfo << "NOT IMPLEMENTED YET 2: Plus, Minus space \n";

    if (knots_mult_1 != knots_mult_2)
        gsInfo << "NOT IMPLEMENTED YET 4: Plus, Minus space \n";

    // todo: fix for general regularity
    index_t m, p;
    p = math::max(kv1.degree(), kv2.degree());
    m = kv1.multiplicityIndex(p+1);

    kv1_result = kv2; // == kv2
    if (m != 1)
        kv1_result.reduceMultiplicity(1);

    kv2_result = kv2; // == kv2
    kv2_result.degreeDecrease(1);
    if (m != 1)
        kv2_result.reduceMultiplicity(1);

/*
* TODO Add geometry inner knot regularity
*
index_t i_3 = 0, i_4 = 0;

std::vector<real_t>::iterator it3 = patch_kv_unique_1.begin();
std::vector<real_t>::iterator it4 = patch_kv_unique_2.begin();

std::vector<real_t>::iterator it2 = knots_unique_2.begin();
for(std::vector<real_t>::iterator it = knots_unique_1.begin(); it != knots_unique_1.end(); ++it)
{
if (*it == *it2)
{
    knot_vector_plus.push_back(*it);
    knot_vector_minus.push_back(*it);
    ++it2;
}
else if (*it < *it2)
{
    //knot_vector_plus.push_back(*it);
    //knot_vector_minus.push_back(*it);
}
else if (*it > *it2)
{
    while (*it > *it2)
    {
        //knot_vector_plus.push_back(*it2);
        //knot_vector_minus.push_back(*it2);
        ++it2;
    }
    knot_vector_plus.push_back(*it2);
    knot_vector_minus.push_back(*it2);
    ++it2;
}
}

// Repeat the first and the last vector p or p-1 times
kv1_result = gsKnotVector<>(knot_vector_plus);
kv1_result.degreeIncrease(p);
if (kv1.multiplicities()[1] > p-2 && knots_unique_1[1] != 1) // TODO Check
kv1_result.increaseMultiplicity(kv1.multiplicities()[1]-2);
kv2_result = gsKnotVector<>(knot_vector_minus);
kv2_result.degreeIncrease(p-1);
if (kv2.multiplicities()[1] > p-2 && knots_unique_2[1] != 1) // TODO Check
kv2_result.increaseMultiplicity(kv2.multiplicities()[1]-2);
*/
}


template<short_t d,class T>
void gsC1SurfSpline<d,T>::createPlusMinusSpace(gsKnotVector<T> & kv1,
                                               gsKnotVector<T> & kv1_patch,
                                               gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result)
{
    std::vector<real_t> knots_unique_1 = kv1.unique();

    std::vector<real_t> patch_kv_unique_1 = kv1_patch.unique();
    std::vector<index_t> patch_kv_mult_1 = kv1_patch.multiplicities();

    index_t m,p;
    p = math::max(kv1.degree(), 0);
    m = kv1.multiplicityIndex(p+1);

    kv1_result = kv1;
    if (m != 1)
        kv1_result.reduceMultiplicity(1);

    kv2_result = kv1;
    kv2_result.degreeDecrease(1);
    if (m != 1)
        kv2_result.reduceMultiplicity(1);


    /*
    * TODO Add geometry inner knot regularity
    *
    index_t i_3 = 0, i_4 = 0;

    std::vector<real_t>::iterator it3 = patch_kv_unique_1.begin();
    std::vector<real_t>::iterator it4 = patch_kv_unique_2.begin();

    std::vector<real_t> knot_vector_plus, knot_vector_minus;

    for(std::vector<real_t>::iterator it = knots_unique_1.begin(); it != knots_unique_1.end(); ++it)
    {
        knot_vector_plus.push_back(*it);
        knot_vector_minus.push_back(*it);
    }

    // Repeat the first and the last vector p or p-1 times
    kv1_result = gsKnotVector<>(knot_vector_plus);
    kv1_result.degreeIncrease(p);
    kv2_result = gsKnotVector<>(knot_vector_minus);
    kv2_result.degreeIncrease(p-1);
     */
}

template<short_t d,class T>
void gsC1SurfSpline<d,T>::createGluingDataSpace(gsKnotVector<T> & kv1, gsKnotVector<T> & kv2,
                                                gsKnotVector<T> & kv1_patch, gsKnotVector<T> & kv2_patch,
                                                gsKnotVector<T> & kv_result)
{
    //index_t p_tilde = math::max(math::max(kv1.degree(), kv2.degree())-1,3); // max(p-1,3)
    //index_t r_tilde = math::max(p_tilde - 2, 1); // C^2 gluing data

    std::vector<real_t> knots_unique_1 = kv1.unique();
    std::vector<real_t> knots_unique_2 = kv2.unique();

    std::vector<real_t> knot_vector;

    /*
     * TODO Add geometry inner knot regularity
     */
    if (knots_unique_1 != knots_unique_2)
        gsInfo << "\n\n ERROR: Interfaces are not matching!!! \n\n";

    knot_vector = knots_unique_2; // = knots_unique_1
/*
std::vector<real_t>::iterator it2 = knots_unique_2.begin();
for(std::vector<real_t>::iterator it = knots_unique_1.begin(); it != knots_unique_1.end(); ++it)
{
if (*it == *it2)
{
    knot_vector.push_back(*it);
    ++it2;
}
else if (*it < *it2)
{
    //knot_vector.push_back(*it);
}
else if (*it > *it2)
{
    while (*it > *it2)
    {
        //knot_vector.push_back(*it2);
        ++it2;
    }
    knot_vector.push_back(*it2);
    ++it2;
}
}
*/

    kv_result = gsKnotVector<>(knot_vector);
    kv_result.degreeIncrease(p_tilde);
    kv_result.increaseMultiplicity(p_tilde-r_tilde-1);
} // createGluingDataSpace


template<short_t d,class T>
void gsC1SurfSpline<d,T>::createLocalEdgeSpace(gsKnotVector<T> & kv_plus, gsKnotVector<T> & kv_minus,
                                               gsKnotVector<T> & kv_gD_1, gsKnotVector<T> & kv_gD_2,
                                               gsKnotVector<T> & kv_patch_1, gsKnotVector<T> & kv_patch_2,
                                               gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result)
{
    index_t p_1 = math::max(kv_plus.degree()+kv_gD_1.degree()-1, kv_minus.degree()+kv_gD_1.degree() );
    //index_t p_2 = math::max(kv_plus.degree()+kv_gD_2.degree()-1, kv_minus.degree()+kv_gD_2.degree() ); == p_1

    std::vector<real_t> knots_unique_plus = kv_plus.unique(); // == kv_minus.unique()

    if (knots_unique_plus != kv_minus.unique())
        gsInfo << "ERROR LOKAL EDGE SPACE \n";

    kv1_result = gsKnotVector<>(knots_unique_plus);
    kv1_result.degreeIncrease(p_1);
    if (knots_unique_plus[1] != 1)
    {
        index_t r_plus = kv_plus.degree() - kv_plus.multiplicities()[1]; // The same for all
        index_t r_minus = kv_minus.degree() - kv_minus.multiplicities()[1]; // The same for all
        index_t r_tilde = kv_gD_1.degree() - kv_gD_1.multiplicities()[1]; // The same for all
        //gsInfo << "R_tilde " << r_tilde << "\n";
        //gsInfo << "r_plus " << r_plus << "\n";
        //gsInfo << "r_minus " << r_minus << "\n";

        index_t r = math::min(r_tilde, math::min(r_plus, r_minus));
        //gsInfo << "r " << r << "\n";

        kv1_result.increaseMultiplicity(p_1-r-1);
    }
    // ==
    kv2_result = kv1_result;

/*
index_t p_plus_diff = p_1 - kv_plus.degree();
index_t p_gD_diff = p_1 - kv_gD_1.degree();
index_t p_patch_diff = p_1 - kv_1.degree();

std::vector<real_t> knots_unique_plus = kv_plus.unique();
std::vector<real_t> knots_unique_gD = kv_gD_1.unique();

std::vector<real_t> knots_unique_1 = kv_1.unique();

knots_unique_1.erase(knots_unique_1.begin()); // First
knots_unique_1.pop_back(); // Last

std::vector<index_t> patch_kv_mult_plus = kv_plus.multiplicities();
std::vector<index_t> patch_kv_mult_gD = kv_gD_1.multiplicities();

std::vector<index_t> patch_kv_mult_1 = kv_1.multiplicities();

if (knots_unique_plus != knots_unique_gD)
gsInfo << "\n\nERROR: TODO \n\n";

std::vector<real_t> knot_vector;

index_t i_plus = 0;
index_t i_1 = 1;
std::vector<real_t>::iterator it_1 = knots_unique_1.begin();
for(std::vector<real_t>::iterator it = knots_unique_plus.begin(); it != knots_unique_plus.end(); ++it, ++i_plus)
{
if (*it_1 == *it && it_1 != knots_unique_1.end())
{
    index_t i_temp = 0;
    while(i_temp < math::max(patch_kv_mult_1[i_1]+p_patch_diff, math::max(patch_kv_mult_plus[i_plus]+p_plus_diff, patch_kv_mult_gD[i_plus]+p_gD_diff)))
    {
        knot_vector.push_back(*it);
        ++i_temp;
    }

    ++it_1;
    ++i_1;
}
else
{
    index_t i_temp = 0;
    while(i_temp < math::max(patch_kv_mult_plus[i_plus]+p_plus_diff, patch_kv_mult_gD[i_plus]+p_gD_diff))
    {
        knot_vector.push_back(*it);
        ++i_temp;
    }
}



}


kv1_result = gsKnotVector<>(knot_vector);
// ==
kv2_result = gsKnotVector<>(knot_vector);
*/
} // createLocalEdgeSpace

template<short_t d,class T>
void gsC1SurfSpline<d,T>::createLocalEdgeSpace(gsKnotVector<T> & kv_plus, gsKnotVector<T> & kv_minus,
                                               gsKnotVector<T> & kv_patch_1,
                                               gsKnotVector<T> & kv1_result)
{
    index_t p_1 = math::max(kv_plus.degree(), kv_minus.degree() );

    std::vector<real_t> knots_unique_plus = kv_plus.unique(); // == kv_minus.unique()

    if (knots_unique_plus != kv_minus.unique())
        gsInfo << "ERROR LOKAL EDGE SPACE \n";

    kv1_result = gsKnotVector<>(knots_unique_plus);
    kv1_result.degreeIncrease(p_1);
    if (knots_unique_plus[1] != 1)
    {
        index_t r_plus = kv_plus.degree() - kv_plus.multiplicities()[1]; // The same for all
        index_t r_minus = kv_minus.degree() - kv_minus.multiplicities()[1]; // The same for all

        index_t r = math::min(r_plus, r_minus);

        kv1_result.increaseMultiplicity(p_1-r-1);
    }
} // createLocalEdgeSpace

template<short_t d,class T>
void gsC1SurfSpline<d,T>::createLocalVertexSpace(gsTensorBSplineBasis<d, T> & basis_vertex, gsTensorBSplineBasis<d, T> & basis_vertex_result)
{
    index_t p_1 = basis_vertex.degree(0); // == basis_vertex.degree(1)
    //index_t p_tilde = math::max(p_1 - 1, 3); // TODO more general


    // todo: fix for general regularity
    index_t r;
    r = p_1 - basis_vertex.knots(0).multiplicityIndex(p_1);

    if (basis_vertex.degree(0) != basis_vertex.degree(1))
        gsInfo << "ERROR LOKAL Vertex SPACE \n";

    basis_vertex_result = basis_vertex;

    //gsInfo << "basis u " << basis_vertex_result.knots(0).asMatrix() << "\n";
    //gsInfo << "basis v " << basis_vertex_result.knots(1).asMatrix() << "\n";

    basis_vertex_result.degreeElevate(p_tilde, 0); // p + \tilde{p} - 1
    basis_vertex_result.degreeElevate(p_tilde, 1); // p + \tilde{p} - 1
    basis_vertex_result.reduceContinuity(r-1);



} // createLocalVertexSpace




} // namespace gismo
