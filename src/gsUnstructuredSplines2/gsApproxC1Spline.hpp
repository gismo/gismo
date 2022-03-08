/** @file gsApproxC1Spline.hpp

    @brief Construct the approx c1 spline space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include<gsUnstructuredSplines2/gsApproxC1Edge.h>
#include<gsUnstructuredSplines2/gsApproxC1Vertex.h>

#include <gsUnstructuredSplines2/gsApproxC1Utils.h>

namespace gismo
{

template<short_t d,class T>
void gsApproxC1Spline<d,T>::defaultOptions()
{
    /*
        to do: general
    */
    gsTensorBSplineBasis<d, T> basis = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(0));
    index_t p = basis.degree(0);

    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        gsTensorBSplineBasis<d, T> basis = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(np));
        if (p != basis.degree(0))
            gsWarn << "Not suitable for different degrees! \n";
    }

    m_options.addSwitch("info","Print debug information",  false );
    m_options.addSwitch("plot","Print debug information",  false );
    m_options.addSwitch("interpolation","Compute the basis with interpolation",  false );
    m_options.addSwitch("second","Compute the second biharmonic problem",  false );

    m_options.addInt("gluingDataDegree","Print debug information",  -1 );
    m_options.addInt("gluingDataSmoothness","Print debug information",  -1 );
}

template<short_t d,class T>
void gsApproxC1Spline<d,T>::init()
{
    index_t row_dofs = 0;

    m_bases.clear();
    m_bases.reserve(m_patches.nPatches()); // For each Patch
    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        // gsContainerBasisBase:
        // # basisContainer
        // - Interior space: [0] : inner,
        // - Edge spaces:    [1] : west, [2] : east, [3] : south, [4] : north,
        // - Vertex spaces:  [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
        // [9] For the initial space, maybe delete later
        gsContainerBasis<d,T> containerBasis(1); // for 9 subspaces
        m_bases.push_back(containerBasis);
    }

    // Create interior spline space
    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        gsTensorBSplineBasis<d, T> basis_inner = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(np));
        //m_bases[np].setBasis(9, basis_inner); // Inner

        // Construct special space for r = p - 1:
        // The first and the last knot (not 0,1) are repeated +1, e.g.
        // deg 3, r = 2: |||| || | [...] | || ||||
        /*index_t m,p;
        for (index_t uv = 0; uv < 2; uv++)
        {
            p = basis_inner.degree(uv);
            m = basis_inner.knots(uv).multiplicityIndex(p+1);
            if (m == 1)
            {
                T knot_u = basis_inner.knot(uv,basis_inner.degree(uv)+1);
                if (knot_u != 1)
                    basis_inner.insertKnot(knot_u,uv,1);

                if (knot_u != 0.5 && knot_u != 1)
                    basis_inner.insertKnot(1-knot_u,uv,1);
            }
        }*/
        index_t dim_u = basis_inner.component(0).size();
        index_t dim_v = basis_inner.component(1).size();

        row_dofs += (dim_u - 4) * (dim_v - 4);

        m_bases[np].setBasis(0, basis_inner); // Inner
    }

    // For loop over Interface to construct the spaces
    for (size_t numInt = 0; numInt < m_patches.interfaces().size(); numInt++)
    {
        const boundaryInterface & item = m_patches.interfaces()[numInt];

        const index_t dir_1 = item.first().side().index() > 2 ? 0 : 1;
        const index_t dir_2 = item.second().side().index() > 2 ? 0 : 1;

        // [!Basis space]
        gsTensorBSplineBasis<d, T> basis_11 = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(item.first().patch));
        gsTensorBSplineBasis<d, T> basis_22 = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(item.second().patch));
/*
        gsTensorBSplineBasis<d, T> basis;
        index_t dir;
        if (basis_11.component(dir_1).numElements() > basis_22.component(dir_2).numElements())
        {
            basis = basis_22;
            dir = dir_2;
        }
        else
        {
            basis = basis_11;
            dir = dir_1;
        }
        // [!Basis space]
*/
        // [!Plus Minus space]
        gsBSplineBasis<T> basis_plus;
        gsBSplineBasis<T> basis_minus;
        createPlusSpace(m_patches.patch(item.first().patch), basis_11, dir_1, basis_plus);
        createMinusSpace(m_patches.patch(item.first().patch), basis_11, dir_1, basis_minus);
        // [!Plus Minus space]

        // [!Gluing data space]
        gsBSplineBasis<T> basis_gluingData;
        createGluingDataSpace(m_patches.patch(item.first().patch), basis_11, dir_1, basis_gluingData,
                              m_options.getInt("gluingDataDegree"), m_options.getInt("gluingDataSmoothness"));
        // [!Gluing data space]

        // [!Edge space]
        gsTensorBSplineBasis<d, T> basis_edge_11, basis_edge_22;
        createEdgeSpace(m_patches.patch(item.first().patch), basis_11, dir_1, basis_plus, basis_minus, basis_gluingData, basis_edge_11);
        createEdgeSpace(m_patches.patch(item.second().patch), basis_22, dir_2, basis_plus, basis_minus, basis_gluingData, basis_edge_22);

        // [!Edge space]
        //m_bases[item.first().patch].setBasis(item.first().side().index(), basis_edge_11);
        //m_bases[item.second().patch].setBasis(item.second().side().index(), basis_edge_22);

        index_t numDofs = math::max(basis_plus.size() + basis_minus.size() - 10, 0);
        row_dofs += numDofs; // The same as side_2
    }

    // For loop over the Edge to construct the spaces
    for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
    {
        const patchSide & bit = m_patches.boundaries()[numBdy];

        index_t dir_1 = m_patches.boundaries()[numBdy].m_index < 3 ? 1 : 0;

        gsTensorBSplineBasis<d, T> basis = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(bit.patch));

        // [!Plus Minus space]
        gsBSplineBasis<T> basis_plus;
        gsBSplineBasis<T> basis_minus;
        createPlusSpace(m_patches.patch(bit.patch), basis, dir_1, basis_plus);
        createMinusSpace(m_patches.patch(bit.patch), basis, dir_1, basis_minus);
        // [!Plus Minus space]

        gsTensorBSplineBasis<d, T> basis_edge_11;
        createEdgeSpace(m_patches.patch(bit.patch), basis, dir_1, basis_plus, basis_minus, basis_edge_11);

        //m_bases[bit.patch].setBasis(bit.side().index(), basis_edge_11);

        index_t numDofs = math::max(basis_plus.size() + basis_minus.size() - 10, 0);
        row_dofs += numDofs;
    }

    // For loop over the Vertex to construct the spaces
    for (size_t numVer = 0; numVer < m_patches.vertices().size(); numVer++)
    {
        std::vector<patchCorner> allcornerLists = m_patches.vertices()[numVer];
        for (size_t j = 0; j < allcornerLists.size(); j++)
        {
            std::vector<patchSide> containingSides;
            allcornerLists[j].getContainingSides(d, containingSides);
            bool isInterface_1 = m_patches.isInterface(patchSide(allcornerLists[j].patch, containingSides.at(0).side()));
            bool isInterface_2 = m_patches.isInterface(patchSide(allcornerLists[j].patch, containingSides.at(1).side()));

            if (containingSides.at(0).side() < 3) // If isInterface_1 == v, then switch
            {
                bool isInterface_temp = isInterface_1;
                isInterface_1 = isInterface_2;
                isInterface_2 = isInterface_temp;
            }

            gsTensorBSplineBasis<d, T> basis_vertex_11;
            createVertexSpace(m_patches.patch(allcornerLists[j].patch), m_multiBasis.basis(allcornerLists[j].patch),
                              isInterface_1, isInterface_2, basis_vertex_11, m_options.getInt("gluingDataDegree"),
                              m_options.getInt("gluingDataSmoothness"));

            //m_bases[allcornerLists[j].patch].setBasis(allcornerLists[j].m_index + 4, basis_vertex_11);
        }
        row_dofs += 6;
    }

    // Initialise the matrix
    m_matrix.clear();
    index_t dim_col = 0;
    for (size_t i = 0; i < m_bases.size(); i++)
    {
        dim_col += m_bases[i].size();
    }

    m_matrix.resize(row_dofs, dim_col);
    const index_t nz = (m_multiBasis.basis(0).degree(0)*2 + 1)*row_dofs; // TODO
    m_matrix.reserve(nz);


    if (m_options.getSwitch("info"))
    {
        for (size_t i = 0; i < m_bases.size(); i++)
        {
            gsInfo << "Patch " << i << ": \n";
            for (index_t j = 0; j < m_bases[i].nPieces(); j++)
            {
                gsInfo << (j == 0 ? "Interior Space: " : ( j < 5 ? "Edge Space: " : "Vertex Space: ") ) << "\n";
                gsTensorBSplineBasis<d, T> basis = dynamic_cast<const gsTensorBSplineBasis<d, T>&>(m_bases[i].piece(j));
                std::vector<index_t> vec_1 = basis.knots(0).multiplicities();
                std::vector<index_t> vec_2 = basis.knots(1).multiplicities();
                gsAsVector<index_t> mult_1(vec_1);
                gsAsVector<index_t> mult_2(vec_2);
                gsInfo << mult_1.transpose() << "\n";
                gsInfo << mult_2.transpose() << "\n";
                gsInfo << "----------------------------------\n";
            }
        }
    }
}   


template<short_t d,class T>
void gsApproxC1Spline<d,T>::compute()
{
    // Compute Inner Basis functions
    index_t shift_row = 0, shift_col = 0;
    for(size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        index_t dim_u = m_bases[np].piece(0).component(0).size();
        index_t dim_v = m_bases[np].piece(0).component(1).size();

        index_t row_i = 0;
        for (index_t j = 2; j < dim_v-2; ++j)
            for (index_t i = 2; i < dim_u-2; ++i)
            {
                m_matrix.insert(shift_row + row_i, shift_col + j*dim_u+i) = 1.0;
                ++row_i;
            }

        shift_row += row_i;
        shift_col += m_bases[np].size();
    }

    // Interfaces
    for (size_t numInt = 0; numInt < m_patches.interfaces().size(); numInt++)
    {
        const boundaryInterface & item = m_patches.interfaces()[numInt];
        index_t side_1 = item.first().side().index();
        index_t side_2 = item.second().side().index();

        index_t patch_1 = item.first().patch;
        index_t patch_2 = item.second().patch;

        gsApproxC1Edge<d, T> approxC1Edge(m_patches, m_bases, item, numInt, m_options);
        std::vector<gsMultiPatch<T>> basisEdge = approxC1Edge.getEdgeBasis();

        index_t begin_col = 0, end_col = 0, shift_col = 0;
        for (index_t np = 0; np < patch_1; ++np)
            shift_col += m_bases[np].size();
//
//        for (index_t ns = 0; ns < side_1; ++ns)
//            begin_col += m_bases[patch_1].piece(ns).size();
//        for (index_t ns = 0; ns < side_1+1; ++ns)
//            end_col += m_bases[patch_1].piece(ns).size();
        end_col = m_bases[patch_1].piece(0).size();

        for (size_t ii = 0; ii < basisEdge[0].nPatches(); ++ii)
        {
            index_t jj = 0;
            for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                if (basisEdge[0].patch(ii).coef(jj, 0) * basisEdge[0].patch(ii).coef(jj, 0) > 1e-25)
                    m_matrix.insert(shift_row + ii, shift_col + j) = basisEdge[0].patch(ii).coef(jj, 0);
            }
        }

        begin_col = 0, end_col = 0, shift_col = 0;
        for (index_t np = 0; np < patch_2; ++np)
            shift_col += m_bases[np].size();

//        for (index_t ns = 0; ns < side_2; ++ns)
//            begin_col += m_bases[patch_2].piece(ns).size();
//        for (index_t ns = 0; ns < side_2+1; ++ns)
//            end_col += m_bases[patch_2].piece(ns).size();
        end_col = m_bases[patch_1].piece(0).size();

        for (size_t ii = 0; ii < basisEdge[1].nPatches(); ++ii)
        {
            index_t jj = 0;
            for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                if (basisEdge[1].patch(ii).coef(jj, 0) * basisEdge[1].patch(ii).coef(jj, 0) > 1e-25)
                    m_matrix.insert(shift_row + ii, shift_col + j) = basisEdge[1].patch(ii).coef(jj, 0);
            }
        }

        shift_row += basisEdge[0].nPatches();
    }

    // Compute Edge Basis functions
    for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
    {
        const patchSide & bit = m_patches.boundaries()[numBdy];
        index_t side_1 = bit.side().index();
        index_t patch_1 = bit.patch;

        gsApproxC1Edge<d, T> approxC1Edge(m_patches, m_bases, bit, numBdy, m_options);
        std::vector<gsMultiPatch<T>> basisEdge = approxC1Edge.getEdgeBasis();

        index_t begin_col = 0, end_col = 0, shift_col = 0;
        for (index_t np = 0; np < patch_1; ++np)
            shift_col += m_bases[np].size();

//        for (index_t ns = 0; ns < side_1; ++ns)
//            begin_col += m_bases[patch_1].piece(ns).size();
//        for (index_t ns = 0; ns < side_1+1; ++ns)
//            end_col += m_bases[patch_1].piece(ns).size();
        end_col = m_bases[patch_1].piece(0).size();

        for (size_t ii = 0; ii < basisEdge[0].nPatches(); ++ii)
        {
            index_t jj = 0;
            for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                if (basisEdge[0].patch(ii).coef(jj, 0) * basisEdge[0].patch(ii).coef(jj, 0) > 1e-20)
                    m_matrix.insert(shift_row + ii, shift_col + j) = basisEdge[0].patch(ii).coef(jj, 0);
            }
        }
        shift_row += basisEdge[0].nPatches();
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

        gsApproxC1Vertex<d, T> approxC1Vertex(m_patches, m_bases, patchIndex, vertIndex, numVer, m_options);
        std::vector<gsMultiPatch<T>> basisVertex = approxC1Vertex.getVertexBasis();

        for (size_t np = 0; np < patchIndex.size(); ++np)
        {
            index_t patch_1 = patchIndex[np];
            index_t corner = vertIndex[np];

            index_t begin_col = 0, end_col = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
                shift_col += m_bases[np].size();

//            for (index_t ns = 0; ns < corner+4; ++ns)
//                begin_col += m_bases[patch_1].piece(ns).size();
//            for (index_t ns = 0; ns < corner+4+1; ++ns)
//                end_col += m_bases[patch_1].piece(ns).size();
            end_col = m_bases[patch_1].piece(0).size();

            for (size_t ii = 0; ii < basisVertex[np].nPatches(); ++ii)
            {
                index_t jj = 0;
                for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                    if (basisVertex[np].patch(ii).coef(jj, 0) * basisVertex[np].patch(ii).coef(jj, 0) > 1e-20)
                        m_matrix.insert(shift_row + ii, shift_col + j) = basisVertex[np].patch(ii).coef(jj, 0);
                }
            }
        }
        shift_row += basisVertex[0].nPatches(); // + 6

    }
    m_matrix.makeCompressed();
}


} // namespace gismo
