/** @file gsApproxC1Spline.hpp

    @brief Construct the approx c1 spline space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include<gsUnstructuredSplines/gsApproxC1Edge.h>
#include<gsUnstructuredSplines/gsApproxC1Vertex.h>

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

    // For the gluing data space
    m_options.addInt("gluingDataDegree","Polynomial degree of the gluing data space", p-1 );
    m_options.addInt("gluingDataRegularity","Regularity of the gluing data space",  p-2 );
    m_options.addSwitch("info","Print debug information",  false );
    m_options.addSwitch("plot","Print debug information",  false );
}

template<short_t d,class T>
void gsApproxC1Spline<d,T>::init()
{
    p_tilde = m_options.getInt("gluingDataDegree");//math::max(m_optionList.getInt("discreteDegree") - 1, 2);
    r_tilde = m_options.getInt("gluingDataRegularity");//p_tilde - 1;

    m_bases.clear();
    m_bases.reserve(m_patches.nPatches()); // For each Patch
    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        // gsContainerBasisBase:
        // # basisContainer
        // - Interior space: [0] : inner,
        // - Edge spaces:    [1] : west, [2] : east, [3] : south, [4] : north,
        // - Vertex spaces:  [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
        //
        // # helperBasisContainer
        // [0] : basis plus, [1] : basis minus, [2] : basis geo, [3] : basis gluing data
        gsContainerBasis<d,T> containerBasis(9, 4); // for 9 subspaces and 4 helper Basis
        m_bases.push_back(containerBasis);
    }

    // Create interior spline space
    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        gsTensorBSplineBasis<d, T> basis_inner = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(np));

        // Construct special space for r = p - 1:
        // The first and the last knot (not 0,1) are repeated +1, e.g.
        // deg 3, r = 2: |||| || | [...] | || ||||

        index_t m,p;
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
        }
        m_bases[np].setBasis(0, basis_inner); // Inner
    }

    // For loop over Interface to construct the spaces
    for (size_t numInt = 0; numInt < m_patches.interfaces().size(); numInt++)
    {
        const boundaryInterface & item = m_patches.interfaces()[numInt];

        const index_t side_1 = item.first().side().index();
        const index_t side_2 = item.second().side().index();
        const index_t patch_1 = item.first().patch;
        const index_t patch_2 = item.second().patch;

        const index_t dir_1 = side_1 > 2 ? 0 : 1;
        const index_t dir_2 = side_2 > 2 ? 0 : 1;

        gsBSplineBasis<T> basis_1 = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(patch_1).component(dir_1));
        gsBSplineBasis<T> basis_2 = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(patch_2).component(dir_2));

        gsBSplineBasis<T> basis_geo_1 = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(patch_1).component(1-dir_1));
        gsBSplineBasis<T> basis_geo_2 = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(patch_2).component(1-dir_2));

        // [!Plus Minus space]
        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        gsBSplineBasis<T> basis_plus(basis_1);
        gsBSplineBasis<T> basis_minus(basis_1);
        if (m != 1)
            basis_plus.elevateContinuity(1);

        basis_minus.degreeDecrease(1);
        if (m != 1)
            basis_minus.reduceContinuity(1);

        m_bases[patch_1].setHelperBasis(side_1-1, 0, basis_plus);
        m_bases[patch_2].setHelperBasis(side_2-1, 0, basis_plus);

        m_bases[patch_1].setHelperBasis(side_1-1, 1, basis_minus);
        m_bases[patch_2].setHelperBasis(side_2-1, 1, basis_minus);
        // [!Plus Minus space]

        // [!Gluing data space]
        //\tilde{p} = max(p-1,2)
        //\tilde{r} = \tilde{p}-1
        if (basis_1.knots().unique() != basis_2.knots().unique())
            gsInfo << "The patches are not matching!!! \n";
        gsKnotVector<T> kv_gluingData(basis_1.knots().unique(), p_tilde, r_tilde);
        gsBSplineBasis<T> basis_gluingData(kv_gluingData); // S(\tilde{p},\tilde{r},h)

        m_bases[patch_1].setHelperBasis(side_1-1, 3, basis_gluingData);
        m_bases[patch_2].setHelperBasis(side_2-1, 3, basis_gluingData);
        // [!Gluing data space]

        // [!Edge space]
        gsBSplineBasis<T> basis_edge(basis_1);
        basis_edge.setDegreePreservingMultiplicity(basis_plus.degree()+basis_gluingData.degree()-1);

        index_t r_plus, r_minus, r_gD, r_edge, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_gD = p_tilde - basis_gluingData.knots().multiplicityIndex(p_tilde+1);
        r_edge = basis_edge.degree() - basis_edge.knots().multiplicityIndex(basis_edge.degree()+1);

        r = math::min(r_gD, math::min(r_plus, r_minus));
        if (r_edge > r)
            basis_edge.reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.elevateContinuity(r - r_edge);

        gsTensorBSplineBasis<d, T> basis_edge_1(dir_1 == 0 ? basis_edge.knots() : basis_geo_1.knots(),
                                                dir_1 == 0 ? basis_geo_1.knots() : basis_edge.knots());
        gsTensorBSplineBasis<d, T> basis_edge_2(dir_2 == 0 ? basis_edge.knots() : basis_geo_2.knots(),
                                                dir_2 == 0 ? basis_geo_2.knots() : basis_edge.knots());
        // [!Edge space]

        m_bases[patch_1].setHelperBasis(side_1-1, 2, basis_geo_1);
        m_bases[patch_2].setHelperBasis(side_2-1, 2, basis_geo_2);

        m_bases[patch_1].setBasis(side_1, basis_edge_1);
        m_bases[patch_2].setBasis(side_2, basis_edge_2);

    }

    // For loop over the Edge to construct the spaces
    for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
    {
        const patchSide & bit = m_patches.boundaries()[numBdy];

        index_t patch_1 = bit.patch;
        index_t side_1 = bit.side().index();

        index_t dir_1 = m_patches.boundaries()[numBdy].m_index < 3 ? 1 : 0;

        gsBSplineBasis<T> basis_1 = dynamic_cast<gsBSplineBasis<> &>(m_multiBasis.basis(patch_1).component(dir_1));
        gsBSplineBasis<T> basis_geo_1 = dynamic_cast<gsBSplineBasis<> &>(m_multiBasis.basis(patch_1).component(1-dir_1));

        // Assume that plus/minus space is the same as the inner space
        gsKnotVector<T> kv_1 = basis_1.knots();

        gsBSplineBasis<T> patch_basis_1 = dynamic_cast<gsBSplineBasis<> &>(m_patches.patch(patch_1).basis().component(dir_1));
        gsKnotVector<T> kv_patch_1 = patch_basis_1.knots();

        // [!Plus Minus space]
        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        gsBSplineBasis<T> basis_plus(basis_1);
        gsBSplineBasis<T> basis_minus(basis_1);
        if (m != 1)
            basis_plus.elevateContinuity(1);

        basis_minus.degreeDecrease(1);
        if (m != 1)
            basis_minus.reduceContinuity(1);

        m_bases[patch_1].setHelperBasis(side_1-1, 0, basis_plus);
        m_bases[patch_1].setHelperBasis(side_1-1, 1, basis_minus);
        // [!Plus Minus space]

        gsBSplineBasis<T> basis_edge(basis_1);
        basis_edge.setDegreePreservingMultiplicity(basis_plus.degree()-1);

        index_t r_plus, r_minus, r_edge, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_edge = basis_edge.degree() - basis_edge.knots().multiplicityIndex(basis_edge.degree()+1);

        r = math::min(r_plus, r_minus);
        if (r_edge > r)
            basis_edge.reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.elevateContinuity(r - r_edge);

        gsTensorBSplineBasis<d, T> basis_edge_1(dir_1 == 0 ? basis_edge.knots() : basis_geo_1.knots(),
                                                dir_1 == 0 ? basis_geo_1.knots() : basis_edge.knots());

        m_bases[patch_1].setHelperBasis(side_1-1, 2, basis_geo_1);
        m_bases[patch_1].setBasis(side_1, basis_edge_1);
    }

    // For loop over the Vertex to construct the spaces
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

        if (patchIndex.size() == 1) // Boundary vertex
        {
            index_t patch_1 = patchIndex[0];
            index_t vertex_1 = vertIndex[0];

            gsTensorBSplineBasis<d, T> basis_vertex_1 = dynamic_cast<gsTensorBSplineBasis<d, real_t> &>(m_multiBasis.basis(patch_1));

            index_t m, p;
            p = basis_vertex_1.degree(0);
            m = basis_vertex_1.knots(0).multiplicityIndex(p+1);
            if (m == 1) // == basis_vertex_1.degree(1)
                basis_vertex_1.reduceContinuity(1); // In the case for the max. smoothness

            m_bases[patch_1].setBasis(vertex_1+4, basis_vertex_1);
        }
        else if (patchIndex.size() > 1)
        {
            gsMultiPatch<T> temp_mp;
            for (size_t j = 0; j < patchIndex.size(); j++)
                temp_mp.addPatch(m_patches.patch(patchIndex[j]));
            temp_mp.computeTopology();

            if (patchIndex.size() == temp_mp.interfaces().size()) // Internal vertex
            {
                for (size_t j = 0; j < patchIndex.size(); j++)
                {
                    index_t patch_1 = patchIndex[j];
                    index_t vertex_1 = vertIndex[j];

                    gsTensorBSplineBasis<d, T> basis_vertex_1 = dynamic_cast<gsTensorBSplineBasis<d, real_t> &>(m_multiBasis.basis(
                            patch_1));

                    basis_vertex_1.degreeElevate(p_tilde-1,0); // Keep smoothness
                    basis_vertex_1.degreeElevate(p_tilde-1,1);

                    // todo: fix for general regularity
                    index_t r, p;
                    p = basis_vertex_1.degree(0);
                    r = p - basis_vertex_1.knots(0).multiplicityIndex(p+1);

                    //if (m_multiBasis.basis(patch_1).degree(0) - r == 1)
                    //    basis_vertex_1.reduceContinuity(r-1); // In the case for the max. smoothness
                    //else if (r > 2)
                    {
                        if (r != 1)
                            basis_vertex_1.reduceContinuity(1); // bcs of minus space
                        if (r_tilde < r-1)
                            basis_vertex_1.reduceContinuity(r-r_tilde-1);
                    }

                    m_bases[patch_1].setBasis(vertex_1+4, basis_vertex_1);
                }
            }
            else if (patchIndex.size() > temp_mp.interfaces().size())// Interface-Boundary vertex
            {
                for (size_t j = 0; j < patchIndex.size(); j++)
                {
                    index_t patch_1 = patchIndex[j];
                    index_t vertex_1 = vertIndex[j];

                    gsTensorBSplineBasis<d, T> basis_vertex_1 = dynamic_cast<gsTensorBSplineBasis<d, real_t> &>(m_multiBasis.basis(patch_1));

                    basis_vertex_1.degreeElevate(p_tilde-1,0); // Keep smoothness
                    basis_vertex_1.degreeElevate(p_tilde-1,1);

                    // todo: fix for general regularity
                    index_t r, p;
                    p = basis_vertex_1.degree(0);
                    r = p - basis_vertex_1.knots(0).multiplicityIndex(p+1);

                    //if (m_multiBasis.basis(patch_1).degree(0) - r == 1)
                    //    basis_vertex_1.reduceContinuity(r-1); // In the case for the max. smoothness
                    //else if (r > 2)
                    {
                        if (r != 1)
                            basis_vertex_1.reduceContinuity(1); // bcs of minus space
                        if (r_tilde < r-1)
                            basis_vertex_1.reduceContinuity(r-r_tilde-1);
                    }

                    m_bases[patch_1].setBasis(vertex_1+4, basis_vertex_1);

                }
            }
        }
    }


    index_t row_dofs = 0;
    // Inner basis
    for (size_t np = 0; np < m_patches.nPatches(); np++)
    {
        index_t dim_u = m_bases[np].getBasis(0).component(0).size();
        index_t dim_v = m_bases[np].getBasis(0).component(1).size();
        row_dofs += (dim_u - 4) * (dim_v - 4);
    }

    // Interfaces
    for (size_t numInt = 0; numInt < m_patches.interfaces().size(); numInt++)
    {
        const boundaryInterface &item = m_patches.interfaces()[numInt];

        const index_t side_1 = item.first().side().index();
        const index_t patch_1 = item.first().patch;

        index_t numDofs = math::max(m_bases[patch_1].getHelperBasis(side_1-1, 0).size() + m_bases[patch_1].getHelperBasis(side_1-1, 1).size() - 10, 0);
        row_dofs += numDofs; // The same as side_2
    }

    // Boundary Edges
    for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
    {
        const patchSide &bit = m_patches.boundaries()[numBdy];

        index_t patch_1 = bit.patch;
        index_t side_1 = bit.side().index();

        index_t numDofs = math::max(m_bases[patch_1].getHelperBasis(side_1-1, 0).size() + m_bases[patch_1].getHelperBasis(side_1-1, 1).size() - 10, 0);
        row_dofs += numDofs;
    }

    // Vertices
    for (size_t numVer = 0; numVer < m_patches.vertices().size(); numVer++)
        row_dofs += 6;



    m_matrix.clear();
    index_t dim_col = 0;
    for (size_t i = 0; i < m_bases.size(); i++)
    {
        dim_col += m_bases[i].size();
    }

    m_matrix.resize(row_dofs, dim_col);
    const index_t nz = (m_multiBasis.basis(0).degree(0)*2 + 1)*row_dofs; // TODO
    m_matrix.reserve(nz);
}   


template<short_t d,class T>
void gsApproxC1Spline<d,T>::compute()
{
    // Compute Inner Basis functions
    index_t shift_row = 0, shift_col = 0;
    for(size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        index_t dim_u = m_bases[np].getBasis(0).component(0).size();
        index_t dim_v = m_bases[np].getBasis(0).component(1).size();

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

        for (index_t ns = 0; ns < side_1; ++ns)
            begin_col += m_bases[patch_1].getBasis(ns).size();
        for (index_t ns = 0; ns < side_1+1; ++ns)
            end_col += m_bases[patch_1].getBasis(ns).size();

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

        for (index_t ns = 0; ns < side_2; ++ns)
            begin_col += m_bases[patch_2].getBasis(ns).size();
        for (index_t ns = 0; ns < side_2+1; ++ns)
            end_col += m_bases[patch_2].getBasis(ns).size();

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

        for (index_t ns = 0; ns < side_1; ++ns)
            begin_col += m_bases[patch_1].getBasis(ns).size();
        for (index_t ns = 0; ns < side_1+1; ++ns)
            end_col += m_bases[patch_1].getBasis(ns).size();

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

            for (index_t ns = 0; ns < corner+4; ++ns)
                begin_col += m_bases[patch_1].getBasis(ns).size();
            for (index_t ns = 0; ns < corner+4+1; ++ns)
                end_col += m_bases[patch_1].getBasis(ns).size();

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
