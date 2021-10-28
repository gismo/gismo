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
        /*
            to do: general
                m_options.addInt("gluingDataDegree","Polynomial degree of the gluing data space", p-1 );
                m_options.addInt("gluingDataRegularity","Regularity of the gluing data space",  p-2 );
                m_options.addSwitch("info","Print debug information",  false );
                m_options.addSwitch("plot","Print debug information",  false );
        */
    }

    template<short_t d,class T>
    void gsC1SurfSpline<d,T>::init()
    {
        m_bases.clear();
        m_bases.reserve(m_patches.nPatches()); // For each Patch
        for (size_t np = 0; np < m_patches.nPatches(); np++)
        {
            // gsContainerBasisBase
            gsContainerBasis<d,T> containerBasis(1, 0); // for 9 subspaces and 4 helper Basis
            m_bases.push_back(containerBasis);
        }


        // Create interior spline space
        for (size_t np = 0; np < m_patches.nPatches(); np++)
        {
            gsTensorBSplineBasis<d, T> basis_patch = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_multiBasis.basis(np));
            gsInfo << "Basis Patch: " << basis_patch.component(0).knots() << "\n";
            m_bases[np].setBasis(0, basis_patch); // Inner
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
            index_t dir = m_patches.interfaces()[numInt].first().m_index < 3 ? 1 : 0;
            gsBSplineBasis<T> basis_edge = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(m_patches.interfaces()[numInt].first().patch).component(dir)); // If the interface matches!!!
            index_t m_p = basis_edge.maxDegree();
            index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
            index_t m_n = basis_edge.numElements();

            index_t numDofs = 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p - 9;
            row_dofs += numDofs;
        }

        // Boundary Edges
        for (size_t numBdy = 0; numBdy < m_patches.boundaries().size(); numBdy++)
        {
            const patchSide &bit = m_patches.boundaries()[numBdy];

            index_t patch_1 = bit.patch;
            index_t side_1 = bit.side().index();

            index_t dir = side_1 < 3 ? 1 : 0;
            gsBSplineBasis<T> basis_edge = dynamic_cast<gsBSplineBasis<T> &>(m_multiBasis.basis(patch_1).component(dir)); // If the interface matches!!!
            index_t m_p = basis_edge.maxDegree();
            index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
            index_t m_n = basis_edge.numElements();

            index_t numDofs = 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p - 9;
            row_dofs += numDofs;
        }

        // Vertices
        for (size_t numVer = 0; numVer < m_patches.vertices().size(); numVer++)
        {
            row_dofs += 6;
        }


        m_matrix.clear();
        index_t dim_col = 0;
        for (size_t i = 0; i < m_bases.size(); i++)
        {
            dim_col += m_bases[i].size();
        }

        m_matrix.resize(row_dofs, dim_col);
        const index_t nz = 7*row_dofs; // TODO
        m_matrix.reserve(nz);

        gsInfo << m_multiBasis << "\n";
        for (size_t np = 0; np < m_patches.nPatches(); np++)
            gsInfo << "Basis " << np << " Size: " << m_multiBasis.basis(np).size() << "\n";
        gsInfo << "Mat dim: (" << row_dofs  << " : " << dim_col << ")\n";
    }


    template<short_t d,class T>
    void gsC1SurfSpline<d,T>::compute()
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
            gsC1SurfEdge<d,T> smoothC1Edge(m_patches,m_patches.interfaces()[numInt]);
            smoothC1Edge.computeG1InterfaceBasis();
            std::vector<gsMultiPatch<T>> basisEdge = smoothC1Edge.getBasis();

            gsWriteParaview(basisEdge[0].patch(0),"test_1",2000);

            const boundaryInterface & item = m_patches.interfaces()[numInt];
            index_t side_1 = item.first().side().index();
            index_t side_2 = item.second().side().index();

            index_t patch_1 = item.first().patch;
            index_t patch_2 = item.second().patch;

            index_t begin_col = 0, end_col = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
                shift_col += m_bases[np].size();

            end_col += m_bases[patch_1].getBasis(0).size();
            for (index_t ii = 0; ii < basisEdge[0].nPatches(); ++ii)
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

            end_col += m_bases[patch_2].getBasis(0).size();

            for (index_t ii = 0; ii < basisEdge[1].nPatches(); ++ii)
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

            index_t patch_1 = bit.patch;
            index_t side_1 = bit.side().index();

            gsC1SurfEdge<d,T> smoothC1Edge(m_patches,bit);
            smoothC1Edge.computeG1BoundaryBasis(side_1);
            std::vector<gsMultiPatch<T>> basisEdge = smoothC1Edge.getBasis();

            index_t begin_col = 0, end_col = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
                shift_col += m_bases[np].size();

            end_col += m_bases[patch_1].getBasis(0).size();
            for (index_t ii = 0; ii < basisEdge[0].nPatches(); ++ii)
            {
                index_t jj = 0;
                for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                    if (basisEdge[0].patch(ii).coef(jj, 0) * basisEdge[0].patch(ii).coef(jj, 0) > 1e-25)
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

            gsC1SurfVertex<d,T> smoothC1Vertex(m_patches, patchIndex, vertIndex);
            smoothC1Vertex.computeG1InternalVertexBasis();
            std::vector<gsMultiPatch<T>> basisVertex = smoothC1Vertex.getBasis();

            for(size_t pInd=0; pInd < patchIndex.size(); pInd++)
            {
                index_t begin_col = 0, end_col = 0, shift_col = 0;
                for (index_t np = 0; np < patchIndex[pInd]; ++np)
                    shift_col += m_bases[np].size();

                end_col += m_bases[patchIndex[pInd]].getBasis(0).size();
                for (index_t ii = 0; ii < basisVertex[pInd].nPatches(); ++ii)
                {
                    index_t jj = 0;
                    for (index_t j = begin_col; j < end_col; ++j, ++jj) {
                        if (basisVertex[pInd].patch(ii).coef(jj, 0) * basisVertex[pInd].patch(ii).coef(jj, 0) > 1e-25)
                            m_matrix.insert(shift_row + ii, shift_col + j) = basisVertex[pInd].patch(ii).coef(jj, 0);
                    }
                }
            }
            shift_row += 6;
        }
        m_matrix.makeCompressed();

        gsInfo << m_matrix << "\n";
    }


} // namespace gismo
