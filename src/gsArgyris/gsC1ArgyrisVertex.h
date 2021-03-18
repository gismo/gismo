/** @file gsC1ArgyrisVertex.h

    @brief Creates the C1 Argyris Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsArgyris/gsC1ArgyrisAuxiliaryPatch.h>
#include <gsArgyris/gsApproxArgyrisVertexBasis.h>



namespace gismo {
template<short_t d, class T>
class gsC1ArgyrisVertex
{

private:
    typedef gsC1ArgyrisBasis<d, T> Basis;
    typedef typename std::vector<Basis> ArgyrisBasisContainer;
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

    /// Shared pointer for gsC1Argyris
    typedef memory::shared_ptr<gsC1ArgyrisVertex> Ptr;

    /// Unique pointer for gsC1Argyris
    typedef memory::unique_ptr<gsC1ArgyrisVertex> uPtr;


public:
    /// Empty constructor
    ~gsC1ArgyrisVertex() { }


    gsC1ArgyrisVertex(gsMultiPatch<T> const & mp,
                ArgyrisBasisContainer & bases,
                const std::vector<size_t> & patchesAroundVertex,
                const std::vector<size_t> & vertexIndices,
                const index_t & numVer,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_patchesAroundVertex(patchesAroundVertex),
                m_vertexIndices(vertexIndices), m_optionList(optionList)
    {
        m_auxPatches.clear();
        basisVertexResult.clear();

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            index_t vertex_1 = m_vertexIndices[i];
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsC1ArgyrisAuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        reparametrizeVertexPatches();

        if (m_optionList.getSwitch("twoPatch") && m_patchesAroundVertex.size() == 1)
        {
            gsMultiPatch<> g1Basis;
            gsTensorBSplineBasis<d, T> basis_edge = m_auxPatches[0].getArygrisBasisRotated().getVertexBasis(
                    m_vertexIndices[0]); // 0 -> u, 1 -> v

            index_t dim_u = basis_edge.component(0).size();
            index_t dim_v = basis_edge.component(1).size();
            for (index_t j = 0; j < 2; j++) // v
            {
                for (index_t i = 0; i < 2; i++) // u
                {
                    gsMatrix<> coefs;
                    coefs.setZero(dim_u * dim_v, 1);

                    coefs(j * dim_u + i, 0) = 1;

                    g1Basis.addPatch(basis_edge.makeGeometry(coefs));
                }
            }
            m_auxPatches[0].parametrizeBasisBack(g1Basis);

            basisVertexResult.push_back(g1Basis);
        }
        else if (!m_optionList.getSwitch("twoPatch"))
        {
            // Compute Sigma
            real_t sigma = computeSigma(m_vertexIndices);

            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
            {
                ArgyrisAuxPatchContainer auxPatchSingle;
                auxPatchSingle.push_back(m_auxPatches[i]);

                std::vector<index_t> sideContainer;
                switch (auxPatchSingle[0].side())
                {
                    case 1:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 2:
                        sideContainer.push_back(2); // u
                        sideContainer.push_back(3); // v
                        break;
                    case 3:
                        sideContainer.push_back(1); // u
                        sideContainer.push_back(4); // v
                        break;
                    case 4:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(2); // v
                        break;
                    default:
                        gsInfo << "Something went wrong\n";
                        break;
                }

                // Compute Gluing data
                gsApproxGluingData<d, T> approxGluingData(auxPatchSingle, m_optionList, sideContainer);
                // Create Basis functions
                gsApproxArgyrisVertexBasis<d, T> approxArgyrisVertexBasis(auxPatchSingle, approxGluingData,
                                                                          m_vertexIndices[i], sideContainer, sigma, m_optionList);

                gsMultiPatch<> result_1;
                approxArgyrisVertexBasis.setBasisVertex(result_1);

                // Store temporary
                basisVertexResult.push_back(result_1);
            }

            if (m_auxPatches[0].getArygrisBasisRotated().getKindOfVertex(m_vertexIndices[0]) != 0) // No internal vertex
            {
                computeKernel();

                gsMatrix<> points;
                points.setZero(2,2);
                points(0,1) = 1.0;


                if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("twoPatch"))
                    for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                        for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                            gsInfo << i << " : " << basisVertexResult[np].patch(i).deriv(points.col(0)) << "\n\n";

                for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                    m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack

                gsInfo << "\n";

                if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("twoPatch"))
                    for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                        for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                            gsInfo << i << " : " << basisVertexResult[np].patch(i).deriv(points.col(1-np)) << "\n\n";

            }
            else // Internal vertex
            {
                for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                    m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
            }


        }

        /*
        gsMatrix<> points;
        points.setZero(2,2);
        points(0,1) = 1.0;


        if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("twoPatch"))
            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                    gsInfo << basisVertexResult[np].patch(i).deriv(points) << "\n\n";
        */

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
            {
                for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                {
                    fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                    gsField<> temp_field(m_mp.patch(m_patchesAroundVertex[np]), basisVertexResult[np].patch(i));
                    gsWriteParaview(temp_field, fileName, 5000);
                    collection.addTimestep(fileName, i, "0.vts");

                }
            }
            collection.save();
        }


    }

    void saveBasisVertex(gsSparseMatrix<T> & system)
    {
        for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
        {
            index_t patch_1 = m_patchesAroundVertex[np];
            index_t corner = m_vertexIndices[np];

            index_t shift_row = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
            {
                shift_row += m_bases[np].size_rows();
                shift_col += m_bases[np].size_cols();
            }

            index_t ii = 0;
            for (index_t i = m_bases[patch_1].rowBegin(corner+4); i < m_bases[patch_1].rowEnd(corner+4); ++i, ++ii)
            {
                index_t jj = 0;
                for (index_t j = m_bases[patch_1].colBegin(corner+4); j < m_bases[patch_1].colEnd(corner+4); ++j, ++jj)
                    if (basisVertexResult[np].patch(ii).coef(jj,0)*basisVertexResult[np].patch(ii).coef(jj,0)>1e-25)
                        system.insert(shift_row+i,shift_col+j) = basisVertexResult[np].patch(ii).coef(jj,0);
            }
        }
    }

    void reparametrizeVertexPatches()
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            switch (m_auxPatches[i].side()) // == vertex
            {
                case 1:
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                    break;
                case 4:
                    m_auxPatches[i].rotateParamAntiClockTwice();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                    break;
                case 2:
                    m_auxPatches[i].rotateParamAntiClock();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                    break;
                case 3:
                    m_auxPatches[i].rotateParamClock();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                    break;
            }
        }
    }

    void checkOrientation(size_t i)
    {
        if (m_auxPatches[i].getPatch().orientation() == -1)
        {
            m_auxPatches[i].swapAxis();
            //gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";

            // Swap vertices index after swapping axis
            if(m_auxPatches[i].side() == 2)
                m_auxPatches[i].setSide(3);
            else if(m_auxPatches[i].side() == 3)
                m_auxPatches[i].setSide(2);
        }
    }

    real_t computeSigma(const std::vector<size_t> & vertexIndices)
    {
        real_t sigma = 0;

        real_t p = 0;
        real_t h_geo = 0;
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            gsTensorBSplineBasis<2, real_t> bsp_temp = m_auxPatches[0].getArygrisBasisRotated().getVertexBasis(vertexIndices[i]);

            real_t p_temp = math::max(bsp_temp.degree(0), bsp_temp.degree(1));

            p = (p < p_temp ? p_temp : p);

            for(index_t j = 0; j < m_auxPatches[i].getPatch().parDim(); j++)
            {
                real_t h_geo_temp = bsp_temp.knot(j,p + 2);
                h_geo = (h_geo < h_geo_temp ? h_geo_temp : h_geo);
            }
        }

        gsMatrix<> zero;
        zero.setZero(2,1);
        for (size_t i = 0; i < m_auxPatches.size(); i++)
            sigma += gsMatrix<> (m_auxPatches[i].getPatch().deriv(zero)).lpNorm<Eigen::Infinity>();
        sigma *= h_geo/(m_auxPatches.size()*p);

        return (1 / sigma);
    }

    void computeKernel()
    {

        // TODO Boundary vertex with valence > 2
        gsMultiPatch<T> mp_vertex;
        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            mp_vertex.addPatch(m_auxPatches[i].getPatch());
        }
        mp_vertex.computeTopology();

        index_t dim_mat = 0;
        std::vector<index_t> dim_u, dim_v, side, patchID;
        gsMatrix<> matrix_det(2,2), points(2,1);
        points.setZero();
        for(size_t np = 0; np < mp_vertex.nPatches(); np++)
        {
            if (mp_vertex.isBoundary(np,3)) // u
            {
                side.push_back(3);
                patchID.push_back(np);
                dim_u.push_back(basisVertexResult[np].basis(0).component(0).size());
                dim_v.push_back(basisVertexResult[np].basis(0).component(1).size());
                dim_mat += basisVertexResult[np].basis(0).component(0).size();

                matrix_det.col(0) = m_auxPatches[np].getPatch().jacobian(points).col(0); // u
            }
            if (mp_vertex.isBoundary(np,1)) // v
            {
                side.push_back(1);
                patchID.push_back(np);
                dim_u.push_back(basisVertexResult[np].basis(0).component(0).size());
                dim_v.push_back(basisVertexResult[np].basis(0).component(1).size());
                dim_mat += basisVertexResult[np].basis(0).component(1).size();

                matrix_det.col(1) = m_auxPatches[np].getPatch().jacobian(points).col(1); // u
            }
        }
        if (patchID.size() != 2)
            gsInfo << "Something went wrong \n";

        index_t dofsCorner = 3;
        if (matrix_det.determinant()*matrix_det.determinant() > 0) // There is a kink
            dofsCorner = 1;

        for(size_t np = 0; np < mp_vertex.nPatches(); np++)
            m_bases[m_patchesAroundVertex[np]].setNumDofsVertex(dofsCorner, m_vertexIndices[np]);

        gsInfo << "Det: " << matrix_det.determinant() << "\n";

        gsMatrix<T> coefs_corner(dim_mat, 6);
        coefs_corner.setZero();

        index_t shift_row = 0;
        for (size_t bdy_index = 0; bdy_index < patchID.size(); ++bdy_index)
        {
            if (side[bdy_index] < 3) // v
            {
                for (index_t i = 0; i < dim_v[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i*dim_u[bdy_index], 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row+i, j) = coef_temp;
                    }
                }
                shift_row += dim_v[bdy_index];
            }
            else // u
            {
                for (index_t i = 0; i < dim_u[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i, 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row+i, j) = coef_temp;
                    }
                }
                shift_row += dim_u[bdy_index];
            }
        }

        real_t threshold = 1e-8;
        Eigen::FullPivLU<gsMatrix<>> KernelCorner(coefs_corner);
        KernelCorner.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelCorner.dimensionOfKernel() < dofsCorner)
        {
            threshold += 1e-8;
            KernelCorner.setThreshold(threshold);
        }
        gsInfo << "Dimension of Kernel: " << KernelCorner.dimensionOfKernel() << " With " << threshold << "\n";

        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);

        gsMatrix<T> kernel = KernelCorner.kernel();

        size_t count = 0;
        while (kernel.cols() < 6) {
            kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
            kernel.col(kernel.cols() - 1) = vertBas.col(count);

            Eigen::FullPivLU<gsMatrix<>> ker_temp(kernel);
            ker_temp.setThreshold(1e-5);
            if (ker_temp.dimensionOfKernel() != 0) {
                kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
            }
            count++;
        }
        gsInfo << "Kernel: " << kernel << "\n";

        for(size_t np = 0; np < m_patchesAroundVertex.size(); np++)
        {
            gsMultiPatch<> temp_result_0 = basisVertexResult[np];

            for (size_t j = 0; j < 6; ++j)
            {
                index_t dim_uv = temp_result_0.basis(j).size();
                gsMatrix<> coef_bf;
                coef_bf.setZero(dim_uv, 1);
                for (index_t i = 0; i < 6; ++i)
                    if (kernel(i, j) * kernel(i, j) > 1e-25)
                        coef_bf += temp_result_0.patch(i).coefs() * kernel(i, j);

                basisVertexResult[np].patch(j).setCoefs(coef_bf);
            }
        }

    }


protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    ArgyrisBasisContainer & m_bases;

    const std::vector<size_t> & m_patchesAroundVertex;
    const std::vector<size_t> & m_vertexIndices;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    ArgyrisAuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

}; // gsC1ArgyrisVertex

} // namespace gismo
