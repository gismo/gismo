/** @file gsApproxC1Vertex.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines2/gsApproxC1Vertex.h>

namespace gismo
{

template<short_t d,class T>
void gsApproxC1Vertex<d, T>::reparametrizeVertexPatches()
{
    for(size_t i = 0; i < m_auxPatches.size(); i++)
    {
        checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

        index_t vertex_1 = m_vertexIndices[i];
        if(m_auxPatches[i].getOrient() == 0) // not switched
            switch (vertex_1) // == vertex
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
        else if (m_auxPatches[i].getOrient() == 1) // switched
            switch (vertex_1) // == vertex
            {
                case 1:
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                    break;
                case 4:
                    m_auxPatches[i].rotateParamAntiClockTwice();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                    break;
                case 3:
                    m_auxPatches[i].rotateParamAntiClock();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                    break;
                case 2:
                    m_auxPatches[i].rotateParamClock();
                    //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                    break;
            }
    }
}

template<short_t d,class T>
real_t gsApproxC1Vertex<d, T>::computeSigma(const std::vector<size_t> &vertexIndices)
{
    real_t sigma = 0;

    real_t p = 0;
    real_t h_geo = 1;
    for(size_t i = 0; i < m_auxPatches.size(); i++)
    {
        const gsTensorBSplineBasis<2, real_t> bsp_temp = dynamic_cast<const gsTensorBSplineBasis<d, T>&>(m_auxPatches[0].getBasisRotated().piece(0));

        real_t p_temp = math::max(bsp_temp.degree(0), bsp_temp.degree(1));

        p = (p < p_temp ? p_temp : p);

        for(index_t j = 0; j < m_auxPatches[i].getPatchRotated().parDim(); j++)
        {
            real_t h_geo_temp = bsp_temp.knot(j,p + 1);
            h_geo = (h_geo > h_geo_temp ? h_geo_temp : h_geo);
        }
    }

    gsMatrix<> zero;
    zero.setZero(2,1);
    for (size_t i = 0; i < m_auxPatches.size(); i++)
        sigma += gsMatrix<> (m_auxPatches[i].getPatchRotated().deriv(zero)).lpNorm<Eigen::Infinity>();
    sigma *= h_geo/(m_auxPatches.size()*p);

    return (1.0 / sigma);
}

template<short_t d,class T>
void gsApproxC1Vertex<d, T>::checkOrientation(size_t i) {
    if (m_auxPatches[i].getPatchRotated().orientation() == -1)
    {
        m_auxPatches[i].swapAxis();
        //gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";
    }
}

template<short_t d,class T>
void gsApproxC1Vertex<d, T>::computeKernel()
{
    // FOR NEUMANN: SOME TESTS:
    std::vector<patchSide> containingSides;
    patchCorner pC(m_patchesAroundVertex[0], m_vertexIndices[0]);
    pC.getContainingSides(d, containingSides);
    // TEST END

    gsMultiPatch<T> mp_vertex;
    for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        mp_vertex.addPatch(m_auxPatches[i].getPatchRotated());

    mp_vertex.computeTopology();

    index_t dim_mat = 0;
    std::vector<index_t> dim_u, dim_v, side, patchID;
    std::vector<index_t> dim_u_iFace, patchID_iFace;
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
            if (!m_optionList.getSwitch("second"))
                dim_mat += basisVertexResult[np].basis(0).component(0).size();

            matrix_det.col(0) = m_auxPatches[np].getPatchRotated().jacobian(points).col(0); // u
        }
        if (mp_vertex.isBoundary(np,1)) // v
        {
            side.push_back(1);
            patchID.push_back(np);
            dim_u.push_back(basisVertexResult[np].basis(0).component(0).size());
            dim_v.push_back(basisVertexResult[np].basis(0).component(1).size());
            dim_mat += basisVertexResult[np].basis(0).component(1).size();
            if (!m_optionList.getSwitch("second"))
                dim_mat += basisVertexResult[np].basis(0).component(1).size();

            matrix_det.col(1) = m_auxPatches[np].getPatchRotated().jacobian(points).col(1); // u
        }
        if(mp_vertex.isInterface(patchSide(np,1)) && mp_vertex.isInterface(patchSide(np,3)))
        {
            if (!m_optionList.getSwitch("second"))
                dim_mat += 4;
            else
                dim_mat += 1;
            patchID_iFace.push_back(np);
            dim_u_iFace.push_back(basisVertexResult[np].basis(0).component(0).size());
        }

    }
    if (patchID.size() != 2)
        gsInfo << "Something went wrong \n";

    index_t dofsCorner = 1;
    if (m_optionList.getSwitch("second"))
        dofsCorner = 3;  // No Neumann

    if (matrix_det.determinant()*matrix_det.determinant() > 1e-15) // There is (numerically) a kink
    {
        if (!m_optionList.getSwitch("second"))
            dofsCorner = 0;  // With Neumann
        else
            dofsCorner = 1;
    }

    if (m_optionList.getSwitch("info"))
        gsInfo << "Det: " << matrix_det.determinant() << "\n";

    gsMatrix<T> coefs_corner(dim_mat, 6);
    coefs_corner.setZero();

    index_t shift_row = 0;
    for (size_t bdy_index = 0; bdy_index < patchID.size(); ++bdy_index)
    {
        if (side[bdy_index] < 3) // v
        {
            index_t shift_row_neumann = dim_v[bdy_index];
            for (index_t i = 0; i < dim_v[bdy_index]; ++i)
            {
                for (index_t j = 0; j < 6; ++j)
                {
                    T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i*dim_u[bdy_index], 0); // v = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_corner(shift_row+i, j) = coef_temp;
                    if (!m_optionList.getSwitch("second"))
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i*dim_u[bdy_index] +1, 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row_neumann + shift_row+i, j) = coef_temp;
                    }
                }
            }
            shift_row += dim_v[bdy_index];
            if (!m_optionList.getSwitch("second"))
                shift_row += dim_v[bdy_index];
        }
        else // u
        {
            index_t shift_row_neumann = dim_u[bdy_index];
            for (index_t i = 0; i < dim_u[bdy_index]; ++i)
            {
                for (index_t j = 0; j < 6; ++j)
                {
                    T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i, 0); // v = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_corner(shift_row+i, j) = coef_temp;
                    if (!m_optionList.getSwitch("second"))
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i+dim_u[bdy_index], 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row_neumann + shift_row+i, j) = coef_temp;
                    }
                }
            }
            shift_row += dim_u[bdy_index];
            if (!m_optionList.getSwitch("second"))
                shift_row += dim_u[bdy_index];
        }
    }
    for (size_t iFace_index = 0; iFace_index < patchID_iFace.size(); ++iFace_index)
    {
        if (!m_optionList.getSwitch("second")) {
            for (index_t i = 0; i < 2; ++i) // Only the first two
            {
                for (index_t j = 0; j < 6; ++j) {
                    T coef_temp = basisVertexResult[patchID[iFace_index]].patch(j).coef(i, 0); // v = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_corner(shift_row + i, j) = coef_temp;
                    coef_temp = basisVertexResult[patchID[iFace_index]].patch(j).coef(i + dim_u[iFace_index],
                                                                                      0); // v = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_corner(shift_row + 2 + i, j) = coef_temp; //  +2 bcs of the previous adding
                }
            }
            shift_row += 4;
        }
        else
        {
            for (index_t j = 0; j < 6; ++j)
            {
                T coef_temp = basisVertexResult[patchID[iFace_index]].patch(j).coef(0, 0); // v = 0
                if (coef_temp * coef_temp > 1e-25)
                    coefs_corner(shift_row, j) = coef_temp;
            }
            shift_row += 1;
        }
    }

    gsMatrix<T> kernel;
    if (dofsCorner > 0)
    {
        real_t threshold = 1e-10;
        Eigen::FullPivLU<gsMatrix<>> KernelCorner(coefs_corner);
        KernelCorner.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelCorner.dimensionOfKernel() < dofsCorner) {
            threshold += 1e-8;
            KernelCorner.setThreshold(threshold);
        }
        if (m_optionList.getSwitch("info"))
            gsInfo << "Dimension of Kernel: " << KernelCorner.dimensionOfKernel() << " With " << threshold << "\n";

        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);

        kernel = KernelCorner.kernel();

        size_t count = 0;
        while (kernel.cols() < 6) {
            kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
            kernel.col(kernel.cols() - 1) = vertBas.col(count);

            Eigen::FullPivLU<gsMatrix<>> ker_temp(kernel);
            ker_temp.setThreshold(1e-6);
            if (ker_temp.dimensionOfKernel() != 0) {
                kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
            }
            count++;
        }
    }
    else
        kernel.setIdentity(6, 6);

    if (m_optionList.getSwitch("info"))
        gsInfo << "NumDofs: " << dofsCorner << " with Kernel: \n" << kernel << "\n";

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

} // namespace gismo