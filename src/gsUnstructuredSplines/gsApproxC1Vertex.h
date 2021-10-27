/** @file gsApproxC1Vertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsContainerBasis.h>
#include <gsUnstructuredSplines/gsPatchReparameterized.h>
#include <gsUnstructuredSplines/gsApproxC1VertexBasisProjection.h>



namespace gismo {
template<short_t d, class T>
class gsApproxC1Vertex
{

private:
    typedef gsContainerBasis<d, T> Basis;
    typedef typename std::vector<Basis> C1BasisContainer;
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Vertex
    typedef memory::shared_ptr<gsApproxC1Vertex> Ptr;

    /// Unique pointer for gsApproxC1Vertex
    typedef memory::unique_ptr<gsApproxC1Vertex> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Vertex() { }


    gsApproxC1Vertex(gsMultiPatch<T> & mp,
                C1BasisContainer & bases,
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

            m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        reparametrizeVertexPatches();

        // Compute Sigma
        real_t sigma = computeSigma(m_vertexIndices);

        std::vector<gsApproxGluingData<d, T>> gD; // delete later

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            C1AuxPatchContainer auxPatchSingle;
            auxPatchSingle.push_back(m_auxPatches[i]);

            std::vector<index_t> sideContainer;
            if (auxPatchSingle[0].getOrient() == 0) // not rotated
            {
                switch (auxPatchSingle[0].side()) // corner
                {
                    case 1:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 2:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(2); // v
                        break;
                    case 3:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 4:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(2); // v
                        break;
                    default:
                        gsInfo << "Something went wrong\n";
                        break;
                }
            }
            else if (auxPatchSingle[0].getOrient() == 1) // rotated
            {
                switch (auxPatchSingle[0].side()) // corner
                {
                    case 1:
                        sideContainer.push_back(1); // u
                        sideContainer.push_back(3); // v
                        break;
                    case 2:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(2); // v
                        break;
                    case 3:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 4:
                        sideContainer.push_back(2); // u
                        sideContainer.push_back(4); // v
                        break;
                    default:
                        gsInfo << "Something went wrong\n";
                        break;
                }
            }

            //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with index: " << m_vertexIndices[i] << "\n";
            std::vector<bool> isInterface;
            isInterface.resize(2);
            isInterface[0] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i],sideContainer[0]));
            isInterface[1] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i],sideContainer[1]));

            // Compute Gluing data
            gsApproxGluingData<d, T> approxGluingData(auxPatchSingle, m_optionList, sideContainer, isInterface);

            // Create Basis functions
            gsMultiPatch<> result_1;

            gsApproxC1VertexBasisProjection<d, T> approxC1VertexBasis(auxPatchSingle, approxGluingData,
                                                                                m_vertexIndices[i], sideContainer, isInterface, sigma, m_optionList);
            approxC1VertexBasis.setBasisVertex(result_1);


            // Store temporary
            basisVertexResult.push_back(result_1);

            gD.push_back(approxGluingData); // delete later
        }

        gsMultiPatch<T> temp_mp;
        for (size_t j = 0; j < m_patchesAroundVertex.size(); j++)
            temp_mp.addPatch(m_mp.patch(m_patchesAroundVertex[j]));
        temp_mp.computeTopology();

        if (m_patchesAroundVertex.size() != temp_mp.interfaces().size()) // No internal vertex
        {
            computeKernel();

            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
        }
        else // Internal vertex
        {
            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
        }

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
            {
                if (basisVertexResult.size() != 0)
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

    gsApproxC1Vertex(gsMultiPatch<T> const & mp,
                      C1BasisContainer & bases,
                      const std::vector<size_t> & patchesAroundVertex,
                      const std::vector<size_t> & vertexIndices,
                      const index_t & numVer,
                      std::vector<std::vector<gsMultiPatch<T>>> & vertex_bf,
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

            m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        // 2 == u,v
        std::vector<std::vector<gsMultiPatch<T>>> vertexAuxiliary_bf(patchesAroundVertex.size(), std::vector<gsMultiPatch<T>>(2));
        reparametrizeVertexBasis(vertex_bf, vertexAuxiliary_bf);

        gsMultiPatch<T> mp_vertex;
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
            mp_vertex.addPatch(m_mp.patch(patchesAroundVertex[i]));
        mp_vertex.computeTopology();

        gsMatrix<T> coeffs_mat;
        coeffs_mat.setZero(mp_vertex.nPatches()*4, mp_vertex.nInterfaces()*5 ); // *8, + mp_vertex.nPatches()*4

        gsMatrix<> points(2, 1);
        points.setZero();

        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            for(index_t dir = 0; dir < 2; dir++)
            {
                index_t i_shift = 0;
                for(size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++)
                {
                    if (mp_vertex.interfaces()[numInt].first().patch == (index_t) i)
                    {
                        if (mp_vertex.interfaces()[numInt].first().side().index() > 2 && dir == 0)
                            i_shift = numInt*5;
                        else if (mp_vertex.interfaces()[numInt].first().side().index() < 3 && dir == 1)
                            i_shift = numInt*5;
                    }
                    if (mp_vertex.interfaces()[numInt].second().patch == (index_t) i)
                    {
                        if (mp_vertex.interfaces()[numInt].second().side().index() > 2 && dir == 0)
                            i_shift = numInt*5;
                        else if (mp_vertex.interfaces()[numInt].second().side().index() < 3 && dir == 1)
                            i_shift = numInt*5;
                    }
                }

                for (size_t ii = 0; ii < vertexAuxiliary_bf[i][dir].nPatches(); ii++) {
                    gsMatrix<> basisVal = vertexAuxiliary_bf[i][dir].patch(ii).eval(points);
                    gsMatrix<> basisDer = vertexAuxiliary_bf[i][dir].patch(ii).deriv(points);
                    gsMatrix<> basisDer2 = vertexAuxiliary_bf[i][dir].patch(ii).deriv2(points);
                    if (basisVal(0, 0) * basisVal(0, 0) > 1e-25)
                        coeffs_mat(i * 4 + 0, i_shift + ii) = (dir == 1 ? -1 : 1) * basisVal(0, 0);
                    if (basisDer(0, 0) * basisDer(0, 0) > 1e-25)
                        coeffs_mat(i * 4 + 1, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer(0, 0); // u
                    if (basisDer(1, 0) * basisDer(1, 0) > 1e-25)
                        coeffs_mat(i * 4 + 2, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer(1, 0); // v
                    if (basisDer2(2, 0) * basisDer2(2, 0) > 1e-25)
                        coeffs_mat(i * 4 + 3, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer2(2, 0); // uv
                }
            }
        }
        // Correction term
/*
        std::vector<size_t> interface_indices;
        index_t patchID_next = 0;
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            index_t patchID = patchID_next;
            index_t i_shift = 0, dir = 0;
            for(size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++)
            {
                bool intexist = false;
                for (size_t i_indizes = 0; i_indizes < interface_indices.size(); i_indizes++)
                    if (interface_indices[i_indizes] == numInt)
                        intexist = true;

                if (mp_vertex.interfaces()[numInt].first().patch == patchID && !intexist)
                {
                    i_shift = numInt;
                    dir = mp_vertex.interfaces()[numInt].first().side().index() > 2 ? 0 : 1;
                    patchID_next = mp_vertex.interfaces()[numInt].second().patch;
                    interface_indices.push_back(numInt);
                    break;
                }
                if (mp_vertex.interfaces()[numInt].second().patch == patchID  && !intexist)
                {
                    i_shift = numInt;
                    dir = mp_vertex.interfaces()[numInt].second().side().index() > 2 ? 0 : 1;
                    patchID_next = mp_vertex.interfaces()[numInt].first().patch;
                    interface_indices.push_back(numInt);
                    break;
                }
            }

            for (size_t ii = 0; ii < vertexAuxiliary_bf[patchID][dir].nPatches(); ii++) {
                gsMatrix<> basisVal = vertexAuxiliary_bf[patchID][dir].patch(ii).eval(points);
                gsMatrix<> basisDer = vertexAuxiliary_bf[patchID][dir].patch(ii).deriv(points);
                gsMatrix<> basisDer2 = vertexAuxiliary_bf[patchID][dir].patch(ii).deriv2(points);
                if (basisVal(0, 0) * basisVal(0, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 0, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisVal(0, 0);
                if (basisDer(0, 0) * basisDer(0, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 1, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer(0, 0); // u
                if (basisDer(1, 0) * basisDer(1, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 2, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer(1, 0); // v
                if (basisDer2(2, 0) * basisDer2(2, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 3, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer2(2, 0); // uv
            }

            gsBasis<> &basis_temp = vertexAuxiliary_bf[patchID][dir].patch(0).basis();
            for (size_t ii = 0; ii < 4; ii++) {
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 0) = (dir == 1 ? 1 : -1) * basis_temp.evalSingle(0,points)(0,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 1) = (dir == 1 ? 1 : -1) * basis_temp.derivSingle(1,points)(0,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 2) = (dir == 1 ? 1 : -1) * basis_temp.derivSingle(basis_temp.component(0).size(),points)(1,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 3) = (dir == 1 ? 1 : -1) * basis_temp.deriv2Single(basis_temp.component(0).size()+1,points)(2,0);
            }


        }
*/


        //gsInfo << "coeffs_mat " << coeffs_mat << "\n";
        real_t threshold = 1e-8;
        Eigen::FullPivLU<gsMatrix<>> KernelVertex(coeffs_mat);
        KernelVertex.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelVertex.dimensionOfKernel() < (index_t) mp_vertex.nInterfaces()+3)
        {
            threshold += 1e-8;
            KernelVertex.setThreshold(threshold);
        }
        gsMatrix<T> kernel = KernelVertex.kernel();

        //gsInfo << "Kernel: " << kernel << "\n";
        //gsInfo << "Kernel dim: " << KernelVertex.dimensionOfKernel() << "\n";

        basisVertexResult.resize(patchesAroundVertex.size());
        for(size_t i = 0; i < patchesAroundVertex.size(); i++) {
            for (index_t kernel_dim = 0; kernel_dim < kernel.cols(); kernel_dim++) {
                gsMatrix<> coef_bf, coef_bf_singleU, coef_bf_singleV;
                coef_bf.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_singleU.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_singleV.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);

                for (index_t dir = 0; dir < 2; dir++) {
                    index_t i_shift = 0;
                    for (size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++) {
                        if (mp_vertex.interfaces()[numInt].first().patch == (index_t) i) {
                            if (mp_vertex.interfaces()[numInt].first().side().index() > 2 && dir == 0)
                                i_shift = numInt * 5;
                            else if (mp_vertex.interfaces()[numInt].first().side().index() < 3 && dir == 1)
                                i_shift = numInt * 5;
                        }
                        if (mp_vertex.interfaces()[numInt].second().patch == (index_t) i) {
                            if (mp_vertex.interfaces()[numInt].second().side().index() > 2 && dir == 0)
                                i_shift = numInt * 5;
                            else if (mp_vertex.interfaces()[numInt].second().side().index() < 3 && dir == 1)
                                i_shift = numInt * 5;
                        }
                    }
                    for (size_t ii = 0; ii < vertexAuxiliary_bf[i][dir].nPatches(); ii++) {
                        coef_bf += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                        if (dir == 0)
                            coef_bf_singleU += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                        if (dir == 1)
                            coef_bf_singleV += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                    }
                }
                gsBasis<> &basis_temp = vertexAuxiliary_bf[i][0].patch(0).basis();
                gsGeometry<>::uPtr geo_tempU = basis_temp.makeGeometry(coef_bf_singleU);
                gsGeometry<>::uPtr geo_tempV = basis_temp.makeGeometry(coef_bf_singleV);
                gsMatrix<> coef_bf_corr;
                coef_bf_corr.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_corr(0,0) = coef_bf_singleU(0,0);
                coef_bf_corr(1,0) = coef_bf_singleU(1,0);
                coef_bf_corr(basis_temp.component(0).size(),0) = coef_bf_singleU(basis_temp.component(0).size(),0);
                coef_bf_corr(basis_temp.component(0).size()+1,0) = coef_bf_singleU(basis_temp.component(0).size()+1,0);

                //gsInfo << "coef_bf_corr: " << coef_bf_corr << "\n";
                //gsInfo << "coef_bf: " << coef_bf << "\n";
                coef_bf -= coef_bf_corr;

/*
                gsGeometry<>::uPtr geo_temp = basis_temp.makeGeometry(coef_bf);
                gsInfo << "Patch: " << i << "\n";
                gsInfo << "eval: " << geo_temp->eval(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_temp->deriv(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_temp->deriv(points)(1,0) << "\n";
                gsInfo << "deriv2: " << geo_temp->deriv2(points)(2,0) << "\n";

                coef_bf_corr.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_corr(0,0) = geo_tempU->eval(points)(0,0);
                coef_bf -= coef_bf_corr;

                gsInfo << "eval: " << geo_tempU->eval(points)(0,0) << " = " << geo_tempV->eval(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_tempU->deriv(points)(0,0) << " = " << geo_tempV->deriv(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_tempU->deriv(points)(1,0) << " = " << geo_tempV->deriv(points)(1,0) << "\n";
                gsInfo << "deriv2: " << geo_tempU->deriv2(points)(2,0) << " = " << geo_tempV->deriv2(points)(2,0) << "\n";

                coef_bf_corr.setZero();
                coef_bf_corr(1,0) = geo_tempU->deriv(points)(0,0) / basis_temp.derivSingle(1,points)(0,0);
                coef_bf -= coef_bf_corr;

                coef_bf_corr.setZero();
                coef_bf_corr(basis_temp.component(0).size(),0) = geo_tempU->deriv(points)(1,0) /
                        basis_temp.derivSingle(basis_temp.component(0).size(), points)(1,0);
                coef_bf -= coef_bf_corr;

                coef_bf_corr.setZero();
                coef_bf_corr(basis_temp.component(0).size()+1,0) = geo_tempU->deriv2(points)(2,0)/
                        basis_temp.deriv2Single(basis_temp.component(0).size()+1,points)(2,0);
                coef_bf -= coef_bf_corr;
*/
                basisVertexResult[i].addPatch(basis_temp.makeGeometry(coef_bf));
            }
        }

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
            m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "NoVerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
            {
                if (basisVertexResult.size() != 0)
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

/*    void saveBasisVertex(gsSparseMatrix<T> & system)
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
    }*/

    void reparametrizeVertexPatches()
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            if(m_auxPatches[i].getOrient() == 0) // not switched
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
            else if (m_auxPatches[i].getOrient() == 1) // switched
                switch (m_auxPatches[i].side()) // == vertex
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

    void reparametrizeVertexBasis(std::vector<std::vector<gsMultiPatch<T>>> & vertex_bf, std::vector<std::vector<gsMultiPatch<T>>> & result) // only for noVertex
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            size_t patch = m_patchesAroundVertex[i];
            size_t vertex = m_vertexIndices[i];
            index_t dir_u = -1, dir_v = -1;
            switch (vertex) // == vertex
            {
                case 1:
                    dir_u = 4;
                    dir_v = 0;
                    break;
                case 2:
                    dir_u = 5;
                    dir_v = 2;
                    break;
                case 3:
                    dir_u = 6;
                    dir_v = 1;
                    break;
                case 4:
                    dir_u = 7;
                    dir_v = 3;
                    break;
                default:
                    break;
            }

            gsMultiPatch<> basisVertex_u = vertex_bf[patch][dir_u];
            gsMultiPatch<> basisVertex_v = vertex_bf[patch][dir_v];

            if (m_auxPatches[i].getPatch().orientation() == -1)
            {
                m_auxPatches[i].swapBasisAxis(basisVertex_u);
                m_auxPatches[i].swapBasisAxis(basisVertex_v);
                gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";
            }

            if(m_auxPatches[i].getOrient() == 0) // not switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(2);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateBasisClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(-1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }
            else if (m_auxPatches[i].getOrient() == 1) // switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(2);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateBasisClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(-1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }

            result[i][0] = basisVertex_u;
            result[i][1] = basisVertex_v;

            //gsWriteParaview(basisVertex_u,"basisUVertexRot"+util::to_string(i),2000);
            //gsWriteParaview(basisVertex_v,"basisVVertexRot"+util::to_string(i),2000);
        }
    }

    void checkOrientation(size_t i)
    {
        if (m_auxPatches[i].getPatch().orientation() == -1)
        {
            m_auxPatches[i].swapAxis();
            //gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";
        }
    }

    real_t computeSigma(const std::vector<size_t> & vertexIndices)
    {
        real_t sigma = 0;

        real_t p = 0;
        real_t h_geo = 1;
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            gsTensorBSplineBasis<2, real_t> bsp_temp = dynamic_cast<gsTensorBSplineBasis<d, T>&>(m_auxPatches[0].getC1BasisRotated().getBasis(vertexIndices[i]+4));

            real_t p_temp = math::max(bsp_temp.degree(0), bsp_temp.degree(1));

            p = (p < p_temp ? p_temp : p);

            for(index_t j = 0; j < m_auxPatches[i].getPatch().parDim(); j++)
            {
                real_t h_geo_temp = bsp_temp.knot(j,p + 1);
                h_geo = (h_geo > h_geo_temp ? h_geo_temp : h_geo);
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
        if (matrix_det.determinant()*matrix_det.determinant() > 1e-15) // There is (numerically) a kink
            dofsCorner = 1;

        //for(size_t np = 0; np < mp_vertex.nPatches(); np++)
        //    m_bases[m_patchesAroundVertex[np]].setNumDofsVertex(dofsCorner, m_vertexIndices[np]);

        if (m_optionList.getSwitch("info"))
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

        real_t threshold = 1e-10;
        Eigen::FullPivLU<gsMatrix<>> KernelCorner(coefs_corner);
        KernelCorner.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelCorner.dimensionOfKernel() < dofsCorner)
        {
            threshold += 1e-8;
            KernelCorner.setThreshold(threshold);
        }
        if (m_optionList.getSwitch("info"))
            gsInfo << "Dimension of Kernel: " << KernelCorner.dimensionOfKernel() << " With " << threshold << "\n";

        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);

        gsMatrix<T> kernel = KernelCorner.kernel();

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
        if (m_optionList.getSwitch("info"))
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

    std::vector<gsMultiPatch<T>> getVertexBasis() { return basisVertexResult; }


protected:

    // Input
    gsMultiPatch<T> & m_mp;
    C1BasisContainer & m_bases;

    const std::vector<size_t> & m_patchesAroundVertex;
    const std::vector<size_t> & m_vertexIndices;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

}; // gsApproxC1Vertex

} // namespace gismo
