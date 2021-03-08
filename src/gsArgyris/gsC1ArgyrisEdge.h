/** @file gsC1Argyris.h

    @brief Creates the C1 Argyris space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsArgyris/gsC1ArgyrisBasis.h>
#include <gsArgyris/gsC1ArgyrisAuxiliaryPatch.h>

#include <gsArgyris/gsGluingData/gsApproxGluingData.h>
#include <gsArgyris/gsApproxArgyrisEdgeBasis.h>

namespace gismo
{
template<short_t d, class T>
class gsC1ArgyrisEdge
{

private:
    typedef gsC1ArgyrisBasis<d, T> Basis;
    typedef typename std::vector<Basis> ArgyrisBasisContainer;
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

    /// Shared pointer for gsC1Argyris
    typedef memory::shared_ptr<gsC1ArgyrisEdge> Ptr;

    /// Unique pointer for gsC1Argyris
    typedef memory::unique_ptr<gsC1ArgyrisEdge> uPtr;


public:
    /// Empty constructor
    ~gsC1ArgyrisEdge() { }


    gsC1ArgyrisEdge(gsMultiPatch<T> const & mp,
                ArgyrisBasisContainer & bases,
                const boundaryInterface & item,
                size_t & numInt,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        side_1 = item.first().side().index();
        side_2 = item.second().side().index();

        patch_1 = item.first().patch;
        patch_2 = item.second().patch;

        //const index_t dir_1 = side_1 > 2 ? 0 : 1;
        //const index_t dir_2 = side_2 > 2 ? 0 : 1;

        m_auxPatches.clear();
        m_auxPatches.push_back(gsC1ArgyrisAuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], side_1));
        m_auxPatches.push_back(gsC1ArgyrisAuxiliaryPatch<d,T>(m_mp.patch(patch_2), m_bases[patch_2], side_2));

        reparametrizeInterfacePatches();

        // Compute GLuing data
        gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList);
        gsApproxArgyrisEdgeBasis<d, T> approxArgyrisEdgeBasis(m_auxPatches, approxGluingData, 0, m_optionList);
        gsApproxArgyrisEdgeBasis<d, T> approxArgyrisEdgeBasis2(m_auxPatches, approxGluingData, 1, m_optionList);
        gsMultiPatch<T> result_1, result_2;
        approxArgyrisEdgeBasis.setG1BasisEdge(result_1);
        approxArgyrisEdgeBasis2.setG1BasisEdge(result_2);

        // Compute Kernel (before parametrizeBack)
        computeKernel(result_1, result_2, side_1);

        // parametrizeBasisBack
        m_auxPatches[0].parametrizeBasisBack(result_1);
        m_auxPatches[1].parametrizeBasisBack(result_2);

        basisEdgeResult.clear();
        basisEdgeResult.push_back(result_1);
        basisEdgeResult.push_back(result_2);
/*
        std::string fileName;
        std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
        gsParaviewCollection collection(basename);

        for (size_t i = 0; i< result_1.nPatches(); i++)
        {
            // First Interface Side
            fileName = basename + "_0_" + util::to_string(i);
            gsField<> temp_field(m_mp.patch(patch_1), result_1.patch(i));
            gsWriteParaview(temp_field, fileName, 5000);
            collection.addTimestep(fileName, i, "0.vts");
            // Second Interface Side
            fileName = basename + "_1_" + util::to_string(i);
            gsField<> temp_field_1(m_mp.patch(patch_2), result_2.patch(i));
            gsWriteParaview(temp_field_1, fileName, 5000);
            collection.addTimestep(fileName, i, "0.vts");
        }
        collection.save();
*/
    }

    gsC1ArgyrisEdge(gsMultiPatch<T> const & mp,
                ArgyrisBasisContainer & bases,
                const patchSide & item,
                size_t & numBdy,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        side_1 = item.side().index();
        patch_1 = item.patch;

        //const index_t dir_1 = side_1 > 2 ? 0 : 1;

        m_auxPatches.clear();
        m_auxPatches.push_back(gsC1ArgyrisAuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], side_1));

        reparametrizeSinglePatch(side_1);

        // Compute GLuing data NO NEED
        // gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList);
        gsApproxArgyrisEdgeBasis<d, T> approxArgyrisEdgeBasis(m_auxPatches, 0, m_optionList);
        gsMultiPatch<T> result_1;
        approxArgyrisEdgeBasis.setG1BasisEdge(result_1);

        // Compute Kernel (before parametrizeBack) NO NEED
        //computeKernel(result_1, result_2, side_1);

        // parametrizeBasisBack
        m_auxPatches[0].parametrizeBasisBack(result_1);

        basisEdgeResult.clear();
        basisEdgeResult.push_back(result_1);
/*
        std::string fileName;
        std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
        gsParaviewCollection collection(basename);

        for (size_t i = 0; i< result_1.nPatches(); i++)
        {
            // First Interface Side
            fileName = basename + "_0_" + util::to_string(i);
            gsField<> temp_field(m_mp.patch(patch_1), result_1.patch(i));
            gsWriteParaview(temp_field, fileName, 5000);
            collection.addTimestep(fileName, i, "0.vts");
        }
        collection.save();
*/
    }

    void saveBasisInterface(gsSparseMatrix<T> & system)
    {

        index_t shift_row = 0, shift_col = 0;
        for (index_t np = 0; np < patch_1; ++np)
        {
            shift_row += m_bases[np].size_rows();
            shift_col += m_bases[np].size_cols();
        }

        index_t ii = 0;
        for (index_t i = m_bases[patch_1].rowBegin("edge",side_1); i < m_bases[patch_1].rowEnd("edge",side_1); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_1].colBegin("edge", side_1);
                 j < m_bases[patch_1].colEnd("edge", side_1); ++j, ++jj) {
                if (basisEdgeResult[0].patch(ii).coef(jj, 0) * basisEdgeResult[0].patch(ii).coef(jj, 0) > 1e-25)
                    system.insert(shift_row + i, shift_col + j) = basisEdgeResult[0].patch(ii).coef(jj, 0);
            }
        }

        shift_row = 0;
        shift_col = 0;
        for (index_t np = 0; np < patch_2; ++np)
        {
            shift_row += m_bases[np].size_rows();
            shift_col += m_bases[np].size_cols();
        }

        gsInfo << "i: " << shift_row << "\n";
        gsInfo << "j: " <<  shift_col << "\n";

        ii = 0;
        for (index_t i = m_bases[patch_2].rowBegin("edge",side_2); i < m_bases[patch_2].rowEnd("edge",side_2); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_2].colBegin("edge", side_2);
                 j < m_bases[patch_2].colEnd("edge", side_2); ++j, ++jj)
                if (basisEdgeResult[1].patch(ii).coef(jj, 0) * basisEdgeResult[1].patch(ii).coef(jj, 0) > 1e-25)
                    system.insert(shift_row + i, shift_col + j) = basisEdgeResult[1].patch(ii).coef(jj, 0);
        }
    }

    void saveBasisBoundary(gsSparseMatrix<T> & system)
    {
        index_t shift_row = 0, shift_col = 0;
        for (index_t np = 0; np < patch_1; ++np)
        {
            shift_row += m_bases[np].size_rows();
            shift_col += m_bases[np].size_cols();
        }

        index_t ii = 0;
        for (index_t i = m_bases[patch_1].rowBegin("edge",side_1); i < m_bases[patch_1].rowEnd("edge",side_1); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_1].colBegin("edge", side_1);
                 j < m_bases[patch_1].colEnd("edge", side_1); ++j, ++jj)
                if (basisEdgeResult[0].patch(ii).coef(jj, 0) * basisEdgeResult[0].patch(ii).coef(jj, 0) > 1e-25)
                    system.insert(shift_row + i, shift_col + j) = basisEdgeResult[0].patch(ii).coef(jj, 0);
        }
    }

protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    ArgyrisBasisContainer & m_bases;

    const gsOptionList & m_optionList;

    index_t patch_1, patch_2, side_1, side_2;

    // Need for rotation, etc.
    ArgyrisAuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    void computeAuxTopology();

    void reparametrizeInterfacePatches();

    void reparametrizeSinglePatch(index_t side);

    void computeKernel(gsMultiPatch<> & result_0, gsMultiPatch<> & result_1, index_t side_0);

}; // Class gsC1ArgyrisEdge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1ArgyrisEdge.hpp)
#endif