/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsC1Basis.h>
#include <gsUnstructuredSplines/gsG1AuxiliaryPatch.h>
#include <gsUnstructuredSplines/gsC1SurfBasisEdge.h>

#include <gsUnstructuredSplines/gsC1SurfGluingData.h>

namespace gismo
{
template<short_t d, class T>
class gsC1SurfEdge
{

private:
    typedef gsC1Basis<d, T> Basis;
    typedef typename std::vector<Basis> C1BasisContainer;
    typedef typename std::vector<gsC1AuxiliaryPatch<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsC1SurfEdge
    typedef memory::shared_ptr<gsC1SurfEdge> Ptr;

    /// Unique pointer for gsC1SurfEdge
    typedef memory::unique_ptr<gsC1SurfEdge> uPtr;


public:
    /// Empty constructor
    ~gsC1SurfEdge() { }

    gsC1SurfEdge(const gsMultiPatch<> & mp, const boundaryInterface & item){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), item.first().patch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), item.second().patch));
    }

    gsC1SurfEdge(const gsMultiPatch<> & sp, const patchSide & item){
        auxGeom.push_back(gsG1AuxiliaryPatch(sp.patch(item.patch), item.patch));
    }

    void computeG1InterfaceBasis(gsG1OptionList g1OptionList)
    {
        gsMultiPatch<> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());// Right -> 0 ====> v along the interface
        mp_init.addPatch(auxGeom[1].getPatch()); // Left -> 1 ====> u along the interface

        gsMultiPatch<> test_mp(this->reparametrizeG1Interface()); // auxGeom contains now the reparametrized geometry
        gsMultiBasis<> test_mb(test_mp);
        gsMultiPatch<> g1Basis_0, g1Basis_1;

        gsG1ASGluingData<real_t> g1BasisEdge(test_mp, test_mp);
        gsG1ASBasisEdge<real_t> g1BasisEdge_0(test_mp.patch(0), test_mb.basis(0), 1, false, g1OptionList, g1BasisEdge);
        gsG1ASBasisEdge<real_t> g1BasisEdge_1(test_mp.patch(1), test_mb.basis(1), 0, false, g1OptionList, g1BasisEdge);
        g1BasisEdge_0.setG1BasisEdge(g1Basis_0);
        g1BasisEdge_1.setG1BasisEdge(g1Basis_1);

//      Patch 0 -> Right
        auxGeom[0].parametrizeBasisBack(g1Basis_0);
//      Patch 1 -> Left
        auxGeom[1].parametrizeBasisBack(g1Basis_1);
    }

    void computeG1BoundaryBasis(gsG1OptionList g1OptionList, const int boundaryInd)
    {
        gsMultiPatch<> test_mp(this->reparametrizeG1Boundary(boundaryInd));
        gsMultiBasis<> test_mb(test_mp);
        gsMultiPatch<> g1Basis_edge;

        gsG1ASGluingData<real_t> bdyGD; // Empty constructor creates the sol and solBeta in a suitable way to manage the GD on the boundary
        gsG1ASBasisEdge<real_t> g1BasisEdge(test_mp, test_mb, 1, true, g1OptionList, bdyGD);
        g1BasisEdge.setG1BasisEdge(g1Basis_edge);

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);
    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i)
    {
        return auxGeom[i];
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
        for (index_t i = m_bases[patch_1].rowBegin(side_1); i < m_bases[patch_1].rowEnd(side_1); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_1].colBegin(side_1);
                 j < m_bases[patch_1].colEnd(side_1); ++j, ++jj) {
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

        ii = 0;
        for (index_t i = m_bases[patch_2].rowBegin(side_2); i < m_bases[patch_2].rowEnd(side_2); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_2].colBegin(side_2);
                 j < m_bases[patch_2].colEnd(side_2); ++j, ++jj)
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
        for (index_t i = m_bases[patch_1].rowBegin(side_1); i < m_bases[patch_1].rowEnd(side_1); ++i, ++ii)
        {
            index_t jj = 0;
            for (index_t j = m_bases[patch_1].colBegin(side_1);
                 j < m_bases[patch_1].colEnd(side_1); ++j, ++jj)
                if (basisEdgeResult[0].patch(ii).coef(jj, 0) * basisEdgeResult[0].patch(ii).coef(jj, 0) > 1e-25)
                    system.insert(shift_row + i, shift_col + j) = basisEdgeResult[0].patch(ii).coef(jj, 0);
        }

    }

protected:

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

    std::vector<gsG1AuxiliaryPatch> auxGeom;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    gsMultiPatch<> computeAuxTopology();

    gsMultiPatch<> reparametrizeInterface();

    gsMultiPatch<> reparametrizeBoundary(index_t side);


}; // Class gsC1SurfEdge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SurfEdge.hpp)
#endif