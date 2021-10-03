/** @file gsApproxC1Edge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsC1Basis.h>
#include <gsUnstructuredSplines/gsC1AuxiliaryPatch.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>
#include <gsUnstructuredSplines/gsApproxC1EdgeBasisProjection.h>

namespace gismo
{
template<short_t d, class T>
class gsApproxC1Edge
{

private:
    typedef gsC1Basis<d, T> Basis;
    typedef typename std::vector<Basis> C1BasisContainer;
    typedef typename std::vector<gsC1AuxiliaryPatch<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Edge
    typedef memory::shared_ptr<gsApproxC1Edge> Ptr;

    /// Unique pointer for gsApproxC1Edge
    typedef memory::unique_ptr<gsApproxC1Edge> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Edge() { }


    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                C1BasisContainer & bases,
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
        m_auxPatches.push_back(gsC1AuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], side_1));
        m_auxPatches.push_back(gsC1AuxiliaryPatch<d,T>(m_mp.patch(patch_2), m_bases[patch_2], side_2));

        reparametrizeInterfacePatches();

        // Compute GLuing data
        gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList);

        gsMultiPatch<T> result_1, result_2;

        gsApproxC1EdgeBasisProjection<d, T> approxEdgeBasis(m_auxPatches, approxGluingData, 0, m_optionList);
        gsApproxC1EdgeBasisProjection<d, T> approxEdgeBasis2(m_auxPatches, approxGluingData, 1, m_optionList);

        approxEdgeBasis.setG1BasisEdge(result_1);
        approxEdgeBasis2.setG1BasisEdge(result_2);


        // parametrizeBasisBack
        m_auxPatches[0].parametrizeBasisBack(result_1);
        m_auxPatches[1].parametrizeBasisBack(result_2);

        basisEdgeResult.clear();
        basisEdgeResult.push_back(result_1);
        basisEdgeResult.push_back(result_2);

        if (m_optionList.getSwitch("plot"))
        {
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
        }
    }

    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                C1BasisContainer & bases,
                const patchSide & item,
                size_t & numBdy,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        side_1 = item.side().index();
        patch_1 = item.patch;

        //const index_t dir_1 = side_1 > 2 ? 0 : 1;

        m_auxPatches.clear();
        m_auxPatches.push_back(gsC1AuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], side_1));

        reparametrizeSinglePatch(side_1);

        gsMultiPatch<> result_1;
        gsApproxC1EdgeBasisProjection<d, T> approxEdgeBasis(m_auxPatches, 0, m_optionList);
        approxEdgeBasis.setG1BasisEdge(result_1);


        // parametrizeBasisBack
        m_auxPatches[0].parametrizeBasisBack(result_1);

        basisEdgeResult.clear();
        basisEdgeResult.push_back(result_1);

        if (m_optionList.getSwitch("plot")) {
            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i < result_1.nPatches(); i++) {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(patch_1), result_1.patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();
        }
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

    void saveBasisVertex(std::vector<std::vector<gsMultiPatch<T>>> & vertex_bf)
    {

        for (index_t i = 0; i < 2; i++) {
            index_t side = i == 0 ? side_1 : side_2;
            index_t patch = i == 0 ? patch_1 : patch_2;

            gsMultiPatch<T> basis_1, basis_2;

            index_t size_plus = m_bases[patch].getBasisPlus(side).size();
            index_t size_minus = m_bases[patch].getBasisMinus(side).size();

            basis_1.addPatch(basisEdgeResult[i].patch(0));
            basis_1.addPatch(basisEdgeResult[i].patch(1));
            basis_1.addPatch(basisEdgeResult[i].patch(2));
            basis_1.addPatch(basisEdgeResult[i].patch(size_plus));
            basis_1.addPatch(basisEdgeResult[i].patch(size_plus + 1));

            basis_2.addPatch(basisEdgeResult[i].patch(size_plus - 1));
            basis_2.addPatch(basisEdgeResult[i].patch(size_plus - 2));
            basis_2.addPatch(basisEdgeResult[i].patch(size_plus - 3));
            basis_2.addPatch(basisEdgeResult[i].patch(size_plus + size_minus - 1));
            basis_2.addPatch(basisEdgeResult[i].patch(size_plus + size_minus - 2));
/*
            for (index_t ii = 0; ii < 5; ii++)
            {
                vertex_bf[patch][(side-1)].addPatch(basis_1.patch(ii)); // -1 bcs of c++ counting
                vertex_bf[patch][(side-1)].addPatch(basis_2.patch(ii));
            }
*/


            for (index_t ii = 0; ii < 5; ii++)
                if (i == 1)
                    switch (side) {
                        case 1:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_2.patch(ii)); // -1 bcs of c++ counting // vertex 1
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_1.patch(ii)); // vertex 3
                            break;
                        case 2:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_1.patch(ii)); // vertex 2
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_2.patch(ii)); // vertex 4
                            break;
                        case 3:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_1.patch(ii)); // vertex 1
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_2.patch(ii)); // vertex 2
                            break;
                        case 4:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_2.patch(ii)); // vertex 3
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_1.patch(ii)); // vertex 4
                            break;
                        default:
                            gsInfo << "Wrong side index\n";
                            break;
                    }
                else if (i == 0)
                    switch (side) {
                        case 1:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_1.patch(ii)); // -1 bcs of c++ counting // vertex 1
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_2.patch(ii)); // vertex 3
                            break;
                        case 2:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_2.patch(ii)); // vertex 2
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_1.patch(ii)); // vertex 4
                            break;
                        case 3:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_2.patch(ii)); // vertex 1
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_1.patch(ii)); // vertex 2
                            break;
                        case 4:
                            vertex_bf[patch][(side-1)*2 + 0].addPatch(basis_1.patch(ii)); // vertex 3
                            vertex_bf[patch][(side-1)*2 + 1].addPatch(basis_2.patch(ii)); // vertex 4
                            break;
                        default:
                            gsInfo << "Wrong side index\n";
                            break;
                    }

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

    void interpolateBasisInterface(gsApproxGluingData<d, T> & approxGluingData, gsMultiPatch<> & result_1, gsMultiPatch<> & result_2)
    {
        for (index_t patchID = 0; patchID < 2; ++patchID)
        {
            index_t dir = patchID == 0 ? 1 : 0;
            index_t side = m_auxPatches[patchID].side();

            gsTensorBSplineBasis<d, T> basis_edge = m_auxPatches[patchID].getC1BasisRotated().getEdgeBasis(side); // 0 -> u, 1 -> v

            gsBSplineBasis<T> basis_plus = m_auxPatches[patchID].getC1BasisRotated().getBasisPlus(side);
            gsBSplineBasis<T> basis_minus = m_auxPatches[patchID].getC1BasisRotated().getBasisMinus(side);
            gsBSplineBasis<T> basis_geo = m_auxPatches[patchID].getC1BasisRotated().getBasisGeo(side);

            index_t n_plus = basis_plus.size();
            index_t n_minus = basis_minus.size();

            // tau/p
            real_t p = basis_geo.degree();
            real_t tau_1 = basis_geo.knots().at(p + 1); // p + 2


            index_t bfID_init = 3;
            for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
            {
                // Points to interpolate at (Greville points):
                gsMatrix<> points = basis_edge.anchors();

                // Evaluate f at the Greville points
                gsMatrix<T> beta, N_0, N_1, N_i_plus, der_N_i_plus;

                approxGluingData.betaS(dir).eval_into(points.row(dir),beta); // 1-dir == PatchID

                basis_geo.evalSingle_into(0,points.row(1-dir),N_0); // u
                basis_geo.evalSingle_into(1,points.row(1-dir),N_1); // u

                basis_plus.evalSingle_into(bfID,points.row(dir),N_i_plus); // v
                basis_plus.derivSingle_into(bfID,points.row(dir),der_N_i_plus);

                gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);
                gsMatrix<T> fValues = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;

                // Returns a geometry with basis = tBasis
                // and coefficients being
                // computed as the interpolant of \a funct
                gsGeometry<>::uPtr interpolant = basis_edge.interpolateAtAnchors(fValues);
                if (patchID == 0)
                    result_1.addPatch(*interpolant);
                else
                    result_2.addPatch(*interpolant);
            }

            bfID_init = 2;
            for (index_t bfID = bfID_init; bfID < n_minus-bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
            {
                // Points to interpolate at (Greville points):
                gsMatrix<T> points = basis_edge.anchors();

                // Evaluate f at the Greville points
                // Evaluate f at the Greville points
                gsMatrix<T> alpha, N_1, N_j_minus;

                approxGluingData.alphaS(dir).eval_into(points.row(dir),alpha); // 1-dir == PatchID

                basis_minus.evalSingle_into(bfID,points.row(dir),N_j_minus); // v
                basis_geo.evalSingle_into(1,points.row(1-dir),N_1); // u

                gsMatrix<T> fValues = (dir == 0 ? -1 : 1) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;

                // Returns a geometry with basis = tBasis
                // and coefficients being
                // computed as the interpolant of \a funct
                gsGeometry<>::uPtr interpolant = basis_edge.interpolateAtAnchors(fValues);
                if (patchID == 0)
                    result_1.addPatch(*interpolant);
                else
                    result_2.addPatch(*interpolant);
            }
        }
    }

    void interpolateBasisBoundary(gsMultiPatch<> & result_1)
    {
        index_t side = m_auxPatches[0].side();
        index_t dir = 1;

        gsTensorBSplineBasis<d, T> basis_edge = m_auxPatches[0].getC1BasisRotated().getEdgeBasis(m_auxPatches[0].side()); // 0 -> u, 1 -> v

        gsBSplineBasis<T> basis_plus = m_auxPatches[0].getC1BasisRotated().getBasisPlus(side);
        gsBSplineBasis<T> basis_minus = m_auxPatches[0].getC1BasisRotated().getBasisMinus(side);
        gsBSplineBasis<T> basis_geo = m_auxPatches[0].getC1BasisRotated().getBasisGeo(side);

        index_t n_plus = basis_plus.size();
        index_t n_minus = basis_minus.size();

        index_t bfID_init = 3;
        for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
        {
            // Points to interpolate at (Greville points):
            gsMatrix<> points = basis_edge.anchors();

            // Evaluate f at the Greville points
            gsMatrix<T> N_0, N_1, N_i_plus;

            basis_geo.evalSingle_into(0,points.row(1-dir),N_0); // u
            basis_geo.evalSingle_into(1,points.row(1-dir),N_1); // u
            basis_plus.evalSingle_into(bfID,points.row(dir),N_i_plus); // v

            gsMatrix<> fValues = N_i_plus.cwiseProduct(N_0 + N_1);

            // Returns a geometry with basis = tBasis
            // and coefficients being
            // computed as the interpolant of \a funct
            gsGeometry<>::uPtr interpolant = basis_edge.interpolateAtAnchors(fValues);
            result_1.addPatch(*interpolant);
        }

        bfID_init = 2;
        for (index_t bfID = bfID_init; bfID < n_minus-bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
        {
            // Points to interpolate at (Greville points):
            gsMatrix<> points = basis_edge.anchors();

            // Evaluate f at the Greville points
            // Evaluate f at the Greville points
            gsMatrix<T> N_0, N_1, N_j_minus;

            basis_minus.evalSingle_into(bfID,points.row(dir),N_j_minus); // v
            basis_geo.evalSingle_into(0,points.row(1-dir),N_0); // u
            basis_geo.evalSingle_into(1,points.row(1-dir),N_1); // u

            gsMatrix<> fValues = (dir == 0 ? -1 : 1) * N_j_minus.cwiseProduct(N_1);

            // Returns a geometry with basis = tBasis
            // and coefficients being
            // computed as the interpolant of \a funct
            gsGeometry<>::uPtr interpolant = basis_edge.interpolateAtAnchors(fValues);
            result_1.addPatch(*interpolant);
        }
    }

protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    C1BasisContainer & m_bases;

    const gsOptionList & m_optionList;

    index_t patch_1, patch_2, side_1, side_2;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

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

}; // Class gsApproxC1Edge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Edge.hpp)
#endif