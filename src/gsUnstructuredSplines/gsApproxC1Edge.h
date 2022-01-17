/** @file gsApproxC1Edge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsContainerBasis.h>
#include <gsUnstructuredSplines/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>

namespace gismo
{


template<short_t d, class T>
class gsApproxC1Edge
{

private:
    typedef gsContainerBasis<d, T> Basis;
    typedef typename std::vector<Basis> BasisContainer;
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Edge
    typedef memory::shared_ptr<gsApproxC1Edge> Ptr;

    /// Unique pointer for gsApproxC1Edge
    typedef memory::unique_ptr<gsApproxC1Edge> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Edge() { }


    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const boundaryInterface & item,
                size_t & numInt,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        index_t side_1 = item.first().side().index();
        index_t side_2 = item.second().side().index();

        index_t patch_1 = item.first().patch;
        index_t patch_2 = item.second().patch;

        //const index_t dir_1 = side_1 > 2 ? 0 : 1;
        //const index_t dir_2 = side_2 > 2 ? 0 : 1;

        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1]));
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_2), m_bases[patch_2]));

        std::vector<index_t> sidesContainer(2);
        sidesContainer[0] = side_1;
        sidesContainer[1] = side_2;

        reparametrizeInterfacePatches();

        // Compute GLuing data
        gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList, sidesContainer);

        //! [Problem setup]
        basisEdgeResult.clear();
        for (index_t patchID = 0; patchID < 2; patchID++) {
            gsMultiPatch<T> result;

            index_t dir = patchID == 0 ? 1 : 0;

            gsBSplineBasis<T> basis_plus = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                    sidesContainer[patchID] - 1, 0));
            gsBSplineBasis<T> basis_minus = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                    sidesContainer[patchID] - 1, 1));
            gsBSplineBasis<T> basis_geo = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                    sidesContainer[patchID] - 1, 2));
            gsGeometry<T> &geo = m_auxPatches[patchID].getPatchRotated();

            gsBSpline<T> beta = approxGluingData.betaS(dir);
            gsBSpline<T> alpha = approxGluingData.alphaS(dir);

            index_t n_plus = basis_plus.size();
            index_t n_minus = basis_minus.size();

            index_t bfID_init = 3;
            for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
            {
                gsSparseSolver<real_t>::LU solver;
                gsExprAssembler<> A(1, 1);

                typedef gsExprAssembler<>::variable variable;
                typedef gsExprAssembler<>::space space;
                typedef gsExprAssembler<>::solution solution;

                // Elements used for numerical integration
                gsMultiBasis<T> edgeSpace(
                        m_auxPatches[patchID].getBasisRotated().getBasis(sidesContainer[patchID]));
                A.setIntegrationElements(edgeSpace);
                gsExprEvaluator<> ev(A);

                // Set the discretization space
                space u = A.getSpace(edgeSpace);

                // Create Mapper
                gsDofMapper map(edgeSpace);
                gsMatrix<index_t> act;
                for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
                map.finalize();

                gsBoundaryConditions<> bc_empty;
                bc_empty.addCondition(dir == 0 ? 1 : 3, condition_type::dirichlet, 0); // Doesn't matter which side
                u.setupMapper(map);
                A.initSystem();

                gsTraceBasis<real_t> traceBasis(geo, basis_plus, basis_geo, beta, false, bfID, dir);
                auto aa = A.getCoeff(traceBasis);

                A.assemble(u * u.tr(), u * aa);

                solver.compute(A.matrix());
                gsMatrix<> solVector = solver.solve(A.rhs());

                solution u_sol = A.getSolution(u, solVector);
                gsMatrix<> sol;
                u_sol.extract(sol);

                result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
            }

            bfID_init = 2;
            for (index_t bfID = bfID_init; bfID < n_minus - bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
            {
                gsSparseSolver<real_t>::LU solver;
                gsExprAssembler<> A(1, 1);

                typedef gsExprAssembler<>::variable variable;
                typedef gsExprAssembler<>::space space;
                typedef gsExprAssembler<>::solution solution;

                // Elements used for numerical integration
                gsMultiBasis<T> edgeSpace(
                        m_auxPatches[patchID].getBasisRotated().getBasis(sidesContainer[patchID]));
                A.setIntegrationElements(edgeSpace);
                gsExprEvaluator<> ev(A);

                // Set the discretization space
                space u = A.getSpace(edgeSpace);

                // Create Mapper
                gsDofMapper map(edgeSpace);
                gsMatrix<index_t> act;
                for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
                map.finalize();

                gsBoundaryConditions<> bc_empty;
                bc_empty.addCondition(dir == 0 ? 1 : 3, condition_type::dirichlet, 0); // Doesn't matter which side
                u.setupMapper(map);
                A.initSystem();

                gsNormalDerivBasis<real_t> normalDerivBasis(geo, basis_minus, basis_geo, alpha, false, bfID, dir);
                auto aa = A.getCoeff(normalDerivBasis);

                A.assemble(u * u.tr(), u * aa);

                solver.compute(A.matrix());
                gsMatrix<> solVector = solver.solve(A.rhs());

                solution u_sol = A.getSolution(u, solVector);
                gsMatrix<> sol;
                u_sol.extract(sol);

                result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
            }

            // parametrizeBasisBack
            m_auxPatches[patchID].parametrizeBasisBack(result);

            basisEdgeResult.push_back(result);
        }

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i< basisEdgeResult[0].nPatches(); i++)
            {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(patch_1), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
                // Second Interface Side
                fileName = basename + "_1_" + util::to_string(i);
                gsField<> temp_field_1(m_mp.patch(patch_2), basisEdgeResult[1].patch(i));
                gsWriteParaview(temp_field_1, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();
        }
    }

    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const patchSide & item,
                size_t & numBdy,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        index_t side_1 = item.side().index();
        index_t patch_1 = item.patch;

        //const index_t dir_1 = side_1 > 2 ? 0 : 1;

        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1]));

        reparametrizeSinglePatch(side_1);

        basisEdgeResult.clear();
        gsMultiPatch<T> result;

        index_t patchID = 0;
        index_t dir = patchID == 0 ? 1 : 0;

        gsBSplineBasis<T> basis_plus = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                side_1 - 1, 0));
        gsBSplineBasis<T> basis_minus = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                side_1 - 1, 1));
        gsBSplineBasis<T> basis_geo = dynamic_cast<gsBSplineBasis<T> &>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(
                side_1 - 1, 2));
        gsGeometry<T> &geo = m_auxPatches[patchID].getPatchRotated();

        gsBSpline<T> beta, alpha;

        index_t n_plus = basis_plus.size();
        index_t n_minus = basis_minus.size();

        index_t bfID_init = 3;
        for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
        {
            gsSparseSolver<real_t>::SimplicialLDLT solver;
            gsExprAssembler<> A(1, 1);

            typedef gsExprAssembler<>::variable variable;
            typedef gsExprAssembler<>::space space;
            typedef gsExprAssembler<>::solution solution;

            // Elements used for numerical integration
            gsMultiBasis<T> edgeSpace(
                    m_auxPatches[patchID].getBasisRotated().getBasis(side_1));
            A.setIntegrationElements(edgeSpace);
            gsExprEvaluator<> ev(A);

            // Set the discretization space
            space u = A.getSpace(edgeSpace);

            // Create Mapper
            gsDofMapper map(edgeSpace);
            gsMatrix<index_t> act;
            for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                 i++) // only the first two u/v-columns are Dofs (0/1)
            {
                act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                map.markBoundary(0, act); // Patch 0
            }
            map.finalize();

            gsBoundaryConditions<> bc_empty;
            bc_empty.addCondition(dir == 0 ? 1 : 3, condition_type::dirichlet, 0); // Doesn't matter which side
            u.setupMapper(map);
            A.initSystem();

            gsTraceBasis<real_t> traceBasis(geo, basis_plus, basis_geo, beta, true, bfID, dir);
            auto aa = A.getCoeff(traceBasis);

            A.assemble(u * u.tr(), u * aa);

            solver.compute(A.matrix());
            gsMatrix<> solVector = solver.solve(A.rhs());

            solution u_sol = A.getSolution(u, solVector);
            gsMatrix<> sol;
            u_sol.extract(sol);

            result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));

            //gsDebugVar(sol-result_1.patch(bfID-3).coefs());
        }

        bfID_init = 2;
        for (index_t bfID = bfID_init; bfID < n_minus - bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
        {
            gsSparseSolver<real_t>::SimplicialLDLT solver;
            gsExprAssembler<> A(1, 1);

            typedef gsExprAssembler<>::variable variable;
            typedef gsExprAssembler<>::space space;
            typedef gsExprAssembler<>::solution solution;

            // Elements used for numerical integration
            gsMultiBasis<T> edgeSpace(
                    m_auxPatches[patchID].getBasisRotated().getBasis(side_1));
            A.setIntegrationElements(edgeSpace);
            gsExprEvaluator<> ev(A);

            // Set the discretization space
            space u = A.getSpace(edgeSpace);

            // Create Mapper
            gsDofMapper map(edgeSpace);
            gsMatrix<index_t> act;
            for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                map.markBoundary(0, act); // Patch 0
            }
            map.finalize();

            gsBoundaryConditions<> bc_empty;
            bc_empty.addCondition(dir == 0 ? 1 : 3, condition_type::dirichlet, 0); // Doesn't matter which side
            u.setupMapper(map);
            A.initSystem();

            gsNormalDerivBasis<real_t> normalDerivBasis(geo, basis_minus, basis_geo, alpha, true, bfID, dir);
            auto aa = A.getCoeff(normalDerivBasis);

            A.assemble(u * u.tr(), u * aa);

            solver.compute(A.matrix());
            gsMatrix<> solVector = solver.solve(A.rhs());

            solution u_sol = A.getSolution(u, solVector);
            gsMatrix<> sol;
            u_sol.extract(sol);

            result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
        }

        // parametrizeBasisBack
        m_auxPatches[patchID].parametrizeBasisBack(result);
        basisEdgeResult.push_back(result);

        if (m_optionList.getSwitch("plot")) {
            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i < basisEdgeResult[0].nPatches(); i++) {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(patch_1), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();
        }
    }

    std::vector<gsMultiPatch<T>> getEdgeBasis() { return basisEdgeResult; };

protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    BasisContainer & m_bases;

    const gsOptionList & m_optionList;

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

}; // Class gsApproxC1Edge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Edge.hpp)
#endif
