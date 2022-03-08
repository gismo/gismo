/** @file gsApproxC1Edge.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines2/gsApproxC1Edge.h>

#include <gsUnstructuredSplines2/gsPatchReparameterized.h>

#include <gsUnstructuredSplines2/gsApproxGluingData.h>

namespace gismo
{
    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::compute(std::vector<patchSide> & sidesContainer) {

        // Compute GLuing data
        gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList, sidesContainer);

        //! [Problem setup]
        basisEdgeResult.clear();
        for (size_t patchID = 0; patchID < sidesContainer.size(); patchID++) {
            gsMultiPatch<T> result;

            index_t dir = patchID == 0 ? 1 : 0;

            gsBSplineBasis<T> basis_plus, basis_minus;

            gsMultiBasis<T> initSpace(m_auxPatches[patchID].getBasisRotated().piece(0));
            createPlusSpace(m_auxPatches[0].getPatchRotated(), initSpace.basis(0), dir, basis_plus);
            createMinusSpace(m_auxPatches[0].getPatchRotated(), initSpace.basis(0), dir, basis_minus);

            gsGeometry<T> &geo = m_auxPatches[patchID].getPatchRotated();

            gsBSpline<T> beta, alpha;
            bool bdy = true;
            if (sidesContainer.size() == 2)
            {
                bdy = false;
                beta = approxGluingData.betaS(dir);
                alpha = approxGluingData.alphaS(dir);
            }

            // [!The same setup for each bf!]
            gsSparseSolver<real_t>::SimplicialLDLT solver;
            gsExprAssembler<T> A(1, 1);

            // Elements used for numerical integration
            gsMultiBasis<T> edgeSpace(
                    m_auxPatches[patchID].getBasisRotated().piece(0));

            A.setIntegrationElements(edgeSpace);
            gsExprEvaluator<T> ev(A);

            // Set the discretization space
            auto u = A.getSpace(edgeSpace);

            // Create Mapper
            gsDofMapper map(edgeSpace);
            if (!m_optionList.getSwitch("interpolation"))
            {
                gsMatrix<index_t> act;
                for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
                map.finalize();

                u.setupMapper(map);

                gsMatrix<T> &fixedDofs = const_cast<expr::gsFeSpace<T> &>(u).fixedPart();
                fixedDofs.setZero(u.mapper().boundarySize(), 1);

                A.initSystem();
                A.assemble(u * u.tr()); // The Matrix is the same for each bf
                solver.compute(A.matrix());
                // [!The same setup for each bf!]
            }

            index_t n_plus = basis_plus.size();
            index_t n_minus = basis_minus.size();

            index_t bfID_init = 3;
            for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
            {
                gsTraceBasis<real_t> traceBasis(geo, beta, basis_plus, initSpace.basis(0), bdy, bfID, dir);

                if (m_optionList.getSwitch("interpolation"))
                {
                    //gsQuasiInterpolate<T>::Schoenberg(edgeSpace.basis(0), traceBasis, sol);
                    //result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                    gsMatrix<> anchors = edgeSpace.basis(0).anchors();
                    gsMatrix<> values = traceBasis.eval(anchors);
                    result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                }
                else
                {
                    A.initVector(); // Just the rhs

                    auto aa = A.getCoeff(traceBasis);

                    A.assemble(u * aa);

                    gsMatrix<T> solVector = solver.solve(A.rhs());

                    auto u_sol = A.getSolution(u, solVector);
                    gsMatrix<T> sol;
                    u_sol.extract(sol);

                    result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
                }
            }

            bfID_init = 2;
            for (index_t bfID = bfID_init; bfID < n_minus - bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
            {
                gsNormalDerivBasis<real_t> normalDerivBasis(geo, alpha, basis_minus, initSpace.basis(0), bdy, bfID,
                                                            dir);
                if (m_optionList.getSwitch("interpolation"))
                {
                    //gsQuasiInterpolate<T>::Schoenberg(edgeSpace.basis(0), traceBasis, sol);
                    //result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                    gsMatrix<> anchors = edgeSpace.basis(0).anchors();
                    gsMatrix<> values = normalDerivBasis.eval(anchors);
                    result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                }
                else
                {
                    A.initVector(); // Just the rhs

                    auto aa = A.getCoeff(normalDerivBasis);

                    A.assemble(u * aa);

                    gsMatrix<T> solVector = solver.solve(A.rhs());

                    auto u_sol = A.getSolution(u, solVector);
                    gsMatrix<T> sol;
                    u_sol.extract(sol);

                    result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
                }
            }

            // parametrizeBasisBack
            m_auxPatches[patchID].parametrizeBasisBack(result);

            basisEdgeResult.push_back(result);
        }
}


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::computeAuxTopology()
    {
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
        {
            if(m_auxPatches[i].getPatchRotated().orientation() == -1)
                m_auxPatches[i].swapAxis();
        }
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeInterfacePatches()
    {
        computeAuxTopology();

        gsMultiPatch<> temp_mp;
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
            temp_mp.addPatch(m_auxPatches[i].getPatchRotated());

        temp_mp.computeTopology();

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (temp_mp.interfaces()[0].second().side().index())
        {
            case 1:
                break;
            case 4: m_auxPatches[0].rotateParamClock();
                break;
            case 3: m_auxPatches[0].rotateParamAntiClock();
                break;
            case 2: m_auxPatches[0].rotateParamAntiClockTwice();
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (temp_mp.interfaces()[0].first().side().index())
        {
            case 3:
                break;
            case 4: m_auxPatches[1].rotateParamAntiClockTwice();
                break;
            case 2: m_auxPatches[1].rotateParamAntiClock();
                break;
            case 1: m_auxPatches[1].rotateParamClock();
                break;
            default:
                break;
        }
    } // reparametrizeInterfacePatches


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeSinglePatch(index_t side)
    {
        computeAuxTopology();

        if(m_auxPatches[0].getOrient())
        {
            switch (side)
            {
                case 3:
                    break;
                case 2:
                    m_auxPatches[0].rotateParamClock();
                    break;
                case 4:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    break;
                case 1:
                    m_auxPatches[0].rotateParamAntiClock();
                    break;
            }
        }
        else
        {
            switch (side)
            {
                case 1:
                    break;
                case 4:
                    m_auxPatches[0].rotateParamClock();
                    break;
                case 2:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    break;
                case 3:
                    m_auxPatches[0].rotateParamAntiClock();
                    break;
            }
        }
    } // reparametrizeSinglePatch

} // namespace gismo
