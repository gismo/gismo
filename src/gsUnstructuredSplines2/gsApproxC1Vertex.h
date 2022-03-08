/** @file gsApproxC1Vertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1Utils.h>

#include <gsUnstructuredSplines/gsContainerBasis.h>
#include <gsUnstructuredSplines/gsPatchReparameterized.h>
#include <gsUnstructuredSplines/gsApproxGluingData.h>


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
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1]));
        }

        reparametrizeVertexPatches();

        // Compute Sigma
        real_t sigma = computeSigma(m_vertexIndices);

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            C1AuxPatchContainer auxPatchSingle;
            auxPatchSingle.push_back(m_auxPatches[i]);

            std::vector<patchSide> containingSides;
            patchCorner pC(m_patchesAroundVertex[i], m_vertexIndices[i]);
            pC.getContainingSides(d, containingSides);

            if (containingSides.at(0).side() < 3) // If isInterface_1 == v, then switch
            {
                patchSide side_temp = containingSides[0];
                containingSides[0] = containingSides[1];
                containingSides[1] = side_temp;
            }

            std::vector<bool> isInterface(2);
            isInterface[0] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(0).side()));
            isInterface[1] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(1).side()));

            // Compute Gluing data
            gsApproxGluingData<d, T> approxGluingData(auxPatchSingle, m_optionList, containingSides, isInterface);

            //Problem setup
            std::vector<gsBSpline<T>> alpha, beta;
            std::vector<gsBSplineBasis<T>> basis_plus, basis_minus;

            std::vector<bool> kindOfEdge;

            alpha.resize(2); beta.resize(2); basis_plus.resize(2); basis_minus.resize(2);

            kindOfEdge.resize(2);

            gsGeometry<T> & geo = auxPatchSingle[0].getPatchRotated();
            gsMultiBasis<T> initSpace(auxPatchSingle[0].getBasisRotated().piece(0));
            for (size_t dir = 0; dir < containingSides.size(); ++dir)
            {
                index_t localdir = auxPatchSingle[0].getMapIndex(containingSides[dir].index()) < 3 ? 1 : 0;

                if (isInterface[localdir])
                {
                    alpha[dir] = approxGluingData.alphaS(dir);
                    beta[dir] = approxGluingData.betaS(dir);
                }

                kindOfEdge[localdir] = isInterface[dir];

                gsBSplineBasis<T> b_plus, b_minus;
                createPlusSpace(geo, initSpace.basis(0), dir, b_plus);
                createMinusSpace(geo, initSpace.basis(0), dir, b_minus);

                basis_plus[dir] = b_plus;
                basis_minus[dir] = b_minus;
            }

            gsSparseSolver<real_t>::SimplicialLDLT solver;
            gsExprAssembler<> A(1, 1);

            // Elements used for numerical integration
            gsMultiBasis<T> vertexSpace(auxPatchSingle[0].getBasisRotated().piece(m_vertexIndices[i] + 4));
            A.setIntegrationElements(vertexSpace);
            gsExprEvaluator<T> ev(A);

            // Set the discretization space
            auto u = A.getSpace(vertexSpace);

            // Create Mapper
            gsDofMapper map(vertexSpace);
            if (!m_optionList.getSwitch("interpolation"))
            {
                gsMatrix<index_t> act;
                for (index_t dir = 0; dir < vertexSpace.basis(0).domainDim(); dir++)
                    for (index_t i = 3 * vertexSpace.basis(0).degree(dir) + 1; i < vertexSpace.basis(0).component(
                            1 - dir).size(); i++) // only the first two u/v-columns are Dofs (0/1)
                    {
                        act = vertexSpace.basis(0).boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                        //map.markBoundary(0, act); // Patch 0
                    }
                map.finalize();
                u.setupMapper(map);

                gsMatrix<T> &fixedDofs = const_cast<expr::gsFeSpace<T> &>(u).fixedPart();
                fixedDofs.setZero(u.mapper().boundarySize(), 1);

                A.initSystem();
                A.assemble(u * u.tr());
                solver.compute(A.matrix());
            }

            // Create Basis functions
            gsMultiPatch<T> result_1;
            for (index_t bfID = 0; bfID < 6; bfID++)
            {
                gsVertexBasis<T> vertexBasis(geo, initSpace.basis(0), alpha, beta, basis_plus, basis_minus, sigma,
                                             kindOfEdge, bfID);
                if (m_optionList.getSwitch("interpolation"))
                {
                    //gsQuasiInterpolate<T>::Schoenberg(edgeSpace.basis(0), traceBasis, sol);
                    //result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                    gsMatrix<> anchors = vertexSpace.basis(0).anchors();
                    gsMatrix<> values = vertexBasis.eval(anchors);
                    result_1.addPatch(vertexSpace.basis(0).interpolateAtAnchors(give(values)));
                }
                else
                {
                    A.initVector();

                    auto aa = A.getCoeff(vertexBasis);
                    A.assemble(u * aa);

                    gsMatrix<T> solVector = solver.solve(A.rhs());

                    auto u_sol = A.getSolution(u, solVector);
                    gsMatrix<T> sol;
                    u_sol.extract(sol);

                    result_1.addPatch(vertexSpace.basis(0).makeGeometry(give(sol)));
                }
            }
            //Problem setup end

            // Store temporary
            basisVertexResult.push_back(result_1);
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
            //if (m_patchesAroundVertex.size() == 2)
            //    gsWriteParaview(basisVertexResult[0], "vertex_basis", 20000);
        }
    }

    void reparametrizeVertexPatches();

    void checkOrientation(size_t i);

    real_t computeSigma(const std::vector<size_t> & vertexIndices);

    void computeKernel();

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

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Vertex.hpp)
#endif
