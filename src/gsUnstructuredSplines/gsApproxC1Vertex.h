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

            index_t vertex_1 = m_vertexIndices[i];
            std::vector<index_t> sideContainer;
            if (auxPatchSingle[0].getOrient() == 0) // not rotated
            {
                switch (vertex_1) // corner
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
                switch (vertex_1) // corner
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

            //Problem setup
            std::vector<gsBSplineBasis<T>> basis_plus;
            std::vector<gsBSplineBasis<T>> basis_minus;
            std::vector<gsBSplineBasis<T>> basis_geo;

            std::vector<gsBSpline<T>> alpha;
            std::vector<gsBSpline<T>> beta;

            std::vector<bool> kindOfEdge;

            basis_plus.resize(2);
            basis_minus.resize(2);
            basis_geo.resize(2);

            alpha.resize(2);
            beta.resize(2);

            kindOfEdge.resize(2);

            for (size_t dir = 0; dir < sideContainer.size(); ++dir)
            {
                index_t localdir = auxPatchSingle[0].getMapIndex(sideContainer[dir]) < 3 ? 1 : 0;

                basis_plus[localdir] = dynamic_cast<gsBSplineBasis<T>&>(auxPatchSingle[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 0));
                basis_minus[localdir] = dynamic_cast<gsBSplineBasis<T>&>(auxPatchSingle[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 1));
                basis_geo[localdir] = dynamic_cast<gsBSplineBasis<T>&>(auxPatchSingle[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 2));

                if (isInterface[localdir])
                {
                    alpha[dir] = approxGluingData.alphaS(dir);
                    beta[dir] = approxGluingData.betaS(dir);
                }


                kindOfEdge[localdir] = isInterface[dir];
            }
            gsGeometry<T> &geo = auxPatchSingle[0].getPatchRotated();

            // Create Basis functions
            gsMultiPatch<> result_1;
            for (index_t bfID = 0; bfID < 6; bfID++) {
                gsSparseSolver<real_t>::SimplicialLDLT solver;
                gsExprAssembler<> A(1, 1);

                typedef gsExprAssembler<>::variable variable;
                typedef gsExprAssembler<>::space space;
                typedef gsExprAssembler<>::solution solution;

                // Elements used for numerical integration
                gsMultiBasis<T> vertexSpace(auxPatchSingle[0].getBasisRotated().getBasis(m_vertexIndices[i] + 4));
                A.setIntegrationElements(vertexSpace);
                gsExprEvaluator<> ev(A);

                // Set the discretization space
                space u = A.getSpace(vertexSpace);

                // Create Mapper
                gsDofMapper map(vertexSpace);
                gsMatrix<index_t> act;
                for (index_t dir = 0; dir < vertexSpace.basis(0).domainDim(); dir++)
                    for (index_t i = 2*vertexSpace.basis(0).degree(dir)+1; i < vertexSpace[0].component(dir).size(); i++) // only the first two u/v-columns are Dofs (0/1)
                    {
                        act = vertexSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                        map.markBoundary(0, act); // Patch 0
                    }
                map.finalize();

                gsBoundaryConditions<> bc_empty;
                bc_empty.addCondition(1, condition_type::dirichlet, 0); // Doesn't matter which side
                u.setupMapper(map);

                A.initSystem();

                gsVertexBasis<real_t> vertexBasis(geo, basis_plus, basis_minus, basis_geo, alpha, beta, sigma,
                                                  kindOfEdge, bfID);
                auto aa = A.getCoeff(vertexBasis);

                A.assemble(u * u.tr(), u * aa);

                solver.compute(A.matrix());
                gsMatrix<> solVector = solver.solve(A.rhs());

                solution u_sol = A.getSolution(u, solVector);
                gsMatrix<> sol;
                u_sol.extract(sol);

                //gsDebugVar(sol - result_1.patch(bfID).coefs());
                result_1.addPatch(vertexSpace.basis(0).makeGeometry(give(sol)));
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
