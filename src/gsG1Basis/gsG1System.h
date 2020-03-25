/** @file gsG1System.h

    @brief Create a G1-System for a Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once


namespace gismo
{
template<class T>
class gsG1System
{
public:

    gsG1System(gsMultiPatch<> & mp,
               gsMultiBasis<> & mb)
    {
        initialize(mp, mb);
    }

    void initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb);
    void finalize(gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsMatrix<> g1);
    gsMatrix<> solve(gsSparseMatrix<real_t> K, gsMatrix<> f);

    void insertInterfaceEdge(gsMultiPatch<> & mp, boundaryInterface item, index_t iID ,index_t bfID);
    void insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID ,size_t bfID);
    void insertVertex(gsMultiPatch<> & mp, std::vector<size_t> patchIndex, index_t vID, index_t nDofs, index_t bfID);

    void constructG1Solution(const gsMatrix<T> &solVector, std::vector<gsMultiPatch<>> & result, const gsMultiPatch<> & geo);

    size_t boundary_size() { return numBoundaryVertexFunctions.last() - numBoundaryEdgeFunctions[0]; }

    size_t sizePlusInterface(index_t i) { return  sizePlusInt[i]; };
    size_t sizePlusBoundary(index_t i) { return  sizePlusBdy[i]; };

    gsVector<> get_numBoundaryEdgeFunctions() {return numBoundaryEdgeFunctions; };
    gsVector<> get_numBoundaryVertexFunctions() {return numBoundaryVertexFunctions; };

    gsMatrix<> getSingleBasis(index_t global_row, index_t patchIdx) { return D_sparse.block(global_row, numBasisFunctions[patchIdx], 1, numBasisFunctions[patchIdx+1] - numBasisFunctions[patchIdx]); };

protected:
    size_t dim_K, dim_G1_Dofs, dim_G1_Bdy;

    gsVector<> numBasisFunctions, numInterfaceFunctions, numEdgeFunctions, numBoundaryEdgeFunctions, numVertexFunctions, numBoundaryVertexFunctions;
    gsVector<> kindOfVertex;
    gsVector<size_t> sizePlusInt, sizePlusBdy;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;
    gsMatrix<> m_g1;

}; // class gsG1System

template<class T>
void gsG1System<T>::initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb)
{
    // Number of the patches
    index_t numPatches = mp.nPatches();

    // Dimension of Basis plus at the edge
    sizePlusInt.setZero(mp.interfaces().size());
    sizePlusBdy.setZero(mp.boundaries().size());

    // Kind of vertex
    // -1 Boundary vertex
    // 0 Internal vertex
    // 1 Interface boundary vertey
    kindOfVertex.setZero(mp.vertices().size());

    // Get the dimension of the basis functions for each patch
    numBasisFunctions.setZero(numPatches+1);

    numInterfaceFunctions.setZero(mp.interfaces().size()+1);

    numEdgeFunctions.setZero(mp.boundaries().size()+1);
    numBoundaryEdgeFunctions.setZero(mp.boundaries().size()+1);

    numVertexFunctions.setZero(mp.vertices().size()+1);
    numBoundaryVertexFunctions.setZero(mp.vertices().size()+1);

    for (size_t i = 0; i < mb.nBases(); i++ )
        numBasisFunctions[i+1] = numBasisFunctions[i] + mb.basis(i).size();

    for (size_t i = 0; i < mp.interfaces().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.interfaces()[i].first().m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.interfaces()[i].first().patch).component(dir)); // If the interface matches!!!
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();

        numInterfaceFunctions[i+1] = numInterfaceFunctions[i] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p - 9;
        sizePlusInt[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }
    for (size_t i = 0; i < mp.boundaries().size(); i++)
    {
        // Get the dimension for the spaces at the edges
        index_t dir = mp.boundaries()[i].m_index < 3 ? 1 : 0;
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(mp.boundaries()[i].patch).component(dir)); // 0 -> u, 1 -> v
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();

        numBoundaryEdgeFunctions[i+1] = numBoundaryEdgeFunctions[i] + (m_p - m_r - 1) * (m_n - 1) + m_p + 1 - 6;
        numEdgeFunctions[i+1] = numEdgeFunctions[i] + (m_p - m_r - 1) * (m_n - 1) + m_p - 4;
        sizePlusBdy[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }

    for (size_t i = 0; i < mp.vertices().size(); i++)
    {
        if (mp.vertices()[i].size() == 1)
        {
            kindOfVertex[i] = -1; // Boundary vertex
            numVertexFunctions[i+1] = numVertexFunctions[i] + 1; // TODO
            numBoundaryVertexFunctions[i+1] = numBoundaryVertexFunctions[i] + 6; // TODO
        }
        else
        {
            gsMultiPatch<> temp_mp;
            for (size_t j = 0; j < mp.vertices()[i].size(); j++)
                temp_mp.addPatch(mp.patch(mp.vertices()[i][j].patch));

            temp_mp.computeTopology();
            if (mp.vertices()[i].size() == temp_mp.interfaces().size())
            {
                kindOfVertex[i] = 0; // Internal vertex
                numVertexFunctions[i+1] = numVertexFunctions[i] + 6;
                numBoundaryVertexFunctions[i+1] = numBoundaryVertexFunctions[i];
            }

            else
            {
                kindOfVertex[i] = 1; // Interface-Boundary vertex
                numVertexFunctions[i+1] = numVertexFunctions[i] + 3; // TODO
                numBoundaryVertexFunctions[i+1] = numBoundaryVertexFunctions[i] + 6; // TODO
            }
        }
    }
    numEdgeFunctions = numEdgeFunctions.array() + numInterfaceFunctions.last();
    numVertexFunctions = numVertexFunctions.array() + numEdgeFunctions.last();
    numBoundaryEdgeFunctions = numBoundaryEdgeFunctions.array() + numVertexFunctions.last();
    numBoundaryVertexFunctions = numBoundaryVertexFunctions.array() + numBoundaryEdgeFunctions.last();

    /*
    gsInfo << "Num Basis Functions " << numBasisFunctions << "\n";
    gsInfo << "Num Interface Functions " << numInterfaceFunctions << "\n";
    gsInfo << "Num Edges Functions " << numEdgeFunctions << "\n";
    gsInfo << "Num Boundary Edges Functions " << numBoundaryEdgeFunctions << "\n";
    gsInfo << "Num Vertex Functions " << numVertexFunctions << "\n";
    gsInfo << "Num Boundary Vertex Functions " << numBoundaryVertexFunctions << "\n";
    gsInfo << "Kind of Vertex Functions " << kindOfVertex << "\n";
    gsInfo << "Size of plus space Bdy  " << sizePlusBdy << "\n";
    gsInfo << "Size of plus space Int  " << sizePlusInt << "\n";
    */

    // Setting the final matrix
    dim_K = numBasisFunctions.last(); // interior basis dimension
    dim_G1_Dofs = numVertexFunctions.last() ; // edges basis dimension
    dim_G1_Bdy = numBoundaryVertexFunctions.last() - numBoundaryEdgeFunctions[0]; // vertex basis dimension

    // Full matrix
    D_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    // Boundary values
    m_g1.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, 1);
    m_g1.setZero();


}


template<class T>
void gsG1System<T>::constructG1Solution(const gsMatrix<T> & solVector, std::vector<gsMultiPatch<>> & result, const gsMultiPatch<> & geo)
{
    gsMultiPatch<T> init_mp;
    std::vector<gsMultiPatch<T>> g1Basis(geo.nPatches(), init_mp);

    for ( size_t rowInt = 0; rowInt < geo.interfaces().size(); rowInt++ ) // each interface edge
    {
        index_t patchIdx = geo.interfaces()[rowInt].first().patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numInterfaceFunctions[rowInt+1] - numInterfaceFunctions[rowInt]; i++)
        {
            index_t ii = numInterfaceFunctions[rowInt] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }

        patchIdx = geo.interfaces()[rowInt].second().patch;
        temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numInterfaceFunctions[rowInt+1] - numInterfaceFunctions[rowInt]; i++)
        {
            index_t ii = numInterfaceFunctions[rowInt] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }
    }


    for ( size_t rowEdge = 0; rowEdge < geo.boundaries().size(); rowEdge++ )
    {
        index_t patchIdx = geo.boundaries()[rowEdge].patch;
        gsTensorBSplineBasis<2,real_t> temp_basis = dynamic_cast<gsTensorBSplineBasis<2,real_t>  &>(geo.basis(patchIdx));
        for (size_t i = 0; i < numBoundaryEdgeFunctions[rowEdge+1] - numBoundaryEdgeFunctions[rowEdge]; i++) // each boundary edge
        {
            index_t ii = numBoundaryEdgeFunctions[rowEdge] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                m_g1.at(ii)));
        }
        for (size_t i = 0; i < numEdgeFunctions[rowEdge+1] - numEdgeFunctions[rowEdge]; i++) // each edge
        {
            index_t ii = numEdgeFunctions[rowEdge] + i;
            g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                solVector.at(ii)));
        }
    }


    for ( size_t rowVertex = 0; rowVertex < geo.vertices().size(); rowVertex++ )
    {
        std::vector<patchCorner> allcornerLists = geo.vertices()[rowVertex];

        for(size_t j = 0; j < allcornerLists.size(); j++)
        {
            index_t patchIdx = geo.vertices()[rowVertex][j].patch;
            gsTensorBSplineBasis<2, real_t>
                temp_basis = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(geo.basis(patchIdx));
            for (size_t i = 0; i < numBoundaryVertexFunctions[rowVertex+1] - numBoundaryVertexFunctions[rowVertex]; i++) // each boundary vertex
            {
                index_t ii = numBoundaryVertexFunctions[rowVertex] + i;
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                    m_g1.at(ii)));
            }
            for (size_t i = 0; i < numVertexFunctions[rowVertex+1] - numVertexFunctions[rowVertex]; i++) // each dofs vertex
            {
                index_t ii = numVertexFunctions[rowVertex] + i;
                g1Basis.at(patchIdx).addPatch(temp_basis.makeGeometry(D_sparse.block(ii,numBasisFunctions[patchIdx],1,temp_basis.size()).transpose() *
                    solVector.at(ii)));
            }
        }
    }

    result = g1Basis;
}


template<class T>
void gsG1System<T>::insertInterfaceEdge(gsMultiPatch<> & mp, boundaryInterface item, index_t iID ,index_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (size_t np = 0; np < 2; ++np) // two interface patches
        for (index_t j = 0; j < mp.patch(np).coefs().size(); j++) // all the coefs
            if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
            {
                index_t jj, ii;
                ii = numInterfaceFunctions[iID] + bfID;
                jj = numBasisFunctions[np == 0 ? item.first().patch : item.second().patch] + j;
                D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
            }

}

template<class T>
void gsG1System<T>::insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID ,size_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (index_t j = 0; j < mp.patch(0).coefs().size(); j++) // all the coefs
        if (mp.patch(0).coefs().at(j) * mp.patch(0).coefs().at(j)  > 10e-25)
        {
            index_t jj, ii;
            if (bfID < sizePlusBdy[bID] - 6)
                ii = numBoundaryEdgeFunctions[bID] + bfID;

            else
                ii = numEdgeFunctions[bID] + bfID - sizePlusBdy[bID] + 6;

            jj = numBasisFunctions[item.patch] + j;
            D_sparse.insert(ii,jj) = mp.patch(0).coefs().at(j);
        }
}

template<class T>
void gsG1System<T>::insertVertex(gsMultiPatch<> & mp, std::vector<size_t> patchIndex, index_t vID, index_t nDofs, index_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (size_t np = 0; np < mp.nPatches(); ++np) // for each patch which has the vertex
        for (index_t j = 0; j < mp.patch(np).coefs().size(); j++) // all the coefs
            if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
            {
                index_t jj, ii = -1;
                if (kindOfVertex[vID] == 0) // interior vertex
                    ii = numVertexFunctions[vID] + bfID; // all six belongs to Dofs
                else if (kindOfVertex[vID] == -1) // boundary
                {
                    if (bfID < nDofs)
                        ii = numVertexFunctions[vID] + bfID; // Dofs
                    else
                        ii = numBoundaryVertexFunctions[vID] + bfID - nDofs; // Boundary
                }
                else if (kindOfVertex[vID] == 1) // interface boundary
                {
                    if (bfID < nDofs)
                        ii = numVertexFunctions[vID] + bfID; // Dofs
                    else
                        ii = numBoundaryVertexFunctions[vID] + bfID - nDofs; // Boundary
                }
                jj = numBasisFunctions[patchIndex[np]] + j;
                D_sparse.insert(ii,jj) = mp.patch(np).coefs().at(j);
            }
}

template<class T>
void gsG1System<T>::finalize(gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsMatrix<> g1)
{

    gsSparseMatrix<> B_0_sparse, B_boundary_sparse, temp;
    B_0_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_G1_Dofs + dim_G1_Bdy + dim_K);
    B_0_sparse.reserve(dim_G1_Dofs + dim_K);
    B_0_sparse.setZero();

    B_boundary_sparse.resize(dim_G1_Dofs + dim_G1_Bdy + dim_K, dim_G1_Dofs + dim_G1_Bdy + dim_K);
    B_boundary_sparse.reserve(dim_G1_Bdy);
    B_boundary_sparse.setZero();

    // Add for the G1 Dofs to 1
    for(size_t i = 0; i < dim_G1_Dofs; i++)
        B_0_sparse.insert(i,i) = 1;

    // Add for the Boundary the edges to 1
    for(size_t i = 0; i < dim_G1_Bdy; i++)
    {
        index_t ii = dim_G1_Dofs + i;
        B_boundary_sparse.insert(ii,ii) = 1;
    }

    // Add the identity to the end of D
    // Construct the internal matrix
    for(size_t np = 0; np < mb.nBases(); np++)
    {
        index_t dim_u = mb.basis(np).component(0).size();
        index_t dim_v = mb.basis(np).component(1).size();
        for(index_t j = 2; j < dim_v - 2; j++)
            for(index_t i = 2; i < dim_u - 2; i++)
            {
                index_t ii = dim_G1_Dofs + dim_G1_Bdy + numBasisFunctions[np] + j*dim_u + i;
                index_t jj = numBasisFunctions[np] + j*dim_u + i;
                B_0_sparse.insert(ii,ii) = 1;
                D_sparse.insert(ii,jj) = 1;
            }
    }
    D_sparse.makeCompressed();

    // Getting D_0
    D_0_sparse = B_0_sparse * D_sparse;
    D_0_sparse.makeCompressed();

    // Compute D_Boundary
    D_boundary_sparse = B_boundary_sparse * D_sparse;
    D_boundary_sparse.makeCompressed();

    // Set up the boundary vector
    m_g1.block(dim_G1_Dofs,0,dim_G1_Bdy,1) = g1;

}

template<class T>
gsMatrix<> gsG1System<T>::solve(gsSparseMatrix<real_t> K, gsMatrix<> f)
{
    gsSparseMatrix<real_t> A = D_0_sparse * K * D_0_sparse.transpose();
    gsVector<real_t> F = D_0_sparse * f - D_0_sparse * K * D_boundary_sparse.transpose() * m_g1;

    gsSparseSolver<real_t>::CGDiagonal solver;
    solver.compute(A);
    gsMatrix<> solVector = solver.solve(F);

    return solVector;
}

} // namespace
