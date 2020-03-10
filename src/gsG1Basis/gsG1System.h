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

    void insertInterfaceEdge(gsMultiPatch<> & mp, boundaryInterface item, index_t iID ,index_t bfID);
    void insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID ,index_t bfID);

    index_t boundarySize() {return numBoundaryEdgeFunctions.last(); };
    index_t sizePlusInterface(index_t i) { return  sizePlusInt[i]; };
    size_t sizePlusBoundary(index_t i) { return  sizePlusBdy[i]; };
    //size_t numInterfaceFcts(index_t i) { return  i == 0 ? numInterfaceFunctions[1] : numInterfaceFunctions[i] - numInterfaceFunctions[i-1]; };

protected:
    index_t dim_K, dim_E, dim_V;

    gsVector<> numBasisFunctions, numInterfaceFunctions, numBoundaryEdgeFunctions;
    gsVector<size_t> sizePlusInt, sizePlusBdy;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;


}; // class gsG1System

template<class T>
void gsG1System<T>::initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb)
{
    // Number of the patches
    index_t numPatches = mp.nPatches();

    sizePlusInt.setZero(mp.interfaces().size());
    sizePlusBdy.setZero(mp.boundaries().size());

    // Get the dimension of the basis functions for each patch
    numBasisFunctions.setZero(numPatches+1);
    numInterfaceFunctions.setZero(mp.interfaces().size()+1);
    numBoundaryEdgeFunctions.setZero(mp.boundaries().size()+1);
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
        // ATTENTION: With the first and the last interface basis functions
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
        // ATTENTION: With the first and the last interface basis functions
        numBoundaryEdgeFunctions[i+1] = numBoundaryEdgeFunctions[i] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p - 9;
        sizePlusBdy[i] = (m_p - m_r - 1) * (m_n - 1) + m_p + 1;
    }

    gsInfo << "Num Basis Functions " << numBasisFunctions << "\n";
    gsInfo << "Num Interface Functions " << numInterfaceFunctions << "\n";
    gsInfo << "Num Boundary Functions " << numBoundaryEdgeFunctions << "\n";

    dim_K = numBasisFunctions.last(); // interior basis
    dim_E = numInterfaceFunctions.last() + numBoundaryEdgeFunctions.last();
    dim_V = 6 * 4 * numPatches;

    // Full matrix
    D_sparse.resize(dim_E + dim_V + dim_K, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_E + dim_V + dim_K, dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_E + dim_V + dim_K , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

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
void gsG1System<T>::insertBoundaryEdge(gsMultiPatch<> & mp, patchSide item, index_t bID ,index_t bfID)
{
    // Insert all coefficients of the g1 Basis at the interface
    for (index_t j = 0; j < mp.patch(0).coefs().size(); j++) // all the coefs
        if (mp.patch(0).coefs().at(j) * mp.patch(0).coefs().at(j)  > 10e-25)
        {
            index_t jj, ii;
            ii = numInterfaceFunctions.last() + numBoundaryEdgeFunctions[bID] + bfID;
            jj = numBasisFunctions[item.patch] + j;
            D_sparse.insert(ii,jj) = mp.patch(0).coefs().at(j);
        }
}

} // namespace
