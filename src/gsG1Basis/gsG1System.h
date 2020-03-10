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
    void insertInterfaceEdge(gsMultiPatch<> & mp, index_t iID, index_t pID_0, index_t pdID_1);

protected:
    index_t dim_K, dim_E, dim_V;

    gsVector<> numBasisFunctions, numEdgeFunctions;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;


}; // class gsG1System

template<class T>
void gsG1System<T>::initialize(gsMultiPatch<> & mp, gsMultiBasis<> mb)
{
    // Number of the patches
    index_t numPatches = mp.nPatches();

    // Get the dimension of the basis functions for each patch
    numBasisFunctions.setZero(numPatches);
    numEdgeFunctions.setZero(4*(numPatches));
    numBasisFunctions[0] = mb.basis(0).size();
    {
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(0).component(0)); // 0 -> u, 1 -> v
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();
        // ATTENTION: With the first and the last interface basis functions
        numEdgeFunctions[0] = 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
        numEdgeFunctions[1] = numEdgeFunctions[0] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;

        basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(0).component(1)); // 0 -> u, 1 -> v
        m_p = basis_edge.maxDegree();
        m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        m_n = basis_edge.numElements();
        numEdgeFunctions[2] = numEdgeFunctions[1] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
        numEdgeFunctions[3] = numEdgeFunctions[2] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
    }
    for (size_t i = 1; i < mb.nBases(); i++ )
    {
        numBasisFunctions[i] = numBasisFunctions[i-1] + mb.basis(i).size();

        // Get the dimension for the spaces at the edges
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(i).component(0)); // 0 -> u, 1 -> v
        index_t m_p = basis_edge.maxDegree();
        index_t m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        index_t m_n = basis_edge.numElements();
        // ATTENTION: With the first and the last interface basis functions
        numEdgeFunctions[4*i + 0] = numEdgeFunctions[4*i - 1] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
        numEdgeFunctions[4*i + 1] = numEdgeFunctions[4*i + 0] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;

        basis_edge = dynamic_cast<gsBSplineBasis<> &>(mb.basis(i).component(1)); // 0 -> u, 1 -> v
        m_p = basis_edge.maxDegree();
        m_r = 1; // Here fixed to 1 TODO MORE GENERAL
        m_n = basis_edge.numElements();
        numEdgeFunctions[4*i + 2] = numEdgeFunctions[4*i + 1] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
        numEdgeFunctions[4*i + 3] = numEdgeFunctions[4*i + 2] + 2 * (m_p - m_r - 1) * (m_n - 1) + 2 * m_p + 1;
    }
    gsInfo << "Num Basis Functions " << numBasisFunctions << "\n";
    gsInfo << "Num Edge Functions " << numEdgeFunctions << "\n";

    dim_K = numBasisFunctions.last(); // interior basis
    dim_E = numEdgeFunctions.last();
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
void gsG1System<T>::insertInterfaceEdge(gsMultiPatch<> & mp, index_t iID ,index_t pID_0, index_t pdID_1)
{


    // Insert all coefficients of the g1 Basis at the interface
    for (size_t np = 0; np < 2; ++np) // two interface patches
    {
        for (size_t j = 0; j < mp.patch(np).coefs().size(); j++) // all the coefs
        {
            if (mp.patch(np).coefs().at(j) * mp.patch(np).coefs().at(j)  > 10e-25)
            {
                index_t jj;
                pID_0 == 0 ? jj = 0 : jj = numBasisFunctions[pID_0 -1 ];
                D_sparse.insert(iID,jj) = mp.patch(np).coefs().at(j); // TODO falsch
            }



        }
    } // edges




}

} // namespace
