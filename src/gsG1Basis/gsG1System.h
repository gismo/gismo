/** @file gsG1System.h

    @brief Create a G1-System for a Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsMatrix/gsMatrix.h>

namespace gismo
{
template<class T>
class gsG1System
{
public:

    gsG1System(gsMultiPatch<T> & g1B_L,
               gsMultiPatch<T> & g1B_R,
               gsMatrix<T> dofs,
               gsDofMapper map_g,
               gsDofMapper map_b,
               gsDofMapper map_i,
               index_t dim,
               index_t n_t,
               index_t n_b)
        : g1Basis_L(g1B_L), g1Basis_R(g1B_R), m_dofs(dofs), map_global(map_g), map_boundary(map_b),
        map_interface(map_i), dim_K(dim), n_tilde(n_t), n_bar(n_b)
    {
        initialize();
    }

    inline void initialize();

    void assemble();

    void constructSolution_G1(const gsMatrix<T>& solVector,
                              gsMultiPatch<T>& g1B_L,
                              gsMultiPatch<T>& g1B_R);

    gsSparseMatrix<T> get_D_sparse() { return D_sparse; };
    gsSparseMatrix<T> get_D_0_sparse() { return D_0_sparse; };
    gsSparseMatrix<T> get_D_boundary_sparse() { return D_boundary_sparse; };

    gsMatrix<T> get_g() { return m_g; };

protected:
    gsMultiPatch<T> g1Basis_L, g1Basis_R;

    gsMatrix<T> m_dofs, m_g;

    gsDofMapper map_global, map_boundary, map_interface;

    index_t dim_K, n_tilde, n_bar;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;

    gsMatrix<T> g1Dofs;

}; // class gsG1System

template<class T>
void gsG1System<T>::constructSolution_G1(const gsMatrix<T>& solVector,
                                         gsMultiPatch<T> & g1B_L,
                                         gsMultiPatch<T> & g1B_R)
{

    for (unsigned i = 0; i < g1B_L.nPatches(); i++)
    {
        if (i == 0)
        {
            g1B_L.patch(i).setCoefs(g1B_L.patch(i).coefs() * g1Dofs.at(0));
            g1B_R.patch(i).setCoefs(g1B_R.patch(i).coefs() * g1Dofs.at(0));
        }
        else if (i == n_tilde - 1)
        {
            g1B_L.patch(i).setCoefs(g1B_L.patch(i).coefs() * g1Dofs.at(1));
            g1B_R.patch(i).setCoefs(g1B_R.patch(i).coefs() * g1Dofs.at(1));
        }
        else if (i == n_tilde)
        {
            g1B_L.patch(i).setCoefs(g1B_L.patch(i).coefs() * g1Dofs.at(2));
            g1B_R.patch(i).setCoefs(g1B_R.patch(i).coefs() * g1Dofs.at(2));
        }
        else if (i == n_tilde + n_bar - 1)
        {
            g1B_L.patch(i).setCoefs(g1B_L.patch(i).coefs() * g1Dofs.at(3));
            g1B_R.patch(i).setCoefs(g1B_R.patch(i).coefs() * g1Dofs.at(3));
        }
/*        else if (i == n_tilde - 2 && north)
        {
            basisG1_L.patch(i).setCoefs(basisG1_L.patch(i).coefs() * g1Dofs.at(5));
            basisG1_R.patch(i).setCoefs(basisG1_R.patch(i).coefs() * g1Dofs.at(5));
        }
        else if (i == 1 && south)
        {
            basisG1_L.patch(i).setCoefs(basisG1_L.patch(i).coefs() * g1Dofs.at(4));
            basisG1_R.patch(i).setCoefs(basisG1_R.patch(i).coefs() * g1Dofs.at(4));
        }*/
        else
        {
            g1B_L.patch(i).setCoefs(g1B_L.patch(i).coefs() * solVector.at(i));
            g1B_R.patch(i).setCoefs(g1B_R.patch(i).coefs() * solVector.at(i));
        }

    }
}

template<class T>
void gsG1System<T>::initialize()
{
    // Full matrix
    D_sparse.resize((g1Basis_L.nPatches() + dim_K) , dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize((g1Basis_L.nPatches() + dim_K) , dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize((g1Basis_L.nPatches() + dim_K) , dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    m_g.resize(g1Basis_L.nPatches() + dim_K,1);
    m_g.setZero();

} // initialize

template<class T>
void gsG1System<T>::assemble()
{
    // First insert all coefficients of the g1 Basis
    for (size_t i = 0; i < g1Basis_L.nPatches() ; i++)
    {
        for (index_t j = 0; j < g1Basis_L.patch(i).coefs().size(); j++)
        {
            if (g1Basis_L.patch(i).coefs().at(j) * g1Basis_L.patch(i).coefs().at(j) > 10e-30)
                D_sparse.insert(i,j) = g1Basis_L.patch(i).coefs().at(j);

            if (g1Basis_R.patch(i).coefs().at(j) * g1Basis_R.patch(i).coefs().at(j) > 10e-30)
                D_sparse.insert(i,j+dim_K/2) = g1Basis_R.patch(i).coefs().at(j);
        }
    }
    //D_sparse.prune(0,10e-15);


    // Set identity in the lower part of D
    for (index_t i = 0; i < dim_K; i++)
        D_sparse.insert(i+g1Basis_L.nPatches(),i) = 1;
    D_sparse.makeCompressed();

    // Create identity matrix for D_0
    gsSparseMatrix<real_t> identity_sparse((g1Basis_L.nPatches() + dim_K),(g1Basis_L.nPatches() + dim_K));
    identity_sparse.reserve((g1Basis_L.nPatches() + dim_K));
    identity_sparse.setIdentity();

    // Setting the normal basis function with NO influence at the interface
    // AND boundary functions
    for (index_t p = 0; p < 2; p++) // Assume that we have only two Patches
    {
        for (index_t i = 0; i < dim_K/2; i++) // /2 bc of two patches
        {
            if (map_interface.is_boundary(i,p) || map_boundary.is_boundary(i,p))
                identity_sparse.coeffRef(map_global.index(i,p)+g1Basis_L.nPatches(),
                    map_global.index(i,p)+g1Basis_L.nPatches()) = 0;
        }
    }

    // Boundary for interface
    identity_sparse.coeffRef(0,0) = 0;
    //identity_sparse.coeffRef(1,1) = 0;
    //identity_sparse.coeffRef(n_tilde-2,n_tilde-2) = 0;
    identity_sparse.coeffRef(n_tilde-1,n_tilde-1) = 0;
    identity_sparse.coeffRef(n_tilde,n_tilde) = 0;
    identity_sparse.coeffRef(n_tilde + n_bar -1,n_tilde + n_bar -1) = 0;

    // Sparse matrix without boundary
    D_0_sparse = identity_sparse * D_sparse;
    D_0_sparse.makeCompressed();

    D_boundary_sparse = D_sparse - D_0_sparse;
    D_boundary_sparse.makeCompressed();

    // Dofs
    g1Dofs = m_dofs.bottomRows(6);
    // for non-interface boundary
    for (index_t i = 0; i < dim_K/2; i++)
    {
        if (map_boundary.is_boundary(i,0))
            m_g.at( g1Basis_L.nPatches() + i) = m_dofs.at(map_boundary.bindex(i,0)); // patch 0

        if (map_boundary.is_boundary(i,1))
            m_g.at( g1Basis_L.nPatches() + dim_K/2 + i) = m_dofs.at(map_boundary.bindex(i,1)); // patch 1
    }

    // for interface boundary
    m_g.at(0) = g1Dofs.at(0);
    m_g.at(n_tilde - 1) = g1Dofs.at(1);
    m_g.at(n_tilde) = g1Dofs.at(2);
    m_g.at(n_tilde + n_bar - 1) = g1Dofs.at(3);

}

} // namespace gismo
