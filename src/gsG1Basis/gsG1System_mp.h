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
class gsG1System_mp
{
public:

    gsG1System_mp(std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>  & g1B,
               gsMatrix<T> dofs,
               gsDofMapper map_g,
               gsDofMapper map_b,
               gsDofMapper map_i,
               index_t dim,
               std::vector<index_t> n_t,
               std::vector<index_t> n_b,
               gsMultiBasis<T> basis)
        : g1Basis(g1B), m_dofs(dofs), map_global(map_g), map_boundary(map_b),
        map_interface(map_i), dim_K(dim), n_tilde(n_t), n_bar(n_b), m_basis(basis)
    {
        // Compute dimension of n_D
        dim_D = dim_K;

        for (std::vector<index_t>::iterator it = n_tilde.begin(); it!=n_tilde.end(); ++it)
            dim_D += *it;
        for (std::vector<index_t>::iterator it = n_bar.begin(); it!=n_bar.end(); ++it) // is the same dim as n_tilde
            dim_D += *it;

        initialize();
    }

    inline void initialize();

    void assemble();

    void constructSolution_G1(const gsMatrix<T>& solVector,
                              std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> & g1B);

    gsSparseMatrix<T> get_D_sparse() { return D_sparse; };
    gsSparseMatrix<T> get_D_0_sparse() { return D_0_sparse; };
    gsSparseMatrix<T> get_D_boundary_sparse() { return D_boundary_sparse; };

    gsMatrix<T> get_g() { return m_g; };

    index_t get_dim_g1() { return dim_D - dim_K; };

protected:
    std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>  g1Basis;

    gsMatrix<T> m_dofs, m_g;

    gsDofMapper map_global, map_boundary, map_interface;

    index_t dim_K, dim_D;

    std::vector<index_t> n_tilde, n_bar;

    gsMultiBasis<T> m_basis;

    gsSparseMatrix<T> D_sparse, D_0_sparse, D_boundary_sparse;

    gsMatrix<T> g1Dofs;

}; // class gsG1System

template<class T>
void gsG1System_mp<T>::constructSolution_G1(const gsMatrix<T>& solVector,
                                            std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> & g1B)
{
    typedef std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> typedef_g1;
    typedef std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>> typedef_iFace;
    for (typedef_g1::iterator it_patchID = g1B.begin(); it_patchID!=g1B.end(); ++it_patchID)
    {
        for (typedef_iFace::iterator it_sideID = it_patchID->second.begin(); it_sideID!=it_patchID->second.end(); ++it_sideID)
        {
            for (std::map<index_t, gsMultiPatch<real_t>>::iterator it_g1Basis = it_sideID->second.begin(); it_g1Basis!=it_sideID->second.end(); ++it_g1Basis)
            {
                gsMultiPatch<T> & mp = it_g1Basis->second;
                index_t i = it_g1Basis->first; // Interface index
                index_t dim_iFace = 0;
                for (index_t jj = 0; jj < i; ++jj)
                    dim_iFace += n_tilde.at(jj) + n_bar.at(jj);

                for (index_t j = 0; j < n_tilde.at(i) + n_bar.at(i); ++j) // basis function index
                {
                    if (j == 0)
                    {
                        mp.patch(j).setCoefs(mp.patch(j).coefs() * g1Dofs.at(i*4));
                    }
                    else if (j == n_tilde.at(i) - 1)
                    {
                        mp.patch(j).setCoefs(mp.patch(j).coefs() * g1Dofs.at(i*4 + 1));

                    }
                    else if (j == n_tilde.at(i))
                    {
                        mp.patch(j).setCoefs(mp.patch(j).coefs() * g1Dofs.at(i*4 + 2));

                    }
                    else if (j == n_tilde.at(i) + n_bar.at(i) - 1)
                    {
                        mp.patch(j).setCoefs(mp.patch(j).coefs() * g1Dofs.at(i*4 + 3));

                    }
                    else
                    {
                        mp.patch(j).setCoefs(mp.patch(j).coefs() * solVector.at(dim_iFace + j));
                    }
                }
            }
        }
    }
}

template<class T>
void gsG1System_mp<T>::initialize()
{
    // Full matrix
    D_sparse.resize(dim_D, dim_K);
    D_sparse.reserve(3*dim_K);
    D_sparse.setZero();

    // Without boundary
    D_0_sparse.resize(dim_D , dim_K);
    D_0_sparse.reserve(3*dim_K);
    D_0_sparse.setZero();

    // Only boundary
    D_boundary_sparse.resize(dim_D, dim_K);
    D_boundary_sparse.reserve(3*dim_K);
    D_boundary_sparse.setZero();

    m_g.resize(dim_D,1);
    m_g.setZero();

} // initialize

template<class T>
void gsG1System_mp<T>::assemble()
{
    typedef std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> typedef_g1;
    typedef std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>> typedef_iFace;

    // First insert all coefficients of the g1 Basis
    for (typedef_g1::iterator it_patchID = g1Basis.begin(); it_patchID!= g1Basis.end(); ++it_patchID)
    {
        index_t patchId = it_patchID->first;

        for (typedef_iFace::iterator it_sideID = it_patchID->second.begin(); it_sideID!=it_patchID->second.end(); ++it_sideID)
        {
            for (std::map<index_t, gsMultiPatch<real_t>>::iterator it_g1Basis = it_sideID->second.begin(); it_g1Basis!=it_sideID->second.end(); ++it_g1Basis)
            {
                gsMultiPatch<T> mp_g1 = it_g1Basis->second;
                for (size_t i = 0; i < mp_g1.nPatches(); ++i)
                {
                    for (index_t j = 0; j < mp_g1.patch(i).coefs().size(); j++)
                    {
                        gsMatrix<unsigned > localDof(1,1), globalDof;
                        localDof << j;
                        map_global.localToGlobal(localDof,patchId,globalDof);

                        index_t dim_iFace = i;
                        for (index_t ii = 0; ii < it_g1Basis->first; ++ii)
                            dim_iFace += n_tilde.at(ii) + n_bar.at(ii);

                        if (mp_g1.patch(i).coefs().at(j) * mp_g1.patch(i).coefs().at(j) > 10e-30)
                            D_sparse.insert(dim_iFace,globalDof.at(0)) = mp_g1.patch(i).coefs().at(j);

                    }
                }
            }
        }
    }
    //D_sparse.prune(0,10e-15);

    // Set identity in the lower part of D
    for (index_t i = 0; i < dim_K; i++)
        D_sparse.insert(i+dim_D - dim_K,i) = 1;
    D_sparse.makeCompressed();

    // Create identity matrix for D_0
    gsSparseMatrix<real_t> identity_sparse(dim_D,dim_D);
    identity_sparse.reserve(dim_D);
    identity_sparse.setIdentity();

    // Setting the normal basis function with NO influence at the interface
    // AND boundary functions
    for (unsigned p = 0; p < m_basis.nBases(); p++) // Assume that we have only two Patches
    {
        for (index_t i = 0; i < m_basis.basis(p).size(); i++) // /2 bc of two patches
        {
            if (map_interface.is_boundary(i,p) || map_boundary.is_boundary(i,p))
                identity_sparse.coeffRef(map_global.index(i,p)+dim_D-dim_K,
                    map_global.index(i,p)+dim_D-dim_K) = 0;
        }
    }

    // Boundary for interface
    index_t dim_iFace = 0;
    for (unsigned i = 0; i < n_tilde.size(); ++i) // each interface
    {
        // also at the vertices the basis functions is deleted
        // maybe rework: we need a vertices mapper

        identity_sparse.coeffRef(dim_iFace,dim_iFace) = 0;
        //identity_sparse.coeffRef(1,1) = 0;
        //identity_sparse.coeffRef(n_tilde-2,n_tilde-2) = 0;
        identity_sparse.coeffRef(dim_iFace+ n_tilde.at(i)-1,dim_iFace+ n_tilde.at(i)-1) = 0;
        identity_sparse.coeffRef(dim_iFace+ n_tilde.at(i), dim_iFace+ n_tilde.at(i)) = 0;
        identity_sparse.coeffRef(dim_iFace+ n_tilde.at(i) + n_bar.at(i) -1, dim_iFace+ n_tilde.at(i) + n_bar.at(i) -1) = 0;

        dim_iFace += n_tilde.at(i) + n_bar.at(i);
    }

    // Sparse matrix without boundary
    D_0_sparse = identity_sparse * D_sparse;
    D_0_sparse.makeCompressed();

    D_boundary_sparse = D_sparse - D_0_sparse;
    D_boundary_sparse.makeCompressed();

    // Dofs
    g1Dofs = m_dofs.bottomRows(4 * n_tilde.size());
    // for non-interface boundary
    for (index_t p = 0; p < m_basis.nBases(); p++) // Assume that we have only two Patches
    {
        for (index_t i = 0; i < m_basis.basis(p).size(); i++) // /2 bc of two patches
        {
            if (map_boundary.is_boundary(i, p))
                m_g.at(dim_D - dim_K + map_global.index(i,p)) = m_dofs.at(map_boundary.bindex(i, p)); // patch 0
        }
    }

    // Boundary for interface
    dim_iFace = 0;
    for (unsigned i = 0; i < n_tilde.size(); ++i) // each interface
    {
        // also at the vertices the basis functions is deleted
        // maybe rework: we need a vertices mapper

        m_g.at(dim_iFace) = g1Dofs.at(i*4);
        //identity_sparse.coeffRef(1,1) = 0;
        //identity_sparse.coeffRef(n_tilde-2,n_tilde-2) = 0;
        m_g.at(dim_iFace + n_tilde.at(i) - 1) = g1Dofs.at(i*4 + 1);
        m_g.at(dim_iFace + n_tilde.at(i)) = g1Dofs.at(i*4 + 2);
        m_g.at(dim_iFace + n_tilde.at(i) + n_bar.at(i) - 1) = g1Dofs.at(i*4 + 3);

        dim_iFace += n_tilde.at(i) + n_bar.at(i);
    }


}

} // namespace gismo
