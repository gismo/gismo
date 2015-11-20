/** @file gsGenericAssembler.h

    @brief Provides an assembler for common IGA matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsVisitorMass.h>
#include <gsAssembler/gsVisitorGradGrad.h>
#include <gsAssembler/gsVisitorMoments.h>

#include <gsPde/gsLaplacePde.h>


namespace gismo
{


template <typename T>
int estimateNonzerosPerRow(const gsBasis<T>& basis)
{
    int nnz = 1;
    for (int i = 0; i < basis.dim(); ++i) // to do: improve
        nnz *= 2 * basis.degree(i) + 1;
    return nnz;
}


template <typename T>
void localToGlobal(const gsMatrix<T>& localStiffness,
        const gsMatrix<unsigned>& localDofs,
        gsSparseMatrix<T>& K,
        bool symmetric)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        for (index_t j = 0; j < numActive; ++j)
        {
            const int jj = localDofs(j,0);
            // if matrix is symmetric, store only lower triangular part
            if (!symmetric || jj <= ii)
                K.coeffRef(ii, jj) += localStiffness(i, j);
        }
    }
}


/** Assemble the mass matrix in the parameter domain for a single basis,
 * i.e., with the geometry mapping being identity.
 */
template<class T>
void assembleParameterMass(const gsBasis<T>& basis, gsSparseMatrix<T>& M)
{
    const int n = basis.size();

    M.resize(n, n);
    M.reserve( gsVector<int>::Constant(n, estimateNonzerosPerRow(basis)) );

    gsMatrix<T> localMass;

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( basis.dim() );
    for (int i = 0; i < basis.dim(); ++i)
        numQuadNodes[i] = basis.degree(i) + 1;
    domIt->computeQuadratureRule( numQuadNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveFunctions().rows();
        domIt->evaluateBasis( 0 );

        localMass.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k];
            localMass.noalias() += weight * (domIt->basisValues().col(k) * domIt->basisValues().col(k).transpose());
        }

        //localToGlobal(localMass, domIt->activeFuncs, M, true);
        localToGlobal(localMass, domIt->activeFuncs, M, false);
    }

    M.makeCompressed();
}


/**
   @brief Assembles the mass, stiffness matrix on a given domain

   
   \ingroup Assembler
 */
template <class T>
class gsGenericAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    /// Constructor with gsMultiBasis
    gsGenericAssembler( const gsMultiPatch<T>    & patches,
                        const gsMultiBasis<T>    & bases,
                        const gsAssemblerOptions & opt = gsAssemblerOptions(),
                        const gsBoundaryConditions<T> * bc = NULL)
    : m_pde(patches)
    {
        if ( bc != NULL)
            m_pde.boundaryConditions() = *bc;

        Base::initialize(m_pde, bases, opt);
        refresh();
    }

    /* TODO (id geometry)
    // Constructor with gsMultiBasis
    gsGenericAssembler( gsMultiBasis<T> const         & bases,
                        bool conforming = false,
                        const gsBoundaryConditions<T> * bc = NULL)
    : Base(patches)
    {
        m_bases.push_back(bases);

        // Init mapper
        m_dofMappers.resize(1);
        if (bc)
            bases.getMapper(conforming, *bc, m_dofMappers.front() );
        else
            bases.getMapper(conforming, m_dofMappers.front() );
        m_dofs = m_dofMappers.front().freeSize();
        m_matrix.setZero();
    }
    */

    void refresh()
    {
        // Setup sparse system
        gsDofMapper mapper;
        m_bases[0].getMapper(m_options.dirStrategy,
                             m_options.intStrategy,
                             this->pde().bc(), mapper);
        m_system = gsSparseSystem<T>(mapper);
        const index_t nz = m_options.numColNz(m_bases[0][0]);
        m_system.reserve(nz, 1);
    }

    /// Mass assembly routine
    const gsSparseMatrix<T> & assembleMass()
    {
        // Clean the sparse system
        m_system.setZero();

        // Assemble mass integrals
        this->template push<gsVisitorMass<T> >();

        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.matrix();
    }

    /// Stiffness assembly routine
    const gsSparseMatrix<T> & assembleStiffness()
    {
        // Clean the sparse system
        m_system.setZero();

        // Assemble stiffness integrals
        this->template push<gsVisitorGradGrad<T> >();

        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.matrix();
    }

    /// Moments assembly routine
    const gsMatrix<T> & assembleMoments(const gsFunction<T> & func)
    {
        // Clean the right-hand side vector
        m_system.rhs().setZero();

        // Assemble moment integrals
        gsVisitorMoments<T> mom(func);
        this->push(mom);

        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.matrix();
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleMass(int patchIndex)
    {
        gsGenericAssembler<T> tmp(m_pde.patches().patch(patchIndex), 
                                  m_bases[patchIndex], m_options);
        tmp.assembleMass();
        m_system.matrix().swap(tmp.m_system.matrix());
        return m_system.matrix();
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleStiffness(int patchIndex)
    {
        gsGenericAssembler<T> tmp(m_pde.patches().patch(patchIndex), 
                                  m_bases[patchIndex],  m_options);
        tmp.assembleStiffness();
        m_system.matrix().swap(tmp.m_system.matrix());
        return m_system.matrix();
    }
    
    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    typename gsSparseMatrix<T>::fullView fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    const typename gsSparseMatrix<T>::constFullView fullMatrix() const
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }
    

private:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

private:

    gsLaplacePde<T> m_pde;
};



} // namespace gismo

