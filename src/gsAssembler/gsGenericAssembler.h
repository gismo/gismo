/** @file gsGenericAssembler.h

    @brief Provides an assembler for common IGA matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsVisitorMass.h>
#include <gsAssembler/gsVisitorGradGrad.h>

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

        localToGlobal(localMass, domIt->activeFuncs, M, true);
    }

    M.makeCompressed();
}


/**
   @brief Assembles the mass, stiffness matrix on a given domain

   
   \ingroup Assembler
 */
template <class T>
class gsGenericAssembler : public gsAssemblerBase<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:

    /// Constructor with gsMultiBasis
    gsGenericAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & bases,
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

    /// Mass assembly routine
    const gsSparseMatrix<T> & assembleMass()
    {
        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;
        m_matrix.resize(m_dofs, m_dofs); // Clean matrices
        m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );

        // Mass visitor
        gsVisitorMass<T> mass;

        for (unsigned np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble mass matrix and rhs for the local patch
            // with index np and add to m_matrix
            this->apply(mass, np);
        }

        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
        return m_matrix;
    }

    /// Stiffness assembly routine
    const gsSparseMatrix<T> & assembleStiffness()
    {
        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;
        m_matrix.resize(m_dofs, m_dofs); // Clean matrices
        m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );

        // Stiffness visitor
        gsVisitorGradGrad<T> stiffness;

        for (unsigned np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble stiffness matrix and rhs for the local patch
            // with index np and add to m_matrix
            this->apply(stiffness, np);
        }

        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
        return m_matrix;
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleMass(int patchIndex)
    {
        const int sz = m_bases.front()[patchIndex].size();

        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front()[patchIndex].degree(i) + 1;
        m_matrix.resize(sz, sz); // Clean matrix
        m_matrix.reserve( gsVector<int>::Constant(sz, nonZerosPerCol) );

        // Mass visitor
        gsVisitorMass<T> mass;

        // Set shift such that "global" indices correspond to patch-local numbering
        // (!) assumes a conforming DoF mapper
        m_dofMappers.front().setShift(-m_dofMappers.front().offset(patchIndex) );
   
        //Assemble stiffness matrix for this patch
        this->apply(mass, patchIndex);

        m_dofMappers.front().setShift(0);

        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
        return m_matrix;
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleStiffness(int patchIndex)
    {
        const int sz = m_bases.front()[patchIndex].size();

        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front()[patchIndex].degree(i) + 1;
        m_matrix.resize(sz, sz); // Clean matrix
        m_matrix.reserve( gsVector<int>::Constant(sz, nonZerosPerCol) );

        // Stiffness visitor
        gsVisitorGradGrad<T> stiffness;

        // Set shift such that "global" indices correspond to patch-local numbering
        // (!) assumes a non-conforming DoF mapper
        m_dofMappers.front().setShift(-m_dofMappers.front().offset(patchIndex) );

        //Assemble stiffness matrix for this patch
        this->apply(stiffness, patchIndex);

        m_dofMappers.front().setShift(0);

        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
        return m_matrix;
    }
    
    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() returns a lower diagonal matrix,
    /// since we exploit symmetry during assembly (whenever possible).
    typename gsSparseMatrix<T>::fullView fullMatrix()
    {
        return m_matrix.template selfadjointView<Lower>();
    }

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() returns a lower diagonal matrix,
    /// since we exploit symmetry during assembly (whenever possible).
    const typename gsSparseMatrix<T>::constFullView fullMatrix() const
    {
        return m_matrix.template selfadjointView<Lower>();
    }
    

private:

    // Members from gsAssemblerBase
    using gsAssemblerBase<T>::m_patches;
    using gsAssemblerBase<T>::m_bases;
    using gsAssemblerBase<T>::m_dofMappers;
    using gsAssemblerBase<T>::m_dofs;
    using gsAssemblerBase<T>::m_matrix;

private:
    // Hiding the rhs
    const gsMatrix<T> & rightHandSide() const;

};



} // namespace gismo

