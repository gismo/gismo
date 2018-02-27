/** @file gsParameterDomainPreconditioners.hpp

    @brief Provides preconditioners that live on the parameter domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsSolver/gsParameterDomainPreconditioners.h>
#include <gsSolver/gsSumOfOperatorsOp.h>
#include <gsSolver/gsProductOfOperatorsOp.h>
#include <gsSolver/gsKroneckerOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsMatrix/gsKronecker.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>
#include <gsTensor/gsTensorBasis.h>

namespace gismo
{

template<typename T>
void gsParameterDomainPreconditioners<T>::init()
{

    GISMO_ASSERT( m_dirichlet == dirichlet::elimination || m_dirichlet == dirichlet::none,
                  "Only the dirichlet strategies dirichlet::elimination and dirichlet::none are currently implemented." );
    // We do not provide:
    //   dirichlet::penalize
    //   dirichlet::nitsche
    //   dirichlet::eliminatNormal

}

namespace {

template <typename T>
int estimateNonzerosPerRow(const gsBasis<T>& basis)
{
    int nnz = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nnz *= 2 * basis.degree(i) + 1;
    return nnz;
}

template <typename T>
void localToGlobal(const gsMatrix<T>& localMat, const gsMatrix<unsigned>& localDofs, gsSparseMatrix<T>& globalMat)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        for (index_t j = 0; j < numActive; ++j)
        {
            const int jj = localDofs(j,0);
            globalMat.coeffRef(ii, jj) += localMat(i, j);
        }
    }
}

} // anonymous namespace

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::assembleMass(const gsBasis<T>& basis)
{
    const int n = basis.size();

    gsSparseMatrix<T> M(n, n);
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

        localToGlobal(localMass, domIt->activeFuncs, M);
    }

    M.makeCompressed();

    return M;
}

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::assembleStiffness(const gsBasis<T>& basis)
{
    const int n = basis.size();
    const int dim = basis.dim();

    gsSparseMatrix<T> K(n, n);
    K.reserve( gsVector<int>::Constant(n, estimateNonzerosPerRow(basis)) );

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( dim );
    for (int i = 0; i < basis.dim(); ++i)
        numQuadNodes[i] = basis.degree(i) + 1;
    domIt->computeQuadratureRule( numQuadNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveFunctions().rows();
        domIt->evaluateBasis( 1 );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k];

            const index_t numGrads = domIt->basisDerivs(1).rows() / dim;
            const gsAsConstMatrix<T> grads_k(domIt->basisDerivs(1).col(k).data(), dim, numGrads);

            localStiffness.noalias() += weight * (grads_k.transpose() * grads_k);
        }  // end loop Gauss nodes

        localToGlobal(localStiffness, domIt->activeFuncs, K);

    } //end loop over all domain elements

    K.makeCompressed();
    return K;
}

namespace {

template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > _assembleTensorMass(const gsTensorBasis<d,T>& basis)
{
    std::vector< gsSparseMatrix<T> > result;
    result.reserve(d);
    //for ( index_t i=0; i<d; ++i )
    for ( index_t i=d-1; i!=-1; --i )
        result.push_back( gsParameterDomainPreconditioners<T>::assembleMass( basis.component(i) ) );
    return result;
}

template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > _assembleTensorStiffness(const gsTensorBasis<d,T>& basis)
{
    std::vector< gsSparseMatrix<T> > result;
    result.reserve(d);
    //for ( index_t i=0; i<d; ++i )
    for ( index_t i=d-1; i!=-1; --i )
        result.push_back( gsParameterDomainPreconditioners<T>::assembleStiffness( basis.component(i) ) );
    return result;
}

} // anonymous namespace

template<typename T>
std::vector< gsSparseMatrix<T> > gsParameterDomainPreconditioners<T>::assembleTensorMass(const gsBasis<T>& basis)
{
    switch (basis.dim()) {
        case 1: return _assembleTensorMass<1,T>(dynamic_cast< const gsTensorBasis<1,T>& >(basis));
        case 2: return _assembleTensorMass<2,T>(dynamic_cast< const gsTensorBasis<2,T>& >(basis));
        case 3: return _assembleTensorMass<3,T>(dynamic_cast< const gsTensorBasis<3,T>& >(basis));
        case 4: return _assembleTensorMass<4,T>(dynamic_cast< const gsTensorBasis<4,T>& >(basis));
        default: GISMO_ENSURE( basis.dim() <= 4, "gsParameterDomainPreconditioners is only instanciated for up to 4 dimensions." );
    }
    return std::vector< gsSparseMatrix<T> >(); // to eliminate warning
}

template<typename T>
std::vector< gsSparseMatrix<T> > gsParameterDomainPreconditioners<T>::assembleTensorStiffness(const gsBasis<T>& basis)
{
    switch (basis.dim()) {
        case 1: return _assembleTensorStiffness<1,T>(dynamic_cast< const gsTensorBasis<1,T>& >(basis));
        case 2: return _assembleTensorStiffness<2,T>(dynamic_cast< const gsTensorBasis<2,T>& >(basis));
        case 3: return _assembleTensorStiffness<3,T>(dynamic_cast< const gsTensorBasis<3,T>& >(basis));
        case 4: return _assembleTensorStiffness<4,T>(dynamic_cast< const gsTensorBasis<4,T>& >(basis));
        default: GISMO_ENSURE( basis.dim() <= 4, "gsParameterDomainPreconditioners is only instanciated for up to 4 dimensions." );
    }
    return std::vector< gsSparseMatrix<T> >(); // to eliminate warning
}

template<typename T>
void gsParameterDomainPreconditioners<T>::handleDirichletConditions(gsSparseMatrix<T>& matrix, const gsBoundaryConditions<T>& bc, const boxSide& west, const boxSide& east)
{
    patchSide mywest(0,west), myeast(0,east);

    int i = 0;

    if (bc.getConditionFromSide( mywest )!=NULL && bc.getConditionFromSide( mywest )->type() == condition_type::dirichlet ) i += 1;
    if (bc.getConditionFromSide( myeast )!=NULL && bc.getConditionFromSide( myeast )->type() == condition_type::dirichlet ) i += 2;

    switch ( i )
    {
        case 0: return;
        case 1: matrix = matrix.block( 1, 1, matrix.rows()-1, matrix.cols()-1 ); return;
        case 2: matrix = matrix.block( 0, 0, matrix.rows()-1, matrix.cols()-1 ); return;
        case 3: matrix = matrix.block( 1, 1, matrix.rows()-2, matrix.cols()-2 ); return;
    }
}

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::getMassMatrix() const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
        for (index_t i=0; i<d; ++i)
            handleDirichletConditions(local_mass[i],m_bc,1+2*i,2+2*i);
    return getKroneckerProduct(local_mass);
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getMassMatrixOp() const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
        for (index_t i=0; i<d; ++i)
            handleDirichletConditions(local_mass[i],m_bc,1+2*i,2+2*i);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeMatrixOp(local_mass[i].moveToPtr());

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getMassMatrixInvOp() const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
        for (index_t i=0; i<d; ++i)
            handleDirichletConditions(local_mass[i],m_bc,1+2*i,2+2*i);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeSparseCholeskySolver(local_mass[i]);

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::getStiffnessMatrix(T a) const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
    {
        for (index_t i=0; i<d; ++i)
        {
            handleDirichletConditions(local_stiff[i],m_bc,1+2*i,2+2*i);
            handleDirichletConditions(local_mass [i],m_bc,1+2*i,2+2*i);
        }
    }
    gsSparseMatrix<T> K = give(local_stiff[0]);
    gsSparseMatrix<T> M = give(local_mass [0]);

    for (index_t i=1; i<d; ++i)
    {
        K  = getKroneckerProduct(K, local_mass[i]);
        K += getKroneckerProduct(M, local_stiff[i]);
        if ( i < d-1 || a != 0 )
            M = getKroneckerProduct(M, local_mass[i]);
    }
    if (a==1)
        K += M;
    else if (a!=0)
        K += a * M;
    return K;
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getStiffnessMatrixOp(T a) const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
    {
        for (index_t i=0; i<d; ++i)
        {
            handleDirichletConditions(local_stiff[i],m_bc,1+2*i,2+2*i);
            handleDirichletConditions(local_mass [i],m_bc,1+2*i,2+2*i);
        }
    }
    std::vector<OpUPtr> local_stiff_op(d);
    std::vector<OpPtr > local_mass_op (d);
    for (index_t i=0; i<d; ++i)
    {
        local_stiff_op[i] = makeMatrixOp(local_stiff[i].moveToPtr());
        local_mass_op [i] = makeMatrixOp(local_mass [i].moveToPtr());
    }
    OpUPtr K = give(local_stiff_op[0]);
    OpPtr  M = give(local_mass_op [0]);

    for (index_t i=1; i<d; ++i)
    {
        K = gsSumOfOperatorsOp<T>::make(
            gsKroneckerOp<T>::make( give(K),      local_mass_op [i]  ),
            gsKroneckerOp<T>::make( M,       give(local_stiff_op[i]) )
        );
        if ( i < d-1 || a != 0 )
            M = gsKroneckerOp<T>::make(M, local_mass_op[i]);
    }
    if (a==1)
        K = gsSumOfOperatorsOp<T>::make( give(K), M );
    else if (a!=0)
        K = gsSumOfOperatorsOp<T>::make( give(K), gsScaledOp<T>::make( M, a ) );

    return K;
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getFastDiagonalizationOp(T a) const
{

    const index_t d = m_basis.dim();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis);
    if (m_dirichlet == dirichlet::elimination)
    {
        for (index_t i=0; i<d; ++i)
        {
            handleDirichletConditions(local_stiff[i],m_bc,1+2*i,2+2*i);
            handleDirichletConditions(local_mass [i],m_bc,1+2*i,2+2*i);
        }
    }

    // Determine overall size
    index_t sz = 1;
    for ( index_t i=0; i<d; ++i )
        sz *= local_stiff[i].rows();

    // Initialize the diagonal with 1
    gsMatrix<T> diag;
    diag.setConstant(sz,1,a); // This is the pure-mass part!

    index_t glob = sz; // Indexing value for setting up the Kronecker product

    typedef typename gsMatrix<T>::GenSelfAdjEigenSolver EVSolver;
    typedef typename EVSolver::EigenvectorsType evMatrix;
    typedef typename EVSolver::RealVectorType evVector;
    EVSolver ges;

    std::vector<OpPtr> Qop(d);
    std::vector<OpPtr> QTop(d);

    // Now, setup the Q's and update the D's
    for ( index_t i=0; i<d; ++i )
    {
        // Solve generalized eigenvalue problem
        ges.compute(local_stiff[i], local_mass [i], Eigen::ComputeEigenvectors);
        // Q^T M Q = I, or M = Q^{-T} Q^{-1}
        // Q^T K Q = D, or K = Q^{-T} D Q^{-1}

        // From the eigenvalues, we setup the matrix D already in an Kroneckerized way.
        const evVector & D = ges.eigenvalues();

        const index_t loc = D.rows();
        glob /= loc;
        const index_t glob2 = sz / loc / glob;

        for ( index_t l=0; l<loc; ++l )
            for ( index_t m=0; m<glob; ++m )
                for ( index_t n=0; n<glob2; ++n )
                    diag( m + l*glob + n*loc*glob, 0 ) += D(l,0);

        // Finally, we store the eigenvectors
        gsMatrix<T> ev;
        ev.swap(const_cast<evMatrix&>(ges.eigenvectors()));

        // These are the operators representing the eigenvectors
        typename gsMatrixOp<  gsMatrix<T> >::Ptr matrOp = makeMatrixOp( ev.moveToPtr() );
        Qop [i] = matrOp;
        // Here we are safe as long as we do not want to apply QTop after Qop got destroyed.
        QTop[i] = makeMatrixOp( matrOp->matrix().transpose() );
    }

    GISMO_ASSERT( glob == 1, "Internal error." );

    for ( index_t l=0; l<sz; ++l )
        diag( l, 0 ) = 1/diag( l, 0 );

    memory::unique_ptr< Eigen::DiagonalMatrix<real_t,Dynamic> > diag_mat( new Eigen::DiagonalMatrix<real_t,Dynamic>( give(diag) ) );

    return gsProductOfOperatorsOp<T>::make(
        gsKroneckerOp<T>::make(QTop),
        makeMatrixOp(give(diag_mat)),
        gsKroneckerOp<T>::make(Qop)
    );
}

//Will be provided in a followup pull request
//template<typename T>
//typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getSubspaceCorrectedMassSmootherOp() const
//{
//
//}

} // namespace gismo
