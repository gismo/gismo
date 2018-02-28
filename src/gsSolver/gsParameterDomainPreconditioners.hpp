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
#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsKroneckerOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsTensor/gsTensorTools.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>
#include <gsTensor/gsTensorBasis.h>
#include <gsAssembler/gsGenericAssembler.h>

namespace gismo
{

namespace {

template<typename T>
gsBoundaryConditions<T> getBoundaryConditionsForDirection( const gsBoundaryConditions<T>& bc, index_t direction )
{
    gsBoundaryConditions<T> result;

    for ( index_t i = 1; i <= 2; ++i)
    {
        patchSide global(0,i+2*direction), local(0,i);
        const boundary_condition<T>* cond = bc.getConditionFromSide(global);
        if (cond!=NULL)
            result.addCondition(local,cond->type(),cond->function());
    }
    return result;
}

template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > _assembleTensorMass(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& options
)
{
    std::vector< gsSparseMatrix<T> > result;
    const gsTensorBasis<d,T> * tb = dynamic_cast< const gsTensorBasis<d,T>* >( &basis );

    if (tb==nullptr)
    {
        gsWarn << "gsParameterDomainPreconditioners: Found a discretization which does not have "
            "tensor-product structure. Therefore, the preconditioner might not be efficient." << std::endl;
        gsGenericAssembler<T> assembler(gsMultiPatch<T>(),gsMultiBasis<T>(basis),options,&bc);
        result.push_back( assembler.assembleMass() );
        return result;
    }
    else
    {
        result.reserve(d);
        for ( index_t i=d-1; i!=-1; --i )
        {
            gsBoundaryConditions<T> local_bc = getBoundaryConditionsForDirection(bc,i);
            gsGenericAssembler<T> assembler(gsMultiPatch<T>(),gsMultiBasis<T>(tb->component(i)),options,&local_bc);
            result.push_back( assembler.assembleMass() );
        }
        return result;
    }
}

template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > _assembleTensorStiffness(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& options
)
{
    std::vector< gsSparseMatrix<T> > result;
    const gsTensorBasis<d,T> * tb = dynamic_cast< const gsTensorBasis<d,T>* >( &basis );

    if (tb==nullptr)
    {
        gsWarn << "gsParameterDomainPreconditioners: Found a discretization which does not have "
            "tensor-product structure. Therefore, the preconditioner might not be efficient." << std::endl;
        gsGenericAssembler<T> assembler(gsMultiPatch<T>(),gsMultiBasis<T>(basis),options,&bc);
        result.push_back( assembler.assembleStiffness() );
        return result;
    }
    else
    {
        result.reserve(d);
        for ( index_t i=d-1; i!=-1; --i )
        {
            gsBoundaryConditions<T> local_bc = getBoundaryConditionsForDirection(bc,i);
            gsGenericAssembler<T> assembler(
                gsMultiPatch<T>(),
                gsMultiBasis<T>(tb->component(i)),
                options,
                &local_bc
            );
            result.push_back( assembler.assembleStiffness() );
        }
        return result;
    }
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorMass(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& options
)
{
    switch (basis.dim()) {
        case 1: return _assembleTensorMass<1,T>(basis, bc, options);
        case 2: return _assembleTensorMass<2,T>(basis, bc, options);
        case 3: return _assembleTensorMass<3,T>(basis, bc, options);
        case 4: return _assembleTensorMass<4,T>(basis, bc, options);
        default: GISMO_ENSURE( basis.dim() <= 4, "gsParameterDomainPreconditioners is only instanciated for up to 4 dimensions." );
    }
    return std::vector< gsSparseMatrix<T> >(); // to eliminate warning
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorStiffness(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& options
)
{
    switch (basis.dim()) {
        case 1: return _assembleTensorStiffness<1,T>(basis, bc, options);
        case 2: return _assembleTensorStiffness<2,T>(basis, bc, options);
        case 3: return _assembleTensorStiffness<3,T>(basis, bc, options);
        case 4: return _assembleTensorStiffness<4,T>(basis, bc, options);
        default: GISMO_ENSURE( basis.dim() <= 4, "gsParameterDomainPreconditioners is only instanciated for up to 4 dimensions." );
    }
    return std::vector< gsSparseMatrix<T> >(); // to eliminate warning
}

} // anonymous namespace

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::getMassMatrix() const
{
    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis, m_bc, m_options);
    return getKroneckerProduct(local_mass);
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getMassMatrixOp() const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis, m_bc, m_options);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeMatrixOp(local_mass[i].moveToPtr());

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getMassMatrixInvOp() const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(m_basis, m_bc, m_options);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeSparseCholeskySolver(local_mass[i]);

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
gsSparseMatrix<T> gsParameterDomainPreconditioners<T>::getStiffnessMatrix(T a) const
{
    const index_t d = m_basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis, m_bc, m_options);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis, m_bc, m_options);

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

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis, m_bc, m_options);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis, m_bc, m_options);

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
        K = gsSumOp<T>::make(
            gsKroneckerOp<T>::make( give(K),      local_mass_op [i]  ),
            gsKroneckerOp<T>::make( M,       give(local_stiff_op[i]) )
        );
        if ( i < d-1 || a != 0 )
            M = gsKroneckerOp<T>::make(M, local_mass_op[i]);
    }
    if (a==1)
        K = gsSumOp<T>::make( give(K), M );
    else if (a!=0)
        K = gsSumOp<T>::make( give(K), gsScaledOp<T>::make( M, a ) );

    return K;
}

template<typename T>
typename gsParameterDomainPreconditioners<T>::OpUPtr gsParameterDomainPreconditioners<T>::getFastDiagonalizationOp(T a) const
{

    const index_t d = m_basis.dim();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(m_basis, m_bc, m_options);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(m_basis, m_bc, m_options);

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

    return gsProductOp<T>::make(
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
