/** @file gsPatchPreconditionersCreator.hpp

    @brief Provides preconditioners that live on the parameter domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsKroneckerOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

namespace {

template<typename T>
gsBoundaryConditions<T> boundaryConditionsForDirection( const gsBoundaryConditions<T>& bc, index_t direction )
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

template<typename T>
void eliminateDirichlet1D(const gsBoundaryConditions<T>& bc,
                          const gsOptionList& opt, gsSparseMatrix<T> & result)
{
    dirichlet::strategy ds = (dirichlet::strategy)opt.askInt("DirichletStrategy",dirichlet::elimination);
    if (ds == dirichlet::elimination)
    {
        patchSide west(0,1), east(0,2);
        index_t i = 0;
        if (bc.getConditionFromSide( west )!=NULL && bc.getConditionFromSide( west )->type() == condition_type::dirichlet ) i += 1;
        if (bc.getConditionFromSide( east )!=NULL && bc.getConditionFromSide( east )->type() == condition_type::dirichlet ) i += 2;
        if (i%2 + i/2 >= result.rows() || i%2 + i/2 >= result.cols())
            result.resize(0,0);
        else switch ( i )
             {
             case 0: break;
             case 1: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
             case 2: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
             case 3: result = result.block( 1, 1, result.rows()-2, result.cols()-2 ); break;
             }
    }
    else
        GISMO_ERROR("Unknown Dirichlet strategy.");
}

template<typename T>
gsSparseMatrix<T> assembleMass(const gsBasis<T>& basis)
{
    gsExprAssembler<T> mass(1,1);
    gsMultiBasis<T> mb(basis);
    mass.setIntegrationElements(mb);
    typename gsExprAssembler<T>::space u = mass.getSpace(mb);
    mass.initMatrix();
    mass.assemble( u * u.tr() );
    gsSparseMatrix<T> result;
    mass.matrix_into(result);
    return result;
}

template<typename T>
gsSparseMatrix<T> assembleStiffness(const gsBasis<T>& basis)
{
    gsExprAssembler<T> stiff(1,1);
    gsMultiBasis<T> mb(basis);
    stiff.setIntegrationElements(mb);
    typename gsExprAssembler<T>::space u = stiff.getSpace(mb);
    stiff.initMatrix();
    stiff.assemble( grad(u) * grad(u).tr() );
    gsSparseMatrix<T> result;
    stiff.matrix_into(result);
    return result;
}


template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > assembleTensorMass_impl(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const gsTensorBasis<d,T> * tb = dynamic_cast< const gsTensorBasis<d,T>* >( &basis );
    GISMO_ENSURE (tb, "gsPatchPreconditionersCreator requires a tensor basis.");

    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = assembleMass(tb->component(d-i-1));
        eliminateDirichlet1D(boundaryConditionsForDirection(bc,i), opt, result[i]);
    }
    return result;
}

template<index_t d, typename T>
std::vector< gsSparseMatrix<T> > assembleTensorStiffness_impl(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const gsTensorBasis<d,T> * tb = dynamic_cast< const gsTensorBasis<d,T>* >( &basis );
    GISMO_ENSURE (tb, "gsPatchPreconditionersCreator requires a tensor basis.");

    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = assembleStiffness(tb->component(d-i-1));
        eliminateDirichlet1D(boundaryConditionsForDirection(bc,i), opt, result[i]);
    }
    return result;
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorMass(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    switch (basis.dim())
    {
    case 1: return assembleTensorMass_impl<1,T>(basis, bc, opt);
    case 2: return assembleTensorMass_impl<2,T>(basis, bc, opt);
    case 3: return assembleTensorMass_impl<3,T>(basis, bc, opt);
    case 4: return assembleTensorMass_impl<4,T>(basis, bc, opt);
    default: GISMO_ERROR("gsPatchPreconditionersCreator is only instanciated for up to 4 dimensions." );
    }
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorStiffness(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    switch (basis.dim())
    {
    case 1: return assembleTensorStiffness_impl<1,T>(basis, bc, opt);
    case 2: return assembleTensorStiffness_impl<2,T>(basis, bc, opt);
    case 3: return assembleTensorStiffness_impl<3,T>(basis, bc, opt);
    case 4: return assembleTensorStiffness_impl<4,T>(basis, bc, opt);
    default: GISMO_ERROR("gsPatchPreconditionersCreator is only instanciated for up to 4 dimensions." );
    }
}

} // anonymous namespace

template<typename T>
gsSparseMatrix<T> gsPatchPreconditionersCreator<T>::massMatrix(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);
    const index_t d = local_mass.size();
    gsSparseMatrix<T> result = local_mass[d-1];
    for (index_t i=d-2; i>-1; --i)
        result = result.kron(local_mass[i]);
    return result;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::massMatrixOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeMatrixOp(local_mass[i].moveToPtr());

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::massMatrixInvOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeSparseCholeskySolver(local_mass[i]);

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
gsSparseMatrix<T> gsPatchPreconditionersCreator<T>::stiffnessMatrix(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T a
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    gsSparseMatrix<T> K = give(local_stiff[d-1]);
    gsSparseMatrix<T> M = give(local_mass [d-1]);

    for (index_t i=d-2; i>-1; --i)
    {
        K  = K.kron(local_mass[i]);
        K += M.kron(local_stiff[i]);
        if ( i != 0 || a != 0 )
            M = M.kron(local_mass[i]);
    }
    if (a==1)
        K += M;
    else if (a!=0)
        K += a * M;
    return K;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::stiffnessMatrixOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T a
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

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
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T a
    )
{

    const index_t d = basis.dim();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

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
    gsMatrix<T> ev;

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

    memory::unique_ptr< Eigen::DiagonalMatrix<T,Dynamic> > diag_mat( new Eigen::DiagonalMatrix<T,Dynamic>( give(diag) ) );

    return gsProductOp<T>::make(
        gsKroneckerOp<T>::make(QTop),
        makeMatrixOp(give(diag_mat)),
        gsKroneckerOp<T>::make(Qop)
        );
}

namespace {

// Get the tilde basis
template<typename T>
void tildeSpaceBasis_oneside(const gsTensorBSplineBasis<1,T>& basis, bool isLeftHandSide, gsMatrix<T>& tildeBasis, gsMatrix<T>& complBasis, index_t b = 0, const bool odd = true)
{
    // b == 0: Neumann (or any other not-eliminating bc), b == 1: Dirichlet

    GISMO_ASSERT ( b == 0 || b == 1, "Not a feasible boundary condition specified." );

    const int p = basis.degree();
    const T h = basis.knots().maxIntervalLength();

    if (p-b<=0)
    {
        tildeBasis.resize(0,0); complBasis.resize(0,0);
        return;
    }

    const T u = (isLeftHandSide ? basis.knots().first() : basis.knots().last());
    gsMatrix<T> U(1,1);
    U(0,0) = u;

    std::vector< gsMatrix<T> > allDerivs;
    basis.evalAllDers_into(U, p-1, allDerivs);

    // Collect all derivatives in matrix (rows: derivatives, columns: basis functions)
    // normalize with h^(deriv)
    // use only odd derivatives if odd = true and only even derivatives if odd = false
    // skip last (on left) or first (on right) basis function since it's always in S-tilde
    gsMatrix<T> derivs;
    // If we have Dirichlet bc, we reduce the number of rows and columns by 1. We basically
    // eliminate the first row and the first column (first basis function, 0th derivative)
    // for the right-side, we have to remove the last basis function.
    derivs.setZero(p-b, p-b);

    const int offset = (isLeftHandSide ? b : 1);

    int i_start;
    if (odd) i_start = 1;
    else i_start = 2*b;

    for (int i = i_start; i < p; i += 2)
        for (int j = 0; j < p-b; ++j)
            derivs(i-b, j) = math::pow(h, i) * allDerivs[i](j+offset);

    typename gsMatrix<T>::JacobiSVD svd = derivs.jacobiSvd(Eigen::ComputeFullV);

    int n_tilde;
    if (odd) n_tilde = (p + 1) / 2 - b;
    else n_tilde = p / 2 - b;

    tildeBasis = svd.matrixV().rightCols(n_tilde);
    complBasis = svd.matrixV().leftCols(p - b - n_tilde);
}

// Compute a basis for S-tilde and one for its orthogonal complement
template<typename T>
void tildeSpaceBasis(const gsTensorBSplineBasis<1,T>& basis, gsSparseMatrix<T>& B_tilde, gsSparseMatrix<T>& B_compl, index_t b, const bool odd = true)
{
    // Boundary conditions left-hand side b%2, right-hand side b/2:
    // 0: Neumann (or any other not-eliminating bc), 1: Dirichlet

    GISMO_ASSERT ( b >= 0 && b < 4, "tildeSpaceBasis: Not a feasible boundary condition specified." );

    gsMatrix<T> b_L, b_compl_L;
    gsMatrix<T> b_R, b_compl_R;

    // Contruct space with vanishing odd derivatives
    tildeSpaceBasis_oneside(basis, true,  b_L, b_compl_L, b%2, odd);
    tildeSpaceBasis_oneside(basis, false, b_R, b_compl_R, b/2, odd);

    const int n = basis.size() - b%2 - b/2;
    const int n_L = b_L.cols();
    const int m_L = b_L.rows();
    const int n_R = b_R.cols();
    const int m_R = b_R.rows();
    const int n_c_L = b_compl_L.cols();
    const int m_c_L = b_compl_L.rows();
    const int n_c_R = b_compl_R.cols();
    const int m_c_R = b_compl_R.rows();
    const int n_I = n - n_L - n_R - n_c_L - n_c_R;

    //GISMO_ENSURE ( n_I >= 0, "tildeSpaceBasis: Too few knots for that spline degree." );
    if ( n_I <= 0 )
    {
        static bool warned = false;
        if (!warned)
        {
            gsWarn << "tildeSpaceBasis was called with too few knots for that spline degree.\n"
                   << "So, we assume that S_tilde is empty.\n";
            warned = true;
        }

        B_tilde.resize(n,0);

        gsSparseEntries<T> E_compl;
        E_compl.reserve(n);
        for (int i = 0; i < n; ++i)
            E_compl.add(i,i,1.0);
        B_compl.resize(n,n);
        B_compl.setFrom(E_compl);
        return;
    }

    gsSparseEntries<T> E_tilde, E_compl;

    // put b_L into upper left block of S-tilde basis
    for (int j = 0; j < n_L; ++j)
    {
        for (int i = 0; i < m_L; ++i)
            E_tilde.add(i, j, b_L(i, j));
    }
    // fill identity matrix into interior part of S-tilde basis
    for (int j = 0; j < n_I; ++j)
    {
        E_tilde.add(m_L + j, n_L + j, 1.0);
    }
    // put b_R into lower right block of S-tilde basis
    for (int j = 0; j < n_R; ++j)
    {
        for (int i = 0; i < m_R; ++i)
            E_tilde.add(m_L + n_I + i, n_L + n_I + j, b_R(i, j));
    }
    B_tilde.resize(n, n_L + n_I + n_R);
    B_tilde.setFrom(E_tilde);

    // put b_compl_L into upper left block of complement basis
    for (int j = 0; j < n_c_L; ++j)
    {
        for (int i = 0; i < m_c_L; ++i)
            E_compl.add(i, j, b_compl_L(i, j));
    }
    // put b_compl_R into lower right block of complement basis
    for (int j = 0; j < n_c_R; ++j)
    {
        for (int i = 0; i < m_c_R; ++i)
            E_compl.add(m_c_L + n_I + i, n_c_L + j, b_compl_R(i, j));
    }
    B_compl.resize(n, n_c_L + n_c_R);
    B_compl.setFrom(E_compl);

}

template<index_t d, typename T>
void constructTildeSpaceBasisTensor(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    std::vector< gsSparseMatrix<T> >& B_tilde,
    std::vector< gsSparseMatrix<T> >& B_l2compl,
    const bool odd = true
    )
{

    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This method only works for gsTensorBSplineBasis.");

    B_tilde.resize(d);
    B_l2compl.resize(d);
    dirichlet::strategy ds = (dirichlet::strategy)opt.askInt("DirichletStrategy",dirichlet::elimination);
    if (ds == dirichlet::elimination)
    {
        for ( index_t i=0; i<d; ++i )
        {
            index_t b = 0;
            // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
            patchSide mywest(0,1+2*i);
            patchSide myeast(0,2+2*i);
            if ( bc.getConditionFromSide( mywest )!=NULL && bc.getConditionFromSide( mywest )->type() == condition_type::dirichlet ) b += 1;
            if ( bc.getConditionFromSide( myeast )!=NULL && bc.getConditionFromSide( myeast )->type() == condition_type::dirichlet ) b += 2;

            tildeSpaceBasis(tb->component(i), B_tilde[i], B_l2compl[i], b, odd);
        }
    }
    else
        GISMO_ERROR("Unknown Dirichlet strategy.");

}

template<typename T>
void constructTildeSpaceBasis(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    std::vector< gsSparseMatrix<T> >& B_tilde,
    std::vector< gsSparseMatrix<T> >& B_l2compl,
    const bool odd = true
    )
{
    switch (basis.dim())
    {
    case 1: constructTildeSpaceBasisTensor<1>( basis, bc, opt, B_tilde, B_l2compl, odd ); return;
    case 2: constructTildeSpaceBasisTensor<2>( basis, bc, opt, B_tilde, B_l2compl, odd ); return;
    case 3: constructTildeSpaceBasisTensor<3>( basis, bc, opt, B_tilde, B_l2compl, odd ); return;
    case 4: constructTildeSpaceBasisTensor<4>( basis, bc, opt, B_tilde, B_l2compl, odd ); return;
    default: GISMO_ERROR ("gsPatchPreconditionersCreator is only instanciated for up to 4 dimensions.");
    }

}

// Constructs a matrix for swapping a tensor product
// from A (x) B (x) C (x) D (x) E  to   A (x) D (x) C (x) B (x) E,
// where only the dimensions of those matrices have to be given.
// Note that also those A, B, etc. could be tensor products (here as dimension just provide the products)
// Note that also those A, B, etc. could vanish. Then just provide a 1 as dimension, i.e., a scalar.
// So, literally every thinkable swap is possible.
template<typename T>
gsSparseMatrix<T> kroneckerSwap( index_t e, index_t d, index_t c, index_t b, index_t a )
{
    const index_t sz = a*b*c*d*e;
    gsSparseMatrix<T> result(sz,sz);
    gsSparseEntries<T> entries;
    entries.reserve(sz);
    for ( index_t i=0; i<a; ++i )
        for ( index_t j=0; j<b; ++j )
            for ( index_t k=0; k<c; ++k )
                for ( index_t l=0; l<d; ++l )
                    for ( index_t m=0; m<e; ++m )
                        entries.add( i+a*(j+b*(k+c*(l+d*m))), i+a*(l+d*(k+c*(j+b*m))), 1. );

    result.setFrom(entries);
    result.makeCompressed();

    return result;
}

} // anonymous namespace

template<typename T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< gsSparseMatrix<T> > > gsPatchPreconditionersCreator<T>::getTildeSpaceBasisTransformation(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    std::vector< gsSparseMatrix<T> > B_tilde, B_l2compl;
    constructTildeSpaceBasis( basis, bc, opt, B_tilde, B_l2compl, true );
    return std::pair< std::vector< gsSparseMatrix<T> >, std::vector< gsSparseMatrix<T> > >( B_tilde, B_l2compl );
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::subspaceCorrectedMassSmootherOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T sigma,
    T alpha
    )
{
    // Get some properties
    const index_t d = basis.dim();
    const T  h = basis.getMinCellLength();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    // Setup of basis
    std::vector< gsSparseMatrix<T> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, opt, B_tilde, B_l2compl);

    std::vector< gsSparseMatrix<T> > M_compl(d), K_compl(d);
    std::vector< typename gsLinearOperator<T>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // Transform the complement
        typename gsLinearOperator<T>::Ptr M_inv = makeSparseCholeskySolver(local_mass[i]);
        gsMatrix<T> B_compl_dense;
        M_inv->apply( B_l2compl[i], B_compl_dense );
        B_compl[i] = B_compl_dense.sparseView();

        // Setup of matrices and corresponding solvers
        gsSparseMatrix<T> M_tilde = B_tilde[i].transpose() * local_mass[i] * B_tilde[i];
        M_tilde_inv[i] = makeSparseCholeskySolver( M_tilde );
        M_compl[i] = B_compl[i].transpose() * local_mass[i] * B_compl[i];
        K_compl[i] = B_compl[i].transpose() * local_stiff[i] * B_compl[i];
    }

    // Setup of final operator
    typename gsSumOp<T>::uPtr result = gsSumOp<T>::make();

    for ( index_t type = 0; type < (1<<d); ++ type )
    {
        std::vector< typename gsLinearOperator<T>::Ptr > correction(0);
        gsSparseMatrix<T> transfer;

        std::vector< gsSparseMatrix<T>* > transfers(d);

        index_t numberInteriors = 0;

        // Setup of transfer
        for ( index_t j = 0; j<d; ++ j )
        {
            if ( type & ( 1 << j ) )
                transfers[j] = &(B_compl[j]);
            else
            {
                transfers[j] = &(B_tilde[j]);
                ++numberInteriors;
            }

            if ( j == 0 )
                transfer = *(transfers[j]);
            else
                transfer = transfers[j]->kron(transfer);

        }

        // If the subspace is not present, ignore it.
        if ( transfer.cols() == 0 )
            continue;

        // Setup of swap, where the boundary part is shifted to the begin

        index_t left = 1, current = transfers[d-1]->cols(), right = 1;
        for ( index_t j = 0; j < d-1; ++j )
            left *= transfers[j]->cols();

        for ( index_t j = d-1; j >= 0; --j )
        {
            if ( type & ( 1 << j ) )
            {
                transfer = transfer * kroneckerSwap<T>( right, current, left, 1, 1 );
                if ( j > 0 )
                {
                    left *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }
            else
            {
                if ( j > 0 )
                {
                    right *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }

        }

        // Setup of interior correction
        for ( index_t j = d-1; j>=0; --j )
        {
            if ( ! ( type & ( 1 << j ) ) )
                correction.push_back( M_tilde_inv[j] );
        }

        // If we are in the interior, we have to do the scaling here as there is no boundary correction.
        if ( numberInteriors == d )
            correction[0] = gsScaledOp<T>::make( correction[0], 1./( alpha + numberInteriors/(sigma*h*h) ) );

        // Setup of bondary correction, like  K(x)M(x)M+M(x)K(x)M+M(x)M(x)K+(alpha + (d-3)/(sigma*h*h)) M(x)M(x)M
        if ( numberInteriors < d )
        {
            gsSparseMatrix<T> bc_matrix;

            {
                gsSparseMatrix<T> s(0,0);
                for ( index_t k = d-1; k>=0; --k )
                {
                    if ( type & ( 1 << k ) )
                    {
                        if ( s.rows() == 0 )
                            s = M_compl[k];
                        else
                            s = M_compl[k].kron(s);
                    }
                }
                bc_matrix = ( alpha + numberInteriors/(sigma*h*h) ) * s;
            }

            for ( index_t j = d-1; j>=0; --j )
            {
                if ( type & ( 1 << j ) )
                {
                    gsSparseMatrix<T> s(0,0);
                    for ( index_t k = d-1; k>=0; --k )
                    {
                        if ( type & ( 1 << k ) )
                        {
                            gsSparseMatrix<T>* chosenMatrix;
                            if ( j == k )
                                chosenMatrix = &(K_compl[k]);
                            else
                                chosenMatrix = &(M_compl[k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                                s = s.kron(*chosenMatrix);
                        }
                    }
                    bc_matrix += s;
                }
            }

            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }

        typename gsMatrixOp< gsSparseMatrix<T> >::Ptr transOp = makeMatrixOp(transfer.moveToPtr());
        // Setup of whole operator
        // The correction is the Kronecker-product of the operators in the vector correction.
        result->addOperator(
            gsProductOp<T>::make(
                makeMatrixOp( transOp->matrix().transpose() ),
                gsKroneckerOp<T>::make( correction ),
                transOp
                )
            );
    }

    return give(result);

}

} // namespace gismo
