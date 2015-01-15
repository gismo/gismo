#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <vector>

namespace gismo {
/**
    \brief Modify the SparseMatrix in place to impose Dirichlet type conditions.

    Change the provided matrix in place in order to impose Dirichlet conditions,
    i.e. it replaces the coefficient (i,j) with
         0 if \f$i\neq j\f$ and i or j is a fixed degree of freedom;
         1 if \f$i=j\f$ and i is a fixed degree of freedom.
    It leaves the other coefficients unchanged.

    While replacing the entries it also compute the required modification to the
    right-hand-side matrix and applies them.

    \param[in,out] M                the system matrix
    \param[in,out] Rhs              the right-hand-side matrix
    \param[in]     dirichletDofs    the indexes of the Dirichlet DoFs
    \param[in]     dirichletValues  the coefficients for the fixed degrees of freedom

    \tparam        T                the type of the coefficients
    \tparam        Major            the storage class of the matrices
    \tparam        RhsM             the type of the right-hand-side matrix
    \tparam        IndexM           the type of the Dirichlet indexes parameter
    
    \ingroup Assembler
*/
template<typename T, int Major, typename RhsM, typename IndexM, typename DirichletM >
void gsForceDirichletConditions (
        gsSparseMatrix<T,Major>      &M,
        RhsM                         &Rhs,
        const IndexM                 &dirichletDofs,
        const DirichletM             &dirichletValues)
{
    if (Rhs.cols()!=dirichletValues.cols())
    {
        GISMO_ERROR("The matrix for the RHS and for the DirichletValues must have the same size");
    }
    if (dirichletDofs.cols()!=1)
    {
        GISMO_ERROR("the list of Dirichlet DoFs must be a column vector");
    }

    int max_d_idx = dirichletDofs.size ();
    int outer_idx = 0;
    int inner_idx = 0;

    for (unsigned k=0; k< (unsigned long) M.outerSize(); ++k)
    {
        if ( k> dirichletDofs(outer_idx) && outer_idx<max_d_idx-1)
            outer_idx++;
        if  (k == (unsigned long) dirichletDofs(outer_idx)) // remove everything
        {
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                if (Major==Eigen::ColMajor)
                    Rhs(it.row())-= it.value()*dirichletValues(outer_idx);
                it.valueRef()           = 0;
            }
        }
        else                              // remove entries corresponding to Dirichlet dofs
        {
            inner_idx=0;
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                while ( (unsigned long)it.index() >  dirichletDofs(inner_idx) && inner_idx<max_d_idx-1)
                    inner_idx++;
                if ((unsigned long) it.index()==  dirichletDofs(inner_idx))
                {
                    if (Major==Eigen::RowMajor)
                        Rhs(it.row())-= it.value()*dirichletValues(inner_idx);
                    it.valueRef()=0;
                }
            }
        }
    }
    // set the rhs on the fixed DoFs, and the unit diagonal of the matrix
    // we can not set the diagonal in the previous loop because it only loop
    // over non zero coefficients
    for (int k=0; k<dirichletDofs.size (); ++k)
    {
        int pos = dirichletDofs(k);
        M(pos,pos)=1;
        Rhs.row(pos).array()=dirichletValues.row(k);
    }
    return;
}



/**
    \brief Modify the SparseMatrix in place to impose Dirichlet type conditions.

    Change the provided matrix in place in order to impose Dirichlet conditions,
    i.e. it replaces the coefficient (i,j) with
        0 if \f$i\neq j\f$ and i or j is a fixed degree of freedom;
        1 if \f$i=j\f$ and i is a fixed degree of freedom.
    It leaves the other coefficients unchanged.

    While replacing the entries it also compute the matrix representing the mapping
    between the coefficients of the fixed degrees of freedom and the corresponding
    modification of the right-hand-side.

    \param[in,out] M                the system matrix
    \param[out]    valueToRhs       the matrix that maps Dirichlet coef to Rhs modifications
    \param[in]     dirichletDofs    the indexes of the Dirichlet DoFs

    \tparam        T                the type of the coefficients
    \tparam        Major            the storage class of the matrices
    \tparam        IndexM           the type of the Dirichlet indexes parameter
    
    \ingroup Assembler
*/
template<typename T, int Major, typename RhsM, typename IndexM >
void gsForceDirichletConditions (
        gsSparseMatrix<T,Major>      &M,
        gsSparseMatrix<T,Major>      &valueToRhs,
        const IndexM                  dirichletDofs)
{
    if (dirichletDofs.cols()!=1)
    {
        GISMO_ERROR("the list of Dirichlet DoFs must be a column vector");
    }

    int max_d_idx = dirichletDofs.size ();
    int outer_idx = 0;
    int inner_idx = 0;

    for (unsigned k=0; k< M.outerSize(); ++k)
    {
        if ( k> dirichletDofs(outer_idx) && outer_idx<max_d_idx-1)
            outer_idx++;
        if  (k == dirichletDofs(outer_idx)) // remove everything
        {
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                if (Major==Eigen::ColMajor)
                    valueToRhs(it.row(),outer_idx)= it.value();
                it.valueRef()           = 0;
            }
        }
        else                              // remove entries corresponding to Dirichlet DoFs
        {
            inner_idx=0;
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                while ( it.index() > dirichletDofs(inner_idx) && inner_idx<max_d_idx-1)
                    inner_idx++;
                if (it.index()== dirichletDofs(inner_idx))
                {
                    if (Major==Eigen::RowMajor)
                        valueToRhs(it.row(),inner_idx)= it.value();
                    it.valueRef()=0;
                }
            }
        }
    }
    // set the rhs on the fixed DoFs, and the unit diagonal of the matrix
    // we can not set the diagonal in the previous loop because it only loop
    // over non zero coefficients
    for (int k=0; k<dirichletDofs.size (); ++k)
    {
        int pos = dirichletDofs(k);
        M(pos,pos)=1;
    }
    return;
}



/**
    \brief Write the Dirichlet data to a right-hand-side matrix.

    Change in place the coefficients of the matrix corresponding to
    fixed degrees of freedom to the assigned values.

    \param[in,out] destination       the matrix to modify
    \param[in]     dirichletDofs     the indexes to set
    \param[in]     dirichletValues   the values
    
    \ingroup Assembler
 */
template <typename destM, typename coefM, typename indexM>
void gsSetDirichletValues (
        destM             &destination,
        const indexM      &dirichletDofs,
        const coefM       &dirichletValues)
{
    for (int k=0; k<dirichletDofs.size (); ++k)
    {
        int pos = dirichletDofs(k);
        destination.row(pos).array()=dirichletValues.row(k);
    }
    return;
}


} // namespace gismo
