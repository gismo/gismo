/** @file gsLinearAlgebra.h

    @brief This is the main header file that collects wrappers of Eigen for linear algebra.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


# pragma once

#include <gsCore/gsMath.h>

// Eigen linear algebra library (http://eigen.tuxfamily.org)

// Plugin provides extra members
#define EIGEN_MATRIXBASE_PLUGIN <gsMatrix/gsMatrixAddons.h>
#include <gsMatrix/gsEigenDeclarations.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#if defined(GISMO_WITH_MPFR)
#include <unsupported/Eigen/MPRealSupport>
#endif

#if defined(GISMO_WITH_MPQ)
#include <unsupported/Eigen/MPGmpqSupport.h>
#endif

// Extra Eigen code
#include <gsMatrix/gsEigenAddons.h>

#ifdef GISMO_WITH_SUPERLU
#include <Eigen/SuperLUSupport>
#endif

#ifdef GISMO_WITH_PARDISO
#include <Eigen/PardisoSupport>
#endif

#ifdef GISMO_WITH_PASTIX
#include <Eigen/PaStiXSupport>
#endif

// sparsesuite
//#include <Eigen/UmfPackSupport>
//#include <Eigen/SPQRSupport>
//#include <Eigen/CholmodSupport>

// METIS
//#include <Eigen/MetisSupport>

// PaStiX
//#include <Eigen/PaStiXSupport>


namespace gismo
{

using Eigen::internal::cast; // from Core/MathFunctions.h

//Global variables related to gsMatrix ( see also
//gsForwardDeclarations,h)
using Eigen::Dynamic ;//=-1

using Eigen::Lower;
using Eigen::Upper;

// Values for matrix align options
using Eigen::RowMajor;//=0
using Eigen::ColMajor;//=1
using Eigen::AutoAlign;//

template<class T, int _Rows, int _Cols> class gsAsMatrix ;
template<class T, int _Rows, int _Cols> class gsAsConstMatrix ;

template<class T, int _Rows> class gsAsVector ;
template<class T, int _Rows> class gsAsConstVector ;

// helper template for changing the dimension of a matrix
template <int Dim, int Change>
struct ChangeDim
{
    enum { D = Dim + Change };
};
template <int Change>
struct ChangeDim<Dynamic, Change>
{
    enum { D = Dynamic };
};


/**
   @brief Adaptor for Eigen types
*/
template<typename T>
struct gsEigenAdaptor
{
public:
    // Note: IncompleteILU is not compatible with
    // Eigen::ConjugateGradient because this preconditionner does not
    // preserve symmetry.

    /// Congugate gradient without preconditioner (identity as preconditioner) 
    /// The matrix is assumed symmetric and only the lower trinagular part is used
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<T,0,index_t>,
            Eigen::Lower, Eigen::IdentityPreconditioner> CGIdentity;

    /// Congugate gradient with diagonal (Jacobi) preconditioner
    /// The matrix is assumed symmetric and only the lower trinagular part is used
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<T,0,index_t>, 
            Eigen::Lower, Eigen::DiagonalPreconditioner<T> > CGDiagonal;

    /// BiCGSTAB with Incomplete LU factorization with dual-threshold strategy
    typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<T,0,index_t>,
                            Eigen::IncompleteLUT<T> > BiCGSTABILUT;

    /// BiCGSTAB with Diagonal (Jacobi) preconditioner
    typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<T,0,index_t>,
                            Eigen::DiagonalPreconditioner<T> > BiCGSTABDiagonal;

    /// BiCGSTAB without preconditioner (identity as preconditioner) 
    typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<T,0,index_t>,
                            Eigen::IdentityPreconditioner > BiCGSTABIdentity;

    /// Direct LDLt factorization
    typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<T,0,index_t> > SimplicialLDLT;

    /// Sparse LU solver
    typedef Eigen::SparseLU<Eigen::SparseMatrix<T,0,index_t>,
                            Eigen::COLAMDOrdering<index_t> > SparseLU;

    /// Sparse QR solver
    typedef Eigen::SparseQR<Eigen::SparseMatrix<T,0,index_t>,
                            Eigen::COLAMDOrdering<index_t> > SparseQR;

    #ifdef GISMO_WITH_SUPERLU
    /// SuperLU (if enabled)
    typedef Eigen::SuperLU<Eigen::SparseMatrix<T,0,index_t> > SuperLU;
    #endif

    #ifdef GISMO_WITH_PARDISO
    /// Pardiso (if enabled)
    typedef Eigen::PardisoLDLT<Eigen::SparseMatrix<T,0,index_t> > PardisoLDLT;
    typedef Eigen::PardisoLLT <Eigen::SparseMatrix<T,0,index_t> > PardisoLLT;
    typedef Eigen::PardisoLU  <Eigen::SparseMatrix<T,0,index_t> > PardisoLU;
    #endif

};

} // namespace gismo


#include <gsMatrix/gsMatrixBlockView.h>
#include <gsMatrix/gsMatrix.h>
#include <gsMatrix/gsVector.h>
#include <gsMatrix/gsAsMatrix.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseVector.h>
#include <gsMatrix/gsSparseSolver.h>
