/** @file gsSpectra.h

    @brief Header file for using Spectra extension

    https://spectralib.org/doc

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris

    Notes:
    -PARDISO notes: https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html

    -Compile on 64bit with MKL machine using (lp64):

    cmake .. -DGISMO_OPTIONAL="gsSpectra" -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx -DGISMO_WITH_PARDISO=ON -DPARDISO_USE_MKL=ON -DEIGEN_USE_MKL_ALL=ON -DGISMO_WITH_OPENMP=ON -DMKL_INTERFACE=lp64 -DTARGET_ARCHITECTURE=none

*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsLinearAlgebra.h>

#define Eigen gsEigen

#include <Spectra/include/Spectra/SymEigsSolver.h>
#include <Spectra/include/Spectra/SymEigsShiftSolver.h>
#include <Spectra/include/Spectra/SymGEigsSolver.h>
#include <Spectra/include/Spectra/SymGEigsShiftSolver.h>
#include <Spectra/include/Spectra/GenEigsSolver.h>
#include <Spectra/include/Spectra/MatOp/SparseSymShiftSolve.h>

namespace gismo {

//using Spectra::SELECT_EIGENVALUE; // forward

// Product operation wrapper
template <class MatrixType> class SpectraMatProd
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
public:
    SpectraMatProd(const MatrixType&&   ) = delete;
    SpectraMatProd(const MatrixType& mat) : m_mat(mat) {}
    index_t rows() const { return m_mat.rows(); }
    index_t cols() const { return m_mat.cols(); }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        GISMO_ASSERT(m_mat.rows()!=0 && m_mat.cols()!=0,"The matrix has zero rows or columns. Is the matrix a temporary (e.g. A-B)?");
        gsAsVector<Scalar>(y_out, m_mat.rows()).noalias() =
            m_mat * gsAsConstVector<Scalar>(x_in,  m_mat.cols());
    }
};

// Product operation wrapper
template <class MatrixType> class SpectraSymMatProd
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
public:
    SpectraSymMatProd(const MatrixType&&   ) = delete;
    SpectraSymMatProd(const MatrixType& mat) : m_mat(mat) {}
    index_t rows() const { return m_mat.rows(); }
    index_t cols() const { return m_mat.cols(); }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        GISMO_ASSERT(m_mat.rows()!=0 && m_mat.cols()!=0,"The matrix has zero rows or columns. Is the matrix a temporary (e.g. A-B)?");
        gsAsVector<Scalar>(y_out, m_mat.rows()).noalias() =
            m_mat.template selfadjointView<gsEigen::Lower>() * gsAsConstVector<Scalar>(x_in,  m_mat.cols());
    }
};

// Shift operation wrapper
template <class MatrixType>
class SpectraMatShiftSolve
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
    const index_t m_n;
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<Scalar>::PardisoLU m_solver;
#elif  GISMO_WITH_SUPERLU
    typename gsSparseSolver<Scalar>::SuperLU m_solver;
#else
    typename gsSparseSolver<Scalar>::LU m_solver;
#endif
public:
    SpectraMatShiftSolve(const MatrixType&&   ) = delete;
    SpectraMatShiftSolve(const MatrixType& mat)
    :
    m_mat(mat), m_n(mat.rows())
    {
        GISMO_ASSERT(m_mat.rows() == m_mat.cols(),"Matrix must be square!");
    }
    index_t rows() const { return m_n; }
    index_t cols() const { return m_n; }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_n).noalias() = m_solver.solve( gsAsConstVector<Scalar>(x_in,  m_n) );
    }

    void set_shift(const Scalar& sigma)
    {
        MatrixType mat(m_mat.template selfadjointView<gsEigen::Lower>());
        MatrixType identity(m_n,m_n);
        identity.setIdentity();

        mat = mat - sigma * identity;
#ifndef GISMO_WITH_PARDISO
#ifndef GISMO_WITH_SUPERLU
        m_solver.isSymmetric(true);
#endif
#endif
        m_solver.compute(mat);

        GISMO_ASSERT(m_solver.info() == gsEigen::Success,"SpectraMatShiftSolve: factorization failed with the given shift");
    }
};

// Shift operation wrapper (dense matrix)
template <class T>
class SpectraMatShiftSolve<gsMatrix<T>>
{
public:
    typedef typename gsMatrix<T>::Nested NestedMatrix;
    typedef typename gsMatrix<T>::Scalar Scalar;
    NestedMatrix m_mat;
    const index_t m_n;
    typename Spectra::BKLDLT<Scalar> m_solver;

public:
    SpectraMatShiftSolve(const gsMatrix<Scalar>&&) = delete;
    SpectraMatShiftSolve(const gsMatrix<Scalar>& mat)
    :
    m_mat(mat), m_n(mat.rows())
    {
        GISMO_ASSERT(m_mat.rows() == m_mat.cols(),"Matrix must be square!");
    }
    index_t rows() const { return m_n; }
    index_t cols() const { return m_n; }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_n).noalias() = m_solver.solve( gsAsConstVector<Scalar>(x_in,  m_n) );
    }

    void set_shift(const Scalar& sigma)
    {
        m_solver.compute(m_mat, gsEigen::Lower, sigma);
        GISMO_ASSERT(m_solver.info() == Spectra::CompInfo::Successful,"SpectraMatShiftSolve: factorization failed with the given shift");
    }
};

// Cholesky wrapper (see Spectra::SparseCholesky)
template <class MatrixType>
class SpectraCholesky
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
    const index_t m_n;
    // NOTE: Does not work for gsSparseSolver<Scalar>::SimplicialLLT!
    gsEigen::SimplicialLLT<gsSparseMatrix<Scalar>, gsEigen::Lower> m_solver;

    Spectra::CompInfo m_info;  // status of the decomposition

public:
    SpectraCholesky(const MatrixType&&   ) = delete;
    SpectraCholesky(const MatrixType& mat)
    :
    m_mat(mat), m_n(mat.rows())
    {
        GISMO_ASSERT(m_mat.rows() == m_mat.cols(),"Matrix A must be square!");

        m_solver.compute(mat);
        m_info = (m_solver.info() == gsEigen::Success) ?
            Spectra::CompInfo::Successful :
            Spectra::CompInfo::NumericalIssue;
    }
    index_t rows() const { return m_n; }
    index_t cols() const { return m_n; }

    /// See Spectra/SparseCholesky.h for help
    Spectra::CompInfo info() const { return m_info; }

    /// See Spectra/SparseCholesky.h for help
    void lower_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsConstVector<Scalar> x(x_in,  m_n);
        gsAsVector<Scalar>      y(y_out, m_n);
        y.noalias() = m_solver.permutationP() * x;
        m_solver.matrixL().solveInPlace(y);
    }

    /// See Spectra/SparseCholesky.h for help
    void upper_triangular_solve(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsConstVector<Scalar> x(x_in,  m_n);
        gsAsVector<Scalar>      y(y_out, m_n);
        y.noalias() = m_solver.matrixU().solve(x);
        y = (m_solver.permutationPinv() * y).eval();
    }
};

// Regular inverse wrapper
template <class MatrixType>
class SpectraRegularInverse
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
    const index_t m_n;
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<Scalar>::PardisoLU m_solver;
#elif  GISMO_WITH_SUPERLU
    typename gsSparseSolver<Scalar>::SuperLU m_solver;
#else
    typename gsSparseSolver<Scalar>::CGDiagonal m_solver;
#endif

    Spectra::CompInfo m_info;  // status of the decomposition

public:
    SpectraRegularInverse(const MatrixType&&   ) = delete;
    SpectraRegularInverse(const MatrixType& mat)
    :
    m_mat(mat), m_n(mat.rows())
    {
        GISMO_ASSERT(m_mat.rows() == m_mat.cols(),"Matrix A must be square!");

        m_solver.compute(mat);
        m_info = (m_solver.info() == gsEigen::Success) ?
            Spectra::CompInfo::Successful :
            Spectra::CompInfo::NumericalIssue;
    }
    index_t rows() const { return m_n; }
    index_t cols() const { return m_n; }

    /// See Spectra/SparseRegularInverse.h for help
    Spectra::CompInfo info() const { return m_info; }

    /// See Spectra/SparseRegularInverse.h for help
    void solve(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_n).noalias() = m_solver.solve( gsAsConstVector<Scalar>(x_in,  m_n) );
    }

    /// See Spectra/SparseRegularInverse.h for help
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_n).noalias() = m_mat.template selfadjointView<gsEigen::Lower>() * gsAsConstVector<Scalar>(x_in,  m_n);
    }
};

// Shift operation wrapper
template <class MatrixType>
class SpectraSymShiftInvert
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_A, m_B;
    const index_t m_n;
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<Scalar>::PardisoLU m_solver;
#elif  GISMO_WITH_SUPERLU
    typename gsSparseSolver<Scalar>::SuperLU m_solver;
#else
    typename gsSparseSolver<Scalar>::LU m_solver;
#endif
public:
    SpectraSymShiftInvert(const MatrixType&&   ,const MatrixType&& ) = delete;
    SpectraSymShiftInvert(const MatrixType& A  ,const MatrixType& B)
    :
    m_A(A), m_B(B), m_n(A.rows())
    {
        GISMO_ASSERT(m_A.rows() == m_A.cols(),"Matrix A must be square!");
        GISMO_ASSERT(m_B.rows() == m_B.cols(),"Matrix B must be square!");
        GISMO_ASSERT(m_B.rows() == m_A.rows(),"Matrix A and B must be the same size!");
    }
    index_t rows() const { return m_n; }
    index_t cols() const { return m_n; }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_n).noalias() = m_solver.solve( gsAsConstVector<Scalar>(x_in,  m_n) );
    }

    void set_shift(const Scalar& sigma)
    {
        MatrixType matA = m_A.template selfadjointView<gsEigen::Lower>();
        MatrixType matB = m_B.template selfadjointView<gsEigen::Lower>();
        MatrixType mat = matA - sigma * matB;

        #ifndef GISMO_WITH_PARDISO
        #ifndef GISMO_WITH_SUPERLU
                m_solver.isSymmetric(true);
        #endif
        #endif
        m_solver.compute(mat);
        GISMO_ASSERT(m_solver.info()==gsEigen::Success,"Factorization failed");
    }
};


/** \brief Eigenvalue solver for general real matrices

    Typical usage:
    \code
    gsSparseMatrix<> A;
    const index_t nEv = 2;
    gsSpectraSymSolver<gsSparseMatrix<> > slv(A, nEv, 2*nEv);
    slv.compute();
    if ( slv.info() == 0 )
    {
    gsInfo << slv.eigenvalues()  <<"\n";
    gsInfo << slv.eigenvectors() <<"\n";
    }
    \endcode
*/
template <class MatrixType>
class gsSpectraSolver : private SpectraMatProd<MatrixType>,
        public Spectra::GenEigsSolver<SpectraMatProd<MatrixType> >
{
    typedef SpectraMatProd<MatrixType> MatOp;
    typedef Spectra::GenEigsSolver<MatOp> Base;
public:
    gsSpectraSolver(const MatrixType &&   , index_t nev_, index_t ncv_) = delete;
    gsSpectraSolver(const MatrixType & mat, index_t nev_, index_t ncv_) :
    MatOp(mat), Base(*this, nev_, ncv_) { Base::init(); }
};

/// Eigenvalue solver for real symmetric matrices
template <class MatrixType>
class gsSpectraSymSolver : private SpectraSymMatProd<MatrixType>,
        public Spectra::SymEigsSolver<SpectraSymMatProd<MatrixType> >
{
    typedef SpectraSymMatProd<MatrixType> MatOp;
    typedef Spectra::SymEigsSolver<MatOp> Base;
public:
    gsSpectraSymSolver(const MatrixType &&   , index_t nev_, index_t ncv_) = delete;
    gsSpectraSymSolver(const MatrixType & mat, index_t nev_, index_t ncv_) :
    MatOp(mat), Base(*this, nev_, ncv_) { Base::init(); }
};

/// Shifted Eigenvalue solver for real symmetric matrices
template <class MatrixType>
class gsSpectraSymShiftSolver :
        private SpectraMatShiftSolve<MatrixType>,
        public Spectra::SymEigsShiftSolver<SpectraMatShiftSolve<MatrixType> >
{
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraMatShiftSolve<MatrixType> Op;
    typedef Spectra::SymEigsShiftSolver<Op> Base;
public:
    gsSpectraSymShiftSolver(const MatrixType &&   , index_t nev_, index_t ncv_, const Scalar& sigma) = delete;
    gsSpectraSymShiftSolver(const MatrixType & mat, index_t nev_, index_t ncv_, const Scalar& sigma) :
    Op(mat), Base(*this, nev_, ncv_,sigma) { Base::init(); }
};

/// SpectraOps is for generalized eigenvalue problems
/// For GEigsMode::Cholesky
template <class MatrixType, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class SpectraOps
{
public:
    typedef SpectraSymMatProd<MatrixType> ProdOp;
    typedef SpectraCholesky<MatrixType>   InvOp;
    SpectraOps(const MatrixType & A, const MatrixType & B) : opA(A), opB(B) { }
    ProdOp opA;
    InvOp  opB;
};

//template<> //compilation fails with this
template <class T> class SpectraOps<gsMatrix<T> >
{
public:
    typedef gsMatrix<T> MatrixType;
    typedef Spectra::DenseCholesky<T>  InvOp;
    typedef SpectraMatProd<MatrixType> ProdOp;
protected:
    SpectraOps(const MatrixType & A, const MatrixType & B) : opA(A), opB(B) { }
    ProdOp opA;
    InvOp  opB;
};

/// For GEigsMode::RegularInverse
/// NOTE: A must be symmetric and B positive semi-definite
template <class MatrixType>
class SpectraOps<MatrixType,Spectra::GEigsMode::RegularInverse>
{
public:
    typedef SpectraRegularInverse<MatrixType> InvOp;
    typedef SpectraSymMatProd<MatrixType>     ProdOp;
    SpectraOps(const MatrixType & A, const MatrixType & B) : opA(A), opB(B) { }
    ProdOp opA;
    InvOp  opB;
};

/// For GEigsMode::ShiftInvert and GEigsMode::Cayley
/// NOTE: A must be symmetric and B positive semi-definite
template <class MatrixType, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::ShiftInvert>
class SpectraShiftOps
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraSymShiftInvert<MatrixType> InvOp;
    typedef SpectraSymMatProd<MatrixType>     ProdOp;

    SpectraShiftOps(const MatrixType & A, const MatrixType & B) : opA(A,B), opB(B) { }
    InvOp  opA;
    ProdOp opB;
};

/// Specialization for GEigsMode::Buckling
/// NOTE: A must be positive semi-definite and B symmetric
template <class MatrixType>
class SpectraShiftOps<MatrixType,Spectra::GEigsMode::Buckling>
{
public:
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraSymShiftInvert<MatrixType> InvOp;
    typedef SpectraSymMatProd<MatrixType>     ProdOp;

    SpectraShiftOps(const MatrixType & A, const MatrixType & B) : opA(A,B), opB(A) { }
    InvOp  opA;
    ProdOp opB;
};

/// GE Solver with shifts. Works for GEigsMode = Cholesky or RegularInverse.
/// See the Spectra Documentation (SymGEigsSolver) for more information
template <class MatrixType, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class gsSpectraGenSymSolver :
    private SpectraOps<MatrixType,GEigsMode>,
    public Spectra::SymGEigsSolver<typename SpectraOps<MatrixType,GEigsMode>::ProdOp, typename SpectraOps<MatrixType,GEigsMode>::InvOp, GEigsMode>
{
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraOps<MatrixType,GEigsMode> OpType;
    typedef typename OpType::ProdOp Op;
    typedef typename OpType::InvOp  BOp;

    typedef Spectra::SymGEigsSolver<Op, BOp,GEigsMode> Base;
public:
    gsSpectraGenSymSolver(const MatrixType &&    , const MatrixType &&    , index_t nev_, index_t ncv_) = delete;
    gsSpectraGenSymSolver(const MatrixType & Amat, const MatrixType & Bmat, index_t nev_, index_t ncv_)
    : OpType(Amat,Bmat), Base(this->opA, this->opB, nev_, math::min(ncv_,Amat.rows()))
    { Base::init(); }
};

/// GE Solver with shifts. Works for GEigsMode = ShiftInvert, Buckling or Cayley
/// See the Spectra Documentation (SymGEigsShiftSolver) for more information
template <class MatrixType, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::ShiftInvert>
class gsSpectraGenSymShiftSolver :
    private SpectraShiftOps<MatrixType,GEigsMode>,
    public Spectra::SymGEigsShiftSolver<typename SpectraShiftOps<MatrixType,GEigsMode>::InvOp, typename SpectraShiftOps<MatrixType,GEigsMode>::ProdOp, GEigsMode>
{
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraShiftOps<MatrixType,GEigsMode> OpType;
    typedef typename OpType::InvOp  Op;
    typedef typename OpType::ProdOp BOp;

    typedef Spectra::SymGEigsShiftSolver<Op, BOp,GEigsMode> Base;
public:
    gsSpectraGenSymShiftSolver(const MatrixType &&    , const MatrixType &&    , index_t nev_, index_t ncv_, const Scalar& sigma) = delete;
    gsSpectraGenSymShiftSolver(const MatrixType & Amat, const MatrixType & Bmat, index_t nev_, index_t ncv_, const Scalar& sigma)
    : OpType(Amat,Bmat), Base(this->opA, this->opB, nev_, math::min(ncv_,Amat.rows()),sigma)
    { Base::init(); }
};

} //namespace gismo


#undef Eigen

