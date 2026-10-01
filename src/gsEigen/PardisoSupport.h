/*
 Copyright (c) 2011, Intel Corporation. All rights reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Intel Corporation nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 ********************************************************************************
 *   Content : Eigen bindings to Intel(R) MKL PARDISO
 ********************************************************************************
*/

/* Modified copy of Eigen 5.0.0 Eigen/src/PardisoSupport/PardisoSupport.h for G+Smo: binds the reference
   (non-MKL) PARDISO library through extern "C" declarations, passes 1-based CSR indices, and renames the classes
   to gsPardiso{Impl,LU,LLT,LDLT} in namespace gsEigen, so they coexist with the MKL-based PardisoLU/LLT/LDLT
   from stock Eigen's PardisoSupport module. */

#include <gsCore/gsLinearAlgebra.h>

#ifndef GISMO_PARDISOSUPPORT_H
#define GISMO_PARDISOSUPPORT_H

#include <cstring>

#ifdef EIGEN_USE_MKL
#error "gsEigen/PardisoSupport.h is the reference-PARDISO backend; with MKL use <Eigen/PardisoSupport>"
#endif

extern "C"
{
void pardiso( void *, int *, int *, int *, int *, int * , void *, int *,
              int * , int *, int *, int *, int *, void *, void *, int *);
void pardiso_chkmatrix (int *, int *, void *, int *, int *, int  *);
void pardiso_chkvec    (int *, int *, void *, int  *);
void pardiso_printstats(int *, int *, void *, int *, int *, int *, void *, int  *);
} // extern "C"

namespace gsEigen {

template <typename MatrixType_>
class gsPardisoLU;
template <typename MatrixType_, int Options = Upper>
class gsPardisoLLT;
template <typename MatrixType_, int Options = Upper>
class gsPardisoLDLT;

namespace internal {
template <typename IndexType>
struct gs_pardiso_run_selector {
  static IndexType run(void* pt, IndexType maxfct, IndexType mnum, IndexType type, IndexType phase,
                       IndexType n, void* a, IndexType* ia, IndexType* ja, IndexType* perm, IndexType nrhs,
                       IndexType* iparm, IndexType msglvl, void* b, void* x) {
    IndexType error = 0;
    ::pardiso(pt, &maxfct, &mnum, &type, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
    return error;
  }
};

template <class Pardiso>
struct gs_pardiso_traits;

template <typename MatrixType_>
struct gs_pardiso_traits<gsPardisoLU<MatrixType_> > {
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType_::Scalar Scalar;
  typedef typename MatrixType_::RealScalar RealScalar;
  typedef typename MatrixType_::StorageIndex StorageIndex;
};

template <typename MatrixType_, int Options>
struct gs_pardiso_traits<gsPardisoLLT<MatrixType_, Options> > {
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType_::Scalar Scalar;
  typedef typename MatrixType_::RealScalar RealScalar;
  typedef typename MatrixType_::StorageIndex StorageIndex;
};

template <typename MatrixType_, int Options>
struct gs_pardiso_traits<gsPardisoLDLT<MatrixType_, Options> > {
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType_::Scalar Scalar;
  typedef typename MatrixType_::RealScalar RealScalar;
  typedef typename MatrixType_::StorageIndex StorageIndex;
};

}  // end namespace internal

template <class Derived>
class gsPardisoImpl : public SparseSolverBase<Derived> {
 protected:
  typedef SparseSolverBase<Derived> Base;
  using Base::derived;
  using Base::m_isInitialized;

  typedef internal::gs_pardiso_traits<Derived> Traits;

 public:
  using Base::_solve_impl;

  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::RealScalar RealScalar;
  typedef typename Traits::StorageIndex StorageIndex;
  typedef SparseMatrix<Scalar, RowMajor, StorageIndex> SparseMatrixType;
  typedef Matrix<Scalar, Dynamic, 1> VectorType;
  typedef Matrix<StorageIndex, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
  typedef Matrix<StorageIndex, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
  typedef Array<StorageIndex, 64, 1, DontAlign> ParameterType;
  enum { ScalarIsComplex = NumTraits<Scalar>::IsComplex, ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };

  gsPardisoImpl() : m_analysisIsOk(false), m_factorizationIsOk(false) {
    m_iparm.setZero();
    m_msglvl = 0;  // No output
    m_isInitialized = false;
  }

  ~gsPardisoImpl() { pardisoRelease(); }

  inline Index cols() const { return m_size; }
  inline Index rows() const { return m_size; }

  /** \brief Reports whether previous computation was successful.
   *
   * \returns \c Success if computation was successful,
   *          \c NumericalIssue if the matrix appears to be negative.
   */
  ComputationInfo info() const {
    GISMO_ASSERT(m_isInitialized, "Decomposition is not initialized.");
    return m_info;
  }

  /** \warning for advanced usage only.
   * \returns a reference to the parameter array controlling PARDISO.
   * See the PARDISO manual to know how to use it. */
  ParameterType& pardisoParameterArray() { return m_iparm; }

  /// Sets the PARDISO control parameter \a i (0-based index into iparm) to \a value.
  void setParam(const int i, const int value) { m_iparm[i] = value; }

  /** Performs a symbolic decomposition on the sparsity of \a matrix.
   *
   * This function is particularly useful when solving for several problems having the same structure.
   *
   * \sa factorize()
   */
  Derived& analyzePattern(const MatrixType& matrix);

  /** Performs a numeric decomposition of \a matrix
   *
   * The given matrix must has the same sparsity than the matrix on which the symbolic decomposition has been performed.
   *
   * \sa analyzePattern()
   */
  Derived& factorize(const MatrixType& matrix);

  Derived& compute(const MatrixType& matrix);

  template <typename Rhs, typename Dest>
  void _solve_impl(const MatrixBase<Rhs>& b, MatrixBase<Dest>& dest) const;

 protected:
  void pardisoRelease() {
    if (m_isInitialized)  // Factorization ran at least once
    {
      internal::gs_pardiso_run_selector<StorageIndex>::run(m_pt, 1, 1, m_type, -1,
                                                        internal::convert_index<StorageIndex>(m_size), 0, 0, 0,
                                                        m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);
      m_isInitialized = false;
    }
  }

  void pardisoInit(int type) {
    m_type = type;
    bool symmetric = std::abs(m_type) < 10;
    m_iparm[0] = 1;                   // No solver default
    m_iparm[1] =                      // 2: METIS nested-dissection ordering, 3: its OpenMP-parallel variant
#ifdef GISMO_WITH_OPENMP
      3;
#else
      2;
#endif
    m_iparm[2] = 0;                   // Reserved. Set to zero. (??Numbers of processors, value of OMP_NUM_THREADS??)
    m_iparm[3] = 0;                   // No iterative-direct algorithm
    m_iparm[4] = 0;                   // No user fill-in reducing permutation
    m_iparm[5] = 0;                   // Write solution into x, b is left unchanged
    m_iparm[6] = 0;                   // Not in use
    m_iparm[7] = 2;                   // Max numbers of iterative refinement steps
    m_iparm[8] = 0;                   // Not in use
    m_iparm[9] = 13;                  // Perturb the pivot elements with 1E-13
    m_iparm[10] = symmetric ? 0 : 1;  // Use nonsymmetric permutation and scaling MPS
    m_iparm[11] = 0;                  // Not in use
    m_iparm[12] = symmetric ? 0 : 1;  // Maximum weighted matching algorithm is switched-off (default for symmetric).
                                      // Try m_iparm[12] = 1 in case of inappropriate accuracy
    m_iparm[13] = 0;                  // Output: Number of perturbed pivots
    m_iparm[14] = 0;                  // Not in use
    m_iparm[15] = 0;                  // Not in use
    m_iparm[16] = 0;                  // Not in use
    m_iparm[17] = 0;                  // Output: Number of nonzeros in the factor LU
    m_iparm[18] = 0;                  // Output: Mflops for LU factorization
    m_iparm[19] = 0;                  // Output: Numbers of CG Iterations

    m_iparm[20] = 0;  // 1x1 pivoting
    m_iparm[26] = 0;  // No matrix checker
    m_iparm[27] = (sizeof(RealScalar) == 4) ? 1 : 0;
    m_iparm[34] = 0;  // iparm[34] (C-style 0-based indexing) is MKL-only; reference PARDISO expects 1-based
                       // (Fortran) CSR, see getMatrix
    m_iparm[36] = 0;  // CSR
    m_iparm[59] = 0;  // 0 - In-Core ; 1 - Automatic switch between In-Core and Out-of-Core modes ; 2 - Out-of-Core

    memset(m_pt, 0, sizeof(m_pt));
  }

 protected:
  // cached data to reduce reallocation, etc.

  void manageErrorCode(Index error) const {
    switch (error) {
      case 0:
        m_info = Success;
        break;
      case -4:
      case -7:
        m_info = NumericalIssue;
        break;
      default:
        m_info = InvalidInput;
    }
  }

  /// Shifts the compressed CSR arrays of m_matrix from 0-based to the 1-based (Fortran) indices reference PARDISO
  /// requires. O(nnz + n). m_matrix is thereafter only valid as PARDISO input; getMatrix() rebuilds it every call.
  void toOneBased()
  {
    const Index nnz   = m_matrix.nonZeros();   // read before the outer index array is modified
    const Index outer = m_matrix.outerSize();
    StorageIndex* ja = m_matrix.innerIndexPtr();
    StorageIndex* ia = m_matrix.outerIndexPtr();
    for (Index k = 0; k < nnz;    ++k) ++ja[k];
    for (Index k = 0; k <= outer; ++k) ++ia[k];
  }

  mutable SparseMatrixType m_matrix;
  mutable ComputationInfo m_info;
  bool m_analysisIsOk, m_factorizationIsOk;
  StorageIndex m_type, m_msglvl;
  mutable void* m_pt[64];
  mutable ParameterType m_iparm;
  mutable IntColVectorType m_perm;
  Index m_size;
};

template <class Derived>
Derived& gsPardisoImpl<Derived>::compute(const MatrixType& a) {
  m_size = a.rows();
  GISMO_ASSERT(a.rows() == a.cols(), "");

  pardisoRelease();
  m_perm.setZero(m_size);
  derived().getMatrix(a);

  Index error;
  error = internal::gs_pardiso_run_selector<StorageIndex>::run(
      m_pt, 1, 1, m_type, 12, internal::convert_index<StorageIndex>(m_size), m_matrix.valuePtr(),
      m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(), m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);
  manageErrorCode(error);
  m_analysisIsOk = m_info == Success;
  m_factorizationIsOk = m_info == Success;
  m_isInitialized = true;
  return derived();
}

template <class Derived>
Derived& gsPardisoImpl<Derived>::analyzePattern(const MatrixType& a) {
  m_size = a.rows();
  GISMO_ASSERT(m_size == a.cols(), "");

  pardisoRelease();
  m_perm.setZero(m_size);
  derived().getMatrix(a);

  Index error;
  error = internal::gs_pardiso_run_selector<StorageIndex>::run(
      m_pt, 1, 1, m_type, 11, internal::convert_index<StorageIndex>(m_size), m_matrix.valuePtr(),
      m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(), m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);

  manageErrorCode(error);
  m_analysisIsOk = m_info == Success;
  m_factorizationIsOk = false;
  m_isInitialized = true;
  return derived();
}

template <class Derived>
Derived& gsPardisoImpl<Derived>::factorize(const MatrixType& a) {
  GISMO_ASSERT(m_analysisIsOk, "You must first call analyzePattern()");
  GISMO_ASSERT(m_size == a.rows() && m_size == a.cols(), "");

  derived().getMatrix(a);

  Index error;
  error = internal::gs_pardiso_run_selector<StorageIndex>::run(
      m_pt, 1, 1, m_type, 22, internal::convert_index<StorageIndex>(m_size), m_matrix.valuePtr(),
      m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(), m_perm.data(), 0, m_iparm.data(), m_msglvl, NULL, NULL);

  manageErrorCode(error);
  m_factorizationIsOk = m_info == Success;
  return derived();
}

template <class Derived>
template <typename BDerived, typename XDerived>
void gsPardisoImpl<Derived>::_solve_impl(const MatrixBase<BDerived>& b, MatrixBase<XDerived>& x) const {
  if (m_iparm[0] == 0)  // Factorization was not computed
  {
    m_info = InvalidInput;
    return;
  }

  // Index n = m_matrix.rows();
  Index nrhs = Index(b.cols());
  GISMO_ASSERT(m_size == b.rows(), "");
  GISMO_ASSERT(((MatrixBase<BDerived>::Flags & RowMajorBit) == 0 || nrhs == 1),
               "Row-major right hand sides are not supported");
  GISMO_ASSERT(((MatrixBase<XDerived>::Flags & RowMajorBit) == 0 || nrhs == 1),
               "Row-major matrices of unknowns are not supported");
  GISMO_ASSERT(((nrhs == 1) || b.outerStride() == b.rows()), "");

  //  switch (transposed) {
  //    case SvNoTrans    : m_iparm[11] = 0 ; break;
  //    case SvTranspose  : m_iparm[11] = 2 ; break;
  //    case SvAdjoint    : m_iparm[11] = 1 ; break;
  //    default:
  //      //std::cerr << "Eigen: transposition  option \"" << transposed << "\" not supported by the PARDISO backend\n";
  //      m_iparm[11] = 0;
  //  }

  Scalar* rhs_ptr = const_cast<Scalar*>(b.derived().data());
  Matrix<Scalar, Dynamic, Dynamic, ColMajor> tmp;

  // Pardiso cannot solve in-place
  if (rhs_ptr == x.derived().data()) {
    tmp = b;
    rhs_ptr = tmp.data();
  }

  Index error;
  error = internal::gs_pardiso_run_selector<StorageIndex>::run(
      m_pt, 1, 1, m_type, 33, internal::convert_index<StorageIndex>(m_size), m_matrix.valuePtr(),
      m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr(), m_perm.data(), internal::convert_index<StorageIndex>(nrhs),
      m_iparm.data(), m_msglvl, rhs_ptr, x.derived().data());

  manageErrorCode(error);
}

/** \ingroup Matrix
 * \class gsPardisoLU
 * \brief A sparse direct LU factorization and solver based on the PARDISO library
 *
 * This class allows to solve for A.X = B sparse linear problems via a direct LU factorization
 * using the reference (non-MKL) PARDISO library. The sparse matrix A must be squared and invertible.
 * The vectors or matrices X and B can be either dense or sparse.
 *
 * By default, it runs in in-core mode. To enable PARDISO's out-of-core feature, set:
 * \code solver.pardisoParameterArray()[59] = 1; \endcode
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 *
 * \implsparsesolverconcept
 *
 * \sa \ref TutorialSparseSolverConcept, class SparseLU
 */
template <typename MatrixType>
class gsPardisoLU : public gsPardisoImpl<gsPardisoLU<MatrixType> > {
 protected:
  typedef gsPardisoImpl<gsPardisoLU> Base;
  using Base::m_matrix;
  using Base::pardisoInit;
  using Base::toOneBased;
  friend class gsPardisoImpl<gsPardisoLU<MatrixType> >;

 public:
  typedef typename Base::Scalar Scalar;
  typedef typename Base::RealScalar RealScalar;

  using Base::compute;
  using Base::solve;

  gsPardisoLU() : Base() { pardisoInit(Base::ScalarIsComplex ? 13 : 11); }

  explicit gsPardisoLU(const MatrixType& matrix) : Base() {
    pardisoInit(Base::ScalarIsComplex ? 13 : 11);
    compute(matrix);
  }

 protected:
  void getMatrix(const MatrixType& matrix) {
    m_matrix = matrix;
    m_matrix.makeCompressed();
    toOneBased();
  }
};

/** \ingroup Matrix
 * \class gsPardisoLLT
 * \brief A sparse direct Cholesky (LLT) factorization and solver based on the PARDISO library
 *
 * This class allows to solve for A.X = B sparse linear problems via a LL^T Cholesky factorization
 * using the reference (non-MKL) PARDISO library. The sparse matrix A must be selfajoint and positive definite.
 * The vectors or matrices X and B can be either dense or sparse.
 *
 * By default, it runs in in-core mode. To enable PARDISO's out-of-core feature, set:
 * \code solver.pardisoParameterArray()[59] = 1; \endcode
 *
 * \tparam MatrixType the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo can be any bitwise combination of Upper, Lower. The default is Upper, meaning only the upper triangular
 * part has to be used. Upper|Lower can be used to tell both triangular parts can be used as input.
 *
 * \implsparsesolverconcept
 *
 * \sa \ref TutorialSparseSolverConcept, class SimplicialLLT
 */
template <typename MatrixType, int UpLo_>
class gsPardisoLLT : public gsPardisoImpl<gsPardisoLLT<MatrixType, UpLo_> > {
 protected:
  typedef gsPardisoImpl<gsPardisoLLT<MatrixType, UpLo_> > Base;
  using Base::m_matrix;
  using Base::pardisoInit;
  using Base::toOneBased;
  friend class gsPardisoImpl<gsPardisoLLT<MatrixType, UpLo_> >;

 public:
  typedef typename Base::Scalar Scalar;
  typedef typename Base::RealScalar RealScalar;
  typedef typename Base::StorageIndex StorageIndex;
  enum { UpLo = UpLo_ };
  using Base::compute;

  gsPardisoLLT() : Base() { pardisoInit(Base::ScalarIsComplex ? 4 : 2); }

  explicit gsPardisoLLT(const MatrixType& matrix) : Base() {
    pardisoInit(Base::ScalarIsComplex ? 4 : 2);
    compute(matrix);
  }

 protected:
  void getMatrix(const MatrixType& matrix) {
    // PARDISO supports only upper, row-major matrices
    m_matrix.resize(matrix.rows(), matrix.cols());
    m_matrix.template selfadjointView<Upper>() = matrix.template selfadjointView<UpLo>();
    m_matrix.makeCompressed();
    toOneBased();
  }
};

/** \ingroup Matrix
 * \class gsPardisoLDLT
 * \brief A sparse direct Cholesky (LDLT) factorization and solver based on the PARDISO library
 *
 * This class allows to solve for A.X = B sparse linear problems via a LDL^T Cholesky factorization
 * using the reference (non-MKL) PARDISO library. The sparse matrix A is assumed to be selfajoint and positive definite.
 * For complex matrices, A can also be symmetric only, see the \a Options template parameter.
 * The vectors or matrices X and B can be either dense or sparse.
 *
 * By default, it runs in in-core mode. To enable PARDISO's out-of-core feature, set:
 * \code solver.pardisoParameterArray()[59] = 1; \endcode
 *
 * \tparam MatrixType the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam Options can be any bitwise combination of Upper, Lower, and Symmetric. The default is Upper, meaning only the
 * upper triangular part has to be used. Symmetric can be used for symmetric, non-selfadjoint complex matrices, the
 * default being to assume a selfadjoint matrix. Upper|Lower can be used to tell both triangular parts can be used as
 * input.
 *
 * \implsparsesolverconcept
 *
 * \sa \ref TutorialSparseSolverConcept, class SimplicialLDLT
 */
template <typename MatrixType, int Options>
class gsPardisoLDLT : public gsPardisoImpl<gsPardisoLDLT<MatrixType, Options> > {
 protected:
  typedef gsPardisoImpl<gsPardisoLDLT<MatrixType, Options> > Base;
  using Base::m_matrix;
  using Base::pardisoInit;
  using Base::toOneBased;
  friend class gsPardisoImpl<gsPardisoLDLT<MatrixType, Options> >;

 public:
  typedef typename Base::Scalar Scalar;
  typedef typename Base::RealScalar RealScalar;
  typedef typename Base::StorageIndex StorageIndex;
  using Base::compute;
  enum { UpLo = Options & (Upper | Lower) };

  gsPardisoLDLT() : Base() { pardisoInit(Base::ScalarIsComplex ? (bool(Options & Symmetric) ? 6 : -4) : -2); }

  explicit gsPardisoLDLT(const MatrixType& matrix) : Base() {
    pardisoInit(Base::ScalarIsComplex ? (bool(Options & Symmetric) ? 6 : -4) : -2);
    compute(matrix);
  }

  void getMatrix(const MatrixType& matrix) {
    // PARDISO supports only upper, row-major matrices
    m_matrix.resize(matrix.rows(), matrix.cols());
    m_matrix.template selfadjointView<Upper>() = matrix.template selfadjointView<UpLo>();
    m_matrix.makeCompressed();
    toOneBased();
  }
};

}  // end namespace gsEigen

#endif  // GISMO_PARDISOSUPPORT_H
