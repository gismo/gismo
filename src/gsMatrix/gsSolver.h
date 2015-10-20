/** @file gsSolver.h

    @brief abstract interfaces for solvers and wrapper around Eigen solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

namespace gismo {

// forward declarations
template<typename T> class gsEigenCGIdentity;
template<typename T> class gsEigenCGDiagonal;
template<typename T> class gsEigenBiCGSTABILUT;
template<typename T> class gsEigenLU;
template<typename T> class gsEigenQR;
template<typename T> class gsEigenSimplicialLDLT;

/** @brief Abstract class for solvers.
 The solver interface is base on 3 methods:
  -compute set the system matrix (possibly compute the factorization or preconditioners)
  -solve solves for a given right hand side
  -succeed returns true if solving succeded according to solver dependent criteria
   (usually tolerance based)
 So in order to solve \f$ A x = b \f$ with a solver \a s two functions must be called:
 s.compute(A) and s.solve(b). The calls can be chained as in  s.compute(A).solve(b).
 \ingroup Solver
*/
template <typename T>
class gsSolver
{
public:
    typedef gsEigenCGIdentity<T> CGIdentity;
    typedef gsEigenCGDiagonal<T> CGDiagonal;
    typedef gsEigenBiCGSTABILUT<T> BiCGSTABILUT;
    typedef gsEigenLU<T> LU;
    typedef gsEigenQR<T> RQ;
    typedef gsEigenSimplicialLDLT<T> SimplicialLDLT;

public:
    typedef gsSparseMatrix<T> MatrixT;
    typedef gsMatrix<T>       VectorT;

public:
    virtual ~gsSolver(){}

    virtual gsSolver& compute (const MatrixT &matrix)       = 0;

    virtual VectorT   solve   (const VectorT &rhs)    const = 0;

    virtual bool      succeed ()                      const = 0;
};


#define GISMO_EIGEN_SPARSE_SOLVER(gsname, eigenName) \
template<typename T=real_t> \
class gsname : public gsSolver<T>, public gsEigenAdaptor<T>::eigenName \
{\
    typedef typename gsSolver<T>::MatrixT MatrixT;\
    typedef typename gsSolver<T>::VectorT VectorT;\
protected:\
    index_t m_rows;\
    index_t m_cols;\
public: \
    gsname() \
        : m_rows(0),m_cols(0)\
    {} \
    gsname(const MatrixT &matrix) \
        : gsEigenAdaptor<T>::eigenName(matrix), m_rows(matrix.rows()),m_cols(matrix.cols()) \
    {} \
    gsname& compute   (const MatrixT &matrix) \
    { \
        m_rows=matrix.rows();\
        m_cols=matrix.cols();\
        gsEigenAdaptor<T>::eigenName::compute(matrix); \
        return *this; \
    } \
    VectorT solve  (const VectorT &rhs) const\
    { \
        return gsEigenAdaptor<T>::eigenName::solve(rhs); \
    } \
    bool succeed() const\
    { \
        return gsEigenAdaptor<T>::eigenName::info()==Eigen::Success; \
    } \
    index_t rows() const {return m_rows;}\
    index_t cols() const {return m_cols;}\
};

GISMO_EIGEN_SPARSE_SOLVER (gsEigenCGIdentity,     CGIdentity)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenCGDiagonal,     CGDiagonal)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenBiCGSTABILUT,   BiCGSTABILUT)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenLU,             LU)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenQR,             QR)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenSimplicialLDLT, SimplicialLDLT)

#undef GISMO_EIGEN_SPARSE_SOLVER

}
