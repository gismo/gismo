/** @file gsSpectra.h

    @brief Header file for using Spectra extension

    https://spectralib.org/doc

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>

#include <Spectra/Spectra/SymEigsSolver.h>
#include <Spectra/Spectra/SymEigsShiftSolver.h>
#include <Spectra/Spectra/SymGEigsSolver.h>
#include <Spectra/Spectra/GenEigsSolver.h>
#include <Spectra/Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/Spectra/MatOp/SparseGenMatProd.h>
//#include <Spectra/Spectra/MatOp/DenseSymMatProd.h> // included by SymEigsSolver.h
#include <Spectra/Spectra/MatOp/SparseCholesky.h>
#include <Spectra/Spectra/MatOp/DenseCholesky.h>

namespace gismo {

//using Spectra::SELECT_EIGENVALUE; // forward

// Product operation wrapper
template <class MatrixType> class SpectraMatProd
{
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Nested NestedMatrix;
    NestedMatrix m_mat;
public:
    SpectraMatProd(const MatrixType& mat) : m_mat(mat) { }
    int rows() const { return m_mat.rows(); }
    int cols() const { return m_mat.cols(); }
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        gsAsVector<Scalar>(y_out, m_mat.rows()).noalias() =
            m_mat * gsAsConstVector<Scalar>(x_in,  m_mat.cols());
    }
};

/** \brief Eigenvalue solver for general real matrices

    Typical usage:
    \code
    gsSparseMatrix<> A;
    const int nEv = 2;
    gsSpectraSymSolver<gsSparseMatrix<> > slv(A, nEv, 2*nEv);
    slv.compute();
    if ( slv.info() == 0 )
    {
    gsInfo << slv.eigenvalues()  <<"\n";
    gsInfo << slv.eigenvectors() <<"\n";
    }
    \code
*/
template <class MatrixType, int SelRule = Spectra::SMALLEST_MAGN>
class gsSpectraSolver : private SpectraMatProd<MatrixType>,
        public Spectra::GenEigsSolver<typename MatrixType::Scalar,
                                      SelRule, SpectraMatProd<MatrixType> >
{
    typedef SpectraMatProd<MatrixType> MatOp;
    typedef Spectra::GenEigsSolver<typename MatrixType::Scalar,
                                   SelRule, MatOp > Base;
public:
    gsSpectraSolver(const MatrixType & mat, int nev_, int ncv_) :
    MatOp(mat), Base(this, nev_, ncv_) { Base::init(); }
};

/// Eigenvalue solver for real symmetric matrices
template <class MatrixType, int SelRule = Spectra::SMALLEST_ALGE>
class gsSpectraSymSolver : private SpectraMatProd<MatrixType>,
        public Spectra::SymEigsSolver<typename MatrixType::Scalar,
                                      SelRule, SpectraMatProd<MatrixType> >
{
    typedef SpectraMatProd<MatrixType> MatOp;
    typedef Spectra::SymEigsSolver<typename MatrixType::Scalar,
                                   SelRule, MatOp> Base;
public:
    gsSpectraSymSolver(const MatrixType & mat, int nev_, int ncv_) :
    MatOp(mat), Base(this, nev_, ncv_) { Base::init(); }
};


template <class MatrixType> class SpectraOps
{
public:
    typedef Spectra::SparseCholesky<typename MatrixType::Scalar> InvOp;
    SpectraOps(const MatrixType & A, const MatrixType & B) : opA(A), opB(B) { }
    SpectraMatProd<MatrixType>                           opA;
    Spectra::SparseCholesky<typename MatrixType::Scalar> opB;
};

//template<> //compilation fails with this
template <class T> class SpectraOps<gsMatrix<T> >
{
public:
    typedef Spectra::DenseCholesky<T> InvOp;
    typedef gsMatrix<T> MatrixType;
protected:
    SpectraOps(const MatrixType & A, const MatrixType & B) : opA(A), opB(B) { }
    SpectraMatProd<MatrixType>                          opA;
    InvOp opB;
};

/// Generalized eigenvalue solver for real symmetric matrices
template <class MatrixType, int SelRule = Spectra::SMALLEST_ALGE, int GEigsMode = Spectra::GEIGS_CHOLESKY>
class gsSpectraGenSymSolver : private SpectraOps<MatrixType>,
public Spectra::SymGEigsSolver<typename MatrixType::Scalar, SelRule,
SpectraMatProd<MatrixType>, typename SpectraOps<MatrixType>::InvOp, GEigsMode>
{
    typedef typename MatrixType::Scalar Scalar;
    typedef SpectraOps<MatrixType> Ops;
    typedef SpectraMatProd<MatrixType> MatOp;

    typedef Spectra::SymGEigsSolver<Scalar,SelRule, MatOp, typename Ops::InvOp,
                                    GEigsMode> Base;
public:
    gsSpectraGenSymSolver(const MatrixType & Amat, const MatrixType & Bmat, int nev_, int ncv_)
    : Ops(Amat,Bmat), Base(&this->opA, &this->opB, nev_, math::min(ncv_,Amat.rows()))
    { Base::init(); }
};

} //namespace gismo
