/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <gismo.h>
#include <gsSpectra/gsSpectra.h>

#include <random>  // Requires C++ 11

using namespace gismo;

template<class T>
gsSparseMatrix<T> sprand(index_t size, T prob = 0.5)
{
    gsSparseMatrix<T> mat(size, size);
    std::default_random_engine gen;
    gen.seed(0);
    std::uniform_real_distribution<T> distr(0.0, 1.0);
    for (index_t i = 0; i < size; i++)
    {
        for (index_t j = 0; j < size; j++)
        {
            if (distr(gen) < prob)
                mat.insert(i, j) = distr(gen) - 0.5;
        }
    }
    return mat;
}

template<class T>
void gen_sparse_data(index_t n, gsSparseMatrix<T> & A, gsSparseMatrix<T> & B, T prob = 0.1)
{
    // Eigen solver only uses the lower triangle of A,
    // so we don't need to make A symmetric here
    A = sprand(n, prob);
    B = A.transpose() * A;
    // To make sure B is positive definite
    for (index_t i = 0; i < n; i++)
        B.coeffRef(i, i) += 0.1;
}

template<class T>
bool check(const gsSparseMatrix<T> & A, const gsVector<T> & evals, const gsMatrix<T> & evecs, T tol = 1e-9)
{
    if (evecs.size()==0)
        return false;
    else
    {
        gsMatrix<T> resid = A.template selfadjointView<gsEigen::Lower>() * evecs - evecs * evals.asDiagonal();
        const T err = resid.array().abs().maxCoeff();
        return math::lessthan(err,tol);
    }
}

template<class T>
bool check(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, const gsVector<T> & evals, const gsMatrix<T> & evecs, T tol = 1e-9)
{
    if (evecs.size()==0)
        return false;
    else
    {
        gsMatrix<T> resid = A.template selfadjointView<gsEigen::Lower>() * evecs -
            B.template selfadjointView<gsEigen::Lower>() * evecs * evals.asDiagonal();
        const T err = resid.array().abs().maxCoeff();
        return math::lessthan(err,tol);
    }
}

template<class T>
bool test_Symm(const gsSparseMatrix<T> & A, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraSymSolver<gsSparseMatrix<T>> eigs(A,k,m);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    return check(A,eigvals,eigvecs);
}

template<class T>
bool test_SymmShift(const gsSparseMatrix<T> & A, T sigma, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraSymShiftSolver<gsSparseMatrix<T>> eigs(A,k,m,sigma);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    return check(A,eigvals,eigvecs);
}

template<class T>
bool test_Cholesky(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraGenSymSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Cholesky> eigs(A,B,k,m);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    gsDebugVar(eigvals.transpose());
    return check(A,B,eigvals,eigvecs);
}

template<class T>
bool test_RegularInverse(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraGenSymSolver<gsSparseMatrix<T>,Spectra::GEigsMode::RegularInverse> eigs(A,B,k,m);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    gsDebugVar(eigvals.transpose());
    return check(A,B,eigvals,eigvecs);
}

template<class T>
bool test_ShiftInvert(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, T sigma, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraGenSymShiftSolver<gsSparseMatrix<T>,Spectra::GEigsMode::ShiftInvert> eigs(A,B,k,m,sigma);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    return check(A,B,eigvals,eigvecs);
}

template<class T>
bool test_Buckling(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, T sigma, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraGenSymShiftSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Buckling> eigs(A,B,k,m,sigma);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    return check(B,A,eigvals,eigvecs);
}

template<class T>
bool test_Cayley(const gsSparseMatrix<T> & A, const gsSparseMatrix<T> & B, T sigma, index_t k, index_t m, Spectra::SortRule sortRule)
{
    gsSpectraGenSymShiftSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Cayley> eigs(A,B,k,m,sigma);
    eigs.compute(sortRule,100);
    if (eigs.info() != Spectra::CompInfo::Successful)
        gsWarn<<"Solver failed\n";

    gsVector<T> eigvals = eigs.eigenvalues();
    gsMatrix<T> eigvecs = eigs.eigenvectors();
    return check(A,B,eigvals,eigvecs);
}

int main(int argc, char *argv[])
{
    gsSparseMatrix<> A, B;
    index_t k, m;
    real_t sigma;

/*
    gen_sparse_data(10, A, B, 0.5);
    k = 3;
    m = 6;
    sigma = 1.0;
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::BothEnds));

    gen_sparse_data(100, A, B, 0.1);
    k = 10;
    m = 20;
    sigma = 10.0;
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::BothEnds));

    gen_sparse_data(1000, A, B, 0.01);
    k = 20;
    m = 50;
    sigma = 100.0;
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Symm(A,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_SymmShift(A,sigma,k,m,Spectra::SortRule::BothEnds));
*/

    gen_sparse_data(10, A, B, 0.5);
    k = 3;
    m = 6;
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::BothEnds));


    gen_sparse_data(100, A, B, 0.1);
    k = 10;
    m = 20;

    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::BothEnds));
    return 1;
/*

    sigma = 1.2345;
    gsDebugVar(test_ShiftInvert(A,B,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_ShiftInvert(A,B,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_ShiftInvert(A,B,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_ShiftInvert(A,B,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_ShiftInvert(A,B,sigma,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_Buckling(A,B,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Buckling(A,B,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Buckling(A,B,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Buckling(A,B,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Buckling(A,B,sigma,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_Cayley(A,B,sigma,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Cayley(A,B,sigma,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Cayley(A,B,sigma,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Cayley(A,B,sigma,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Cayley(A,B,sigma,k,m,Spectra::SortRule::BothEnds));

    gen_sparse_data(1000, A, B, 0.01);
    k = 20;
    m = 50;
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_Cholesky(A,B,k,m,Spectra::SortRule::BothEnds));

    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::LargestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestMagn));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::SmallestAlge));
    gsDebugVar(test_RegularInverse(A,B,k,m,Spectra::SortRule::BothEnds));
*/
    return 1;
}
