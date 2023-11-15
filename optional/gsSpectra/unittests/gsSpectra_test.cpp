/** @file gsSpectra_test.cpp

    @brief Provides unittests for the gsSpectra class

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): H.M.Verhelst (2019 - ..., TU Delft, 2023 - ... UniFi)
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsSpectra/gsSpectra.h>
#include <random>  // Requires C++ 11

SUITE(gsSpectra_test)                 // The suite should have the same name as the file
{

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SparseSymmSolver_CHECK(const gsSparseMatrix<real_t> & A,
                                index_t nvalues,
                                index_t nfac,
                                index_t niter,
                                Spectra::SortRule sortRule,
                                bool allowFail = false);

    void SparseSymmShiftSolver_CHECK(const gsSparseMatrix<real_t> & A,
                                     real_t sigma,
                                     index_t nvalues,
                                     index_t nfac,
                                     index_t niter,
                                     Spectra::SortRule sortRule,
                                     bool allowFail = false);

    void SparseCholesky_CHECK(const gsSparseMatrix<real_t> & A,
                              const gsSparseMatrix<real_t> & B,
                              index_t nvalues,
                              index_t nfac,
                              index_t niter,
                              Spectra::SortRule sortRule,
                              bool allowFail = false);

    void SparseRegularInverse_CHECK(const gsSparseMatrix<real_t> & A,
                                    const gsSparseMatrix<real_t> & B,
                                    index_t nvalues,
                                    index_t nfac,
                                    index_t niter,
                                    Spectra::SortRule sortRule,
                                    bool allowFail = false);

    void SparseShiftInvert_CHECK(const gsSparseMatrix<real_t> & A,
                                 const gsSparseMatrix<real_t> & B,
                                 real_t sigma,
                                 index_t nvalues,
                                 index_t nfac,
                                 index_t niter,
                                 Spectra::SortRule sortRule,
                                 bool allowFail = false);

    void SparseBuckling_CHECK(const gsSparseMatrix<real_t> & A,
                              const gsSparseMatrix<real_t> & B,
                              real_t sigma,
                              index_t nvalues,
                              index_t nfac,
                              index_t niter,
                              Spectra::SortRule sortRule,
                              bool allowFail = false);

    void SparseCayley_CHECK(const gsSparseMatrix<real_t> & A,
                            const gsSparseMatrix<real_t> & B,
                            real_t sigma,
                            index_t nvalues,
                            index_t nfac,
                            index_t niter,
                            Spectra::SortRule sortRule,
                            bool allowFail = false);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Auxiliary functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsSparseMatrix<real_t> sprand(index_t size, real_t prob = 0.5)
    {
        gsSparseMatrix<real_t> mat(size, size);
        std::default_random_engine gen;
        gen.seed(0);
        std::uniform_real_distribution<real_t> distr(0.0, 1.0);
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

    void gen_sparse_data(index_t n, gsSparseMatrix<real_t> & A, gsSparseMatrix<real_t> & B, real_t prob = 0.1)
    {
        // Eigen solver only uses the lower triangle of A,
        // so we don't need to make A symmetric here
        A = sprand(n, prob);
        B = A.transpose() * A;
        // To make sure B is positive definite
        for (index_t i = 0; i < n; i++)
            B.coeffRef(i, i) += 0.1;
    }

    bool check(const gsSparseMatrix<real_t> & A, const gsVector<real_t> & evals, const gsMatrix<real_t> & evecs, real_t tol = 1e-9)
    {
        if (evecs.size()==0)
            return false;
        else
        {
            gsMatrix<real_t> resid = A.template selfadjointView<gsEigen::Lower>() * evecs - evecs * evals.asDiagonal();
            const real_t err = resid.array().abs().maxCoeff();
            return math::lessthan(err,tol);
        }
    }

    bool check(const gsSparseMatrix<real_t> & A, const gsSparseMatrix<real_t> & B, const gsVector<real_t> & evals, const gsMatrix<real_t> & evecs, real_t tol = 1e-9)
    {
        if (evecs.size()==0)
            return false;
        else
        {
            gsMatrix<real_t> resid = A.template selfadjointView<gsEigen::Lower>() * evecs -
                B.template selfadjointView<gsEigen::Lower>() * evecs * evals.asDiagonal();
            const real_t err = resid.array().abs().maxCoeff();
            return math::lessthan(err,tol);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseSymm_10)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(10,A,B,0.5);
        index_t k = 3;
        index_t m = 6;
        index_t niter = 1000;
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseSymm_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 1000;
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseSymm_1000)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(1000,A,B,0.01);
        index_t k = 20;
        index_t m = 50;
        index_t niter = 1000;
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestMagn);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmSolver_CHECK(A,k,m,niter,Spectra::SortRule::BothEnds);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseSymmShift_10)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(10,A,B,0.5);
        index_t k = 3;
        index_t m = 6;
        index_t niter = 100;
        real_t sigma = 1.0;
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseSymmShift_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 100;
        real_t sigma = 10.0;
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseSymmShift_1000)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(1000,A,B,0.01);
        index_t k = 20;
        index_t m = 50;
        index_t niter = 100;
        real_t sigma = 100.0;
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseSymmShiftSolver_CHECK(A,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseCholesky_10)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(10,A,B,0.5);
        index_t k = 3;
        index_t m = 6;
        index_t niter = 300;
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseCholesky_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 300;
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseCholesky_1000)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(1000,A,B,0.01);
        index_t k = 20;
        index_t m = 50;
        index_t niter = 300;
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseCholesky_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseRegularInverse_10)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(10,A,B,0.5);
        index_t k = 3;
        index_t m = 6;
        index_t niter = 100;
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseRegularInverse_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 100;
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }

    TEST(SparseRegularInverse_1000)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(1000,A,B,0.01);
        index_t k = 20;
        index_t m = 50;
        index_t niter = 100;
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseRegularInverse_CHECK(A,B,k,m,niter,Spectra::SortRule::BothEnds);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseShiftInvert_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 100;
        real_t sigma = 1.2345;
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseBuckling_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 100;
        real_t sigma = 1.2345;
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TEST(SparseCayley_100)
    {
        gsSparseMatrix<real_t> A,B;
        gen_sparse_data(100,A,B,0.1);
        index_t k = 10;
        index_t m = 20;
        index_t niter = 100;
        real_t sigma = 1.2345;
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestMagn);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::LargestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestMagn,true);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::SmallestAlge);
        SparseShiftInvert_CHECK(A,B,sigma,k,m,niter,Spectra::SortRule::BothEnds);
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SparseSymmSolver_CHECK(const gsSparseMatrix<real_t> & A,
                                index_t nvalues,
                                index_t nfac,
                                index_t niter,
                                Spectra::SortRule sortRule,
                                bool allowFail)
    {
        gsSpectraSymSolver<gsSparseMatrix<real_t>> eigs(A,nvalues,nfac);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");
        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,eigvals,eigvecs) || allowFail,true);
    }

    void SparseSymmShiftSolver_CHECK(const gsSparseMatrix<real_t> & A,
                                     real_t sigma,
                                     index_t nvalues,
                                     index_t nfac,
                                     index_t niter,
                                     Spectra::SortRule sortRule,
                                     bool allowFail)
    {
        gsSpectraSymShiftSolver<gsSparseMatrix<real_t>> eigs(A,nvalues,nfac,sigma);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,eigvals,eigvecs) || allowFail,true);
    }

    void SparseCholesky_CHECK(const gsSparseMatrix<real_t> & A,
                              const gsSparseMatrix<real_t> & B,
                              index_t nvalues,
                              index_t nfac,
                              index_t niter,
                              Spectra::SortRule sortRule,
                              bool allowFail)
    {
        gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> eigs(A,B,nvalues,nfac);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,B,eigvals,eigvecs) || allowFail,true);
    }

    void SparseRegularInverse_CHECK(const gsSparseMatrix<real_t> & A,
                                    const gsSparseMatrix<real_t> & B,
                                    index_t nvalues,
                                    index_t nfac,
                                    index_t niter,
                                    Spectra::SortRule sortRule,
                                    bool allowFail)
    {
        gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::RegularInverse> eigs(A,B,nvalues,nfac);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,B,eigvals,eigvecs) || allowFail,true);
    }

    void SparseShiftInvert_CHECK(const gsSparseMatrix<real_t> & A,
                                 const gsSparseMatrix<real_t> & B,
                                 real_t sigma,
                                 index_t nvalues,
                                 index_t nfac,
                                 index_t niter,
                                 Spectra::SortRule sortRule,
                                 bool allowFail)
    {
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> eigs(A,B,nvalues,nfac,sigma);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,B,eigvals,eigvecs) || allowFail,true);
    }

    void SparseBuckling_CHECK(const gsSparseMatrix<real_t> & A,
                              const gsSparseMatrix<real_t> & B,
                              real_t sigma,
                              index_t nvalues,
                              index_t nfac,
                              index_t niter,
                              Spectra::SortRule sortRule,
                              bool allowFail)
    {
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Buckling> eigs(A,B,nvalues,nfac,sigma);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,B,eigvals,eigvecs) || allowFail,true);
    }

    void SparseCayley_CHECK(const gsSparseMatrix<real_t> & A,
                            const gsSparseMatrix<real_t> & B,
                            real_t sigma,
                            index_t nvalues,
                            index_t nfac,
                            index_t niter,
                            Spectra::SortRule sortRule,
                            bool allowFail)
    {
        gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cayley> eigs(A,B,nvalues,nfac,sigma);
        eigs.compute(sortRule,niter);
        GISMO_ENSURE(eigs.info() == Spectra::CompInfo::Successful || allowFail,"Solver failed, but it is not supposed to do so");        gsVector<real_t> eigvals = eigs.eigenvalues();
        gsMatrix<real_t> eigvecs = eigs.eigenvectors();
        CHECK_EQUAL(check(A,B,eigvals,eigvecs) || allowFail,true);
    }
}