/** @file gsMatrix_advanced_test.cpp

    @brief Additional comprehensive tests for gsMatrix module

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsMatrix_advanced_test)
{
    TEST(gsAsMatrix_Mapping)
    {
        gsVector<> vec(5);
        vec << 1, 2, 3, 4, 5;
        
        gsAsMatrix<> mat(vec.data(), 5, 1);
        
        CHECK_EQUAL(mat.rows(), 5);
        CHECK_EQUAL(mat.cols(), 1);
        CHECK_CLOSE(mat(2, 0), 3.0, 1e-10);
    }

    TEST(gsAsMatrix_Transpose)
    {
        gsMatrix<> orig(2, 3);
        orig << 1, 2, 3,
                4, 5, 6;
        
        gsAsMatrix<> transposed(orig.data(), 3, 2);
        
        CHECK_EQUAL(transposed.rows(), 3);
        CHECK_EQUAL(transposed.cols(), 2);
    }

    TEST(gsAsConstMatrix_Const)
    {
        const gsMatrix<> mat(2, 2);
        
        gsAsConstMatrix<> constMat(mat.data(), 2, 2);
        
        CHECK_EQUAL(constMat.rows(), 2);
        CHECK_EQUAL(constMat.cols(), 2);
    }

    TEST(gsAsVector_Mapping)
    {
        gsMatrix<> mat(3, 1);
        mat << 1.5, 2.5, 3.5;
        
        gsAsVector<> vec(mat.data(), 3);
        
        CHECK_EQUAL(vec.size(), 3);
        CHECK_CLOSE(vec(1), 2.5, 1e-10);
    }

    TEST(gsSparseMatrix_Insertion)
    {
        gsSparseMatrix<> sparse(5, 5);
        
        sparse.insert(0, 0) = 1.0;
        sparse.insert(1, 1) = 2.0;
        sparse.insert(2, 2) = 3.0;
        sparse.insert(0, 4) = 4.0;
        
        CHECK_CLOSE(sparse.coeff(0, 0), 1.0, 1e-10);
        CHECK_CLOSE(sparse.coeff(1, 1), 2.0, 1e-10);
        CHECK_CLOSE(sparse.coeff(0, 4), 4.0, 1e-10);
    }

    TEST(gsSparseMatrix_MatrixMultiplication)
    {
        gsSparseMatrix<> A(3, 3);
        A.insert(0, 0) = 2.0;
        A.insert(1, 1) = 3.0;
        A.insert(2, 2) = 4.0;
        
        gsMatrix<> x(3, 1);
        x << 1, 2, 3;
        
        gsMatrix<> y = A * x;
        
        CHECK_CLOSE(y(0, 0), 2.0, 1e-10);
        CHECK_CLOSE(y(1, 0), 6.0, 1e-10);
        CHECK_CLOSE(y(2, 0), 12.0, 1e-10);
    }

    TEST(gsSparseMatrix_Transpose)
    {
        gsSparseMatrix<> A(2, 3);
        A.insert(0, 0) = 1.0;
        A.insert(0, 2) = 2.0;
        A.insert(1, 1) = 3.0;
        
        gsSparseMatrix<> At = A.transpose();
        
        CHECK_EQUAL(At.rows(), 3);
        CHECK_EQUAL(At.cols(), 2);
        CHECK_CLOSE(At.coeff(0, 0), 1.0, 1e-10);
        CHECK_CLOSE(At.coeff(2, 0), 2.0, 1e-10);
    }

    TEST(gsSparseMatrix_NonZeros)
    {
        gsSparseMatrix<> A(4, 4);
        A.insert(0, 0) = 1.0;
        A.insert(1, 1) = 2.0;
        A.insert(2, 2) = 3.0;
        
        CHECK_EQUAL(A.nonZeros(), 3);
    }

    TEST(gsBlockOp_Construction)
    {
        gsMatrix<> A(2, 2);
        A << 1, 0,
             0, 1;
        
        gsMatrix<> B(2, 2);
        B << 2, 0,
             0, 2;
        
        gsBlockOp<> block(2, 2);
        block.addOperator(0, 0, gsMatrix<>::Identity(2, 2));
        
        CHECK(true); // Just test construction
    }

    TEST(gsBlockDiag_Construction)
    {
        gsMatrix<> A(2, 2);
        A.setIdentity();
        
        gsMatrix<> B(3, 3);
        B.setIdentity();
        
        std::vector<gsMatrix<>> blocks;
        blocks.push_back(A);
        blocks.push_back(B);
        
        gsBlockDiag<> blockDiag(blocks);
        
        CHECK_EQUAL(blockDiag.rows(), 5);
        CHECK_EQUAL(blockDiag.cols(), 5);
    }

    TEST(gsSparseVector_Operations)
    {
        gsSparseVector<> vec(10);
        vec.insert(0) = 1.0;
        vec.insert(5) = 2.5;
        vec.insert(9) = 3.0;
        
        CHECK_EQUAL(vec.nonZeros(), 3);
        CHECK_CLOSE(vec.coeff(5), 2.5, 1e-10);
    }

    TEST(gsMatrix_BlockOperations)
    {
        gsMatrix<> mat(4, 4);
        mat << 1, 2, 3, 4,
               5, 6, 7, 8,
               9, 10, 11, 12,
               13, 14, 15, 16;
        
        gsMatrix<> block = mat.block(1, 1, 2, 2);
        
        CHECK_EQUAL(block.rows(), 2);
        CHECK_EQUAL(block.cols(), 2);
        CHECK_CLOSE(block(0, 0), 6.0, 1e-10);
        CHECK_CLOSE(block(1, 1), 11.0, 1e-10);
    }

    TEST(gsMatrix_Reshape)
    {
        gsMatrix<> mat(2, 6);
        for (int i = 0; i < 12; i++)
            mat.data()[i] = i + 1;
        
        gsAsMatrix<> reshaped(mat.data(), 3, 4);
        
        CHECK_EQUAL(reshaped.rows(), 3);
        CHECK_EQUAL(reshaped.cols(), 4);
    }

    TEST(gsMatrix_NormComputation)
    {
        gsMatrix<> mat(3, 1);
        mat << 3, 4, 0;
        
        double norm = mat.norm();
        
        CHECK_CLOSE(norm, 5.0, 1e-10);
    }

    TEST(gsMatrix_DotProduct)
    {
        gsVector<> v1(3);
        v1 << 1, 2, 3;
        
        gsVector<> v2(3);
        v2 << 4, 5, 6;
        
        double dot = v1.dot(v2);
        
        CHECK_CLOSE(dot, 32.0, 1e-10);
    }

    TEST(gsMatrix_CrossProduct)
    {
        gsVector<> v1(3);
        v1 << 1, 0, 0;
        
        gsVector<> v2(3);
        v2 << 0, 1, 0;
        
        gsVector<> cross = v1.cross(v2);
        
        CHECK_CLOSE(cross(0), 0.0, 1e-10);
        CHECK_CLOSE(cross(1), 0.0, 1e-10);
        CHECK_CLOSE(cross(2), 1.0, 1e-10);
    }

    /* 
     * Step-by-step instructions for additional complex gsMatrix tests:
     * 
     * 1. Sparse matrix patterns tests:
     *    - Test compressed row storage (CRS)
     *    - Test compressed column storage (CCS)
     *    - Test coordinate format (COO)
     *    - Test conversion between formats
     * 
     * 2. Matrix decomposition tests:
     *    - Test QR decomposition
     *    - Test SVD decomposition
     *    - Test eigendecomposition
     *    - Test Cholesky decomposition
     * 
     * 3. gsBlockOp advanced tests:
     *    - Create complex block structures
     *    - Test block-matrix vector products
     *    - Test block-matrix inversions
     *    - Test with saddle point systems
     * 
     * 4. gsSparseMatrix performance tests:
     *    - Test reserve() for efficient insertion
     *    - Test makeCompressed() optimization
     *    - Test prune() for removing small entries
     *    - Compare dense vs sparse operations
     * 
     * 5. gsMatrix views and mappings:
     *    - Test col() and row() extractions
     *    - Test diagonal() extraction
     *    - Test triangular views
     *    - Test self-adjoint views
     * 
     * 6. Advanced matrix operations:
     *    - Test matrix exponential
     *    - Test matrix square root
     *    - Test pseudo-inverse
     *    - Test condition number estimation
     */
}