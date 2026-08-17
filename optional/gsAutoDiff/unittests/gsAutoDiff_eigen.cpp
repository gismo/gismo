#include "gismo_unittest.h"
#include <gsAutoDiff/gsAutoDiff.h>
#include <gsCore/gsLinearAlgebra.h>

using namespace gismo;

SUITE(gsAutoDiff_Eigen)
{

// ============================================================================
// Test 1: Basic element access and assignment
// ============================================================================
TEST(ElementAccess)
{
    using T = var_t;
    
    gsMatrix<T> M(3, 3);
    for (index_t i = 0; i < 9; ++i)
        M(i/3, i%3) = T(static_cast<double>(i));
    
    // Test 1a: Reading elements preserves value
    CHECK_CLOSE(static_cast<double>(M(0,0)), 0.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(M(1,1)), 4.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(M(2,2)), 8.0, 1e-10);
    
    // Test 1b: Copy element and verify it's connected to computation graph
    T x = M(1, 1);
    T y = x * 2.0;
    auto grads = autodiff::derivatives(y, autodiff::reverse::detail::wrt(M(1, 1)));
    CHECK_CLOSE(grads[0], 2.0, 1e-10);
}

// ============================================================================
// Test 2: Matrix arithmetic operations
// ============================================================================
TEST(MatrixArithmetic)
{
    using T = var_t;
    
    gsMatrix<T> A(2, 2);
    A(0,0) = T(1.0); A(0,1) = T(2.0);
    A(1,0) = T(3.0); A(1,1) = T(4.0);
    
    // Test 2a: Scalar multiplication
    gsMatrix<T> B = A * 2.0;
    CHECK_CLOSE(static_cast<double>(B(0,0)), 2.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(B(1,1)), 8.0, 1e-10);
    
    // Test 2b: Addition
    gsMatrix<T> C = A + A;
    CHECK_CLOSE(static_cast<double>(C(0,0)), 2.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(C(1,1)), 8.0, 1e-10);
    
    // Test 2c: Derivatives through scalar multiplication
    T result = B(0,0) * B(0,1);  // (1*2) * (2*2) = 4 * 4 = 16
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(0,0)));
    // d( (2*A00) * (2*A01) ) / d(A00) = 2 * (2*A01) = 2 * 4 = 8
    CHECK_CLOSE(grads[0], 8.0, 1e-10);
}

// ============================================================================
// Test 3: Matrix multiplication (THE PROBLEMATIC CASE)
// ============================================================================
TEST(MatrixMultiplication)
{
    using T = var_t;
    
    gsMatrix<T> A(2, 2);
    A(0,0) = T(1.0); A(0,1) = T(0.0);
    A(1,0) = T(0.0); A(1,1) = T(1.0);
    
    gsMatrix<T> B(2, 2);
    B(0,0) = T(2.0); B(0,1) = T(3.0);
    B(1,0) = T(4.0); B(1,1) = T(5.0);
    
    // Test 3a: Standard multiplication gives correct values
    gsMatrix<T> C = A * B;
    CHECK_CLOSE(static_cast<double>(C(0,0)), 2.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(C(0,1)), 3.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(C(1,0)), 4.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(C(1,1)), 5.0, 1e-10);
    
    // Test 3b: Derivatives through multiplication
    T result = C(0,0);  // Should be A(0,0)*B(0,0) + A(0,1)*B(1,0) = 1*2 + 0*4 = 2
    
    // Check if derivative w.r.t. B(0,0) exists
    auto grads_B = autodiff::derivatives(result, autodiff::reverse::detail::wrt(B(0,0)));
    double deriv_B = grads_B[0];
    
    // Check if derivative w.r.t. A(0,0) exists
    auto grads_A = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(0,0)));
    double deriv_A = grads_A[0];
    
    // These will fail if Eigen multiplication breaks derivatives
    CHECK_CLOSE(deriv_B, 1.0, 1e-10);  // C[0,0] = A[0,0]*B[0,0] + ... so dC/dB[0,0] = A[0,0] = 1
    CHECK_CLOSE(deriv_A, 2.0, 1e-10);  // dC/dA[0,0] = B[0,0] = 2
}

// ============================================================================
// Test 4: Transpose operations
// ============================================================================
TEST(TransposeOperations)
{
    using T = var_t;
    
    gsMatrix<T> A(2, 3);
    A(0,0) = T(1.0); A(0,1) = T(2.0); A(0,2) = T(3.0);
    A(1,0) = T(4.0); A(1,1) = T(5.0); A(1,2) = T(6.0);
    
    // Test 4a: Transpose values
    gsMatrix<T> At = A.transpose();
    CHECK_CLOSE(static_cast<double>(At(0,0)), 1.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(At(1,0)), 2.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(At(0,1)), 4.0, 1e-10);
    
    // Test 4b: Derivatives through transpose
    T result = At(0,1) * At(1,0);  // 4.0 * 2.0 = 8.0
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(1,0)));  // A(1,0) = At(0,1) = 4
    // d(At[0,1]*At[1,0])/d(A[1,0]) = d(4*2)/d(4) = 2 (derivative of second factor)
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 2.0, 1e-10);
}

// ============================================================================
// Test 5: Block operations
// ============================================================================
TEST(BlockOperations)
{
    using T = var_t;
    
    gsMatrix<T> A(4, 4);
    for (index_t i = 0; i < 4; ++i)
        for (index_t j = 0; j < 4; ++j)
            A(i,j) = T(static_cast<double>(i*4 + j));
    
    // Test 5a: Block extraction
    auto block = A.block(1, 1, 2, 2);
    CHECK_CLOSE(static_cast<double>(block(0,0)), 5.0, 1e-10);  // A(1,1)
    CHECK_CLOSE(static_cast<double>(block(1,1)), 10.0, 1e-10); // A(2,2)
    
    // Test 5b: Derivatives through block
    T result = block(0,0) * block(1,1);  // 5 * 10 = 50
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(1,1)));
    // d(5*10)/d(5) = 10
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 10.0, 1e-10);
}

// ============================================================================
// Test 6: Row and column operations
// ============================================================================
TEST(RowColumnOperations)
{
    using T = var_t;
    
    gsMatrix<T> A(3, 3);
    for (index_t i = 0; i < 9; ++i)
        A(i/3, i%3) = T(static_cast<double>(i));
    
    // Test 6a: Row access
    auto row = A.row(1);
    CHECK_CLOSE(static_cast<double>(row(0)), 3.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(row(1)), 4.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(row(2)), 5.0, 1e-10);
    
    // Test 6b: Column access
    auto col = A.col(1);
    CHECK_CLOSE(static_cast<double>(col(0)), 1.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(col(1)), 4.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(col(2)), 7.0, 1e-10);
    
    // Test 6c: Derivatives through row/column
    T result = row(0) + col(1);  // 3 + 4 = 7
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(1,0)));
    // d(3+4)/d(3) = 1
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 1.0, 1e-10);
}

// ============================================================================
// Test 7: Array operations vs Matrix operations
// ============================================================================
TEST(ArrayVsMatrix)
{
    using T = var_t;
    
    gsMatrix<T> A(2, 2);
    A(0,0) = T(1.0); A(0,1) = T(2.0);
    A(1,0) = T(3.0); A(1,1) = T(4.0);
    
    gsMatrix<T> B(2, 2);
    B(0,0) = T(2.0); B(0,1) = T(1.0);
    B(1,0) = T(4.0); B(1,1) = T(3.0);
    
    // Test 7a: Array-wise multiplication (element-wise)
    auto C_array = A.array() * B.array();
    CHECK_CLOSE(static_cast<double>(C_array(0,0)), 2.0, 1e-10);  // 1*2
    CHECK_CLOSE(static_cast<double>(C_array(1,1)), 12.0, 1e-10); // 4*3
    
    // Test 7b: Derivatives through array ops
    T result = C_array(0,0) + C_array(1,1);  // 2 + 12 = 14
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(A(0,0)));
    // d(2 + 12)/d(1) = d(1*2)/d(1) + d(4*3)/d(1) = 2 + 0 = 2
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 2.0, 1e-10);
}

// ============================================================================
// Test 8: Cwise operations
// ============================================================================
TEST(CwiseOperations)
{
    using T = var_t;
    
    gsMatrix<T> A(2, 2);
    A(0,0) = T(1.0); A(0,1) = T(2.0);
    A(1,0) = T(3.0); A(1,1) = T(4.0);
    
    // Test 8a: Cwise unary operations
    auto B = A.cwiseAbs();
    CHECK_CLOSE(static_cast<double>(B(0,0)), 1.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(B(1,1)), 4.0, 1e-10);
    
    // Test 8b: Derivatives through cwise operations
    gsMatrix<T> C(2, 2);
    C(0,0) = T(2.0); C(0,1) = T(-1.0);
    C(1,0) = T(3.0); C(1,1) = T(-2.0);
    
    auto D = C.cwiseAbs();
    T result = D(0,0) + D(0,1);  // 2 + 1 = 3
    auto grads = autodiff::derivatives(result, autodiff::reverse::detail::wrt(C(0,0)));
    // d(2 + 1)/d(2) = 1
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 1.0, 1e-10);
}

// ============================================================================
// Test 9: Matrix-vector multiplication
// ============================================================================
TEST(MatrixVectorMultiplication)
{
    using T = var_t;
    
    gsMatrix<T> A(3, 3);
    A(0,0) = T(1.0); A(0,1) = T(0.0); A(0,2) = T(0.0);
    A(1,0) = T(0.0); A(1,1) = T(1.0); A(1,2) = T(0.0);
    A(2,0) = T(0.0); A(2,1) = T(0.0); A(2,2) = T(1.0);
    
    gsVector<T> v(3);
    v(0) = T(1.0); v(1) = T(2.0); v(2) = T(3.0);
    
    // Test 9a: Matrix-vector product
    gsVector<T> result = A * v;
    CHECK_CLOSE(static_cast<double>(result(0)), 1.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(result(1)), 2.0, 1e-10);
    CHECK_CLOSE(static_cast<double>(result(2)), 3.0, 1e-10);
    
    // Test 9b: Derivatives
    T sum = result(0) + result(1);  // 1 + 2 = 3
    auto grads = autodiff::derivatives(sum, autodiff::reverse::detail::wrt(v(0)));
    // d(1 + 2)/d(v[0]) = d(A[0,0]*v[0] + A[1,1]*v[1])/d(v[0]) = A[0,0] = 1
    double deriv = grads[0];
    CHECK_CLOSE(deriv, 1.0, 1e-10);
}

// ============================================================================
// Test 10: Comparison with dual_t (forward mode)
// ============================================================================
TEST(ComparisonWithDual)
{
    // Forward mode should work better since derivatives are stored in the type
    using T = dual_t;
    
    gsMatrix<T> A(2, 2);
    A(0,0) = T(1.0); A(0,0).grad = 1.0;
    A(0,1) = T(2.0);
    A(1,0) = T(3.0);
    A(1,1) = T(4.0);
    
    gsMatrix<T> B(2, 2);
    B(0,0) = T(2.0);
    B(0,1) = T(1.0);
    B(1,0) = T(4.0);
    B(1,1) = T(3.0);
    
    // Matrix multiplication with dual_t
    gsMatrix<T> C = A * B;
    
    // Check derivative in C[0,0]
    // C[0,0] = A[0,0]*B[0,0] + A[0,1]*B[1,0]
    // dC[0,0]/dA[0,0] = B[0,0] = 2.0
    double deriv = C(0,0).grad;
    CHECK_CLOSE(deriv, 2.0, 1e-10);
}

}
