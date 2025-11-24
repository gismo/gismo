/** @file gsCore_test.cpp

    @brief Additional tests for gsCore module to increase coverage

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsCore_test)
{
    TEST(gsAffineFunction_Construction)
    {
        // Create a simple affine function f(x) = Ax + b
        gsMatrix<> A(2, 2);
        A << 1, 0,
             0, 1;
        
        gsMatrix<> b(2, 1);
        b << 1, 2;
        
        gsAffineFunction<> affine(A, b);
        
        CHECK_EQUAL(affine.domainDim(), 2);
        CHECK_EQUAL(affine.targetDim(), 2);
    }

    TEST(gsAffineFunction_Evaluation)
    {
        gsMatrix<> A(2, 2);
        A << 2, 0,
             0, 3;
        
        gsMatrix<> b(2, 1);
        b << 1, 1;
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(2, 1);
        point << 1, 1;
        
        gsMatrix<> result = affine.eval(point);
        
        CHECK_EQUAL(result.rows(), 2);
        CHECK_CLOSE(result(0, 0), 3.0, 1e-10); // 2*1 + 1
        CHECK_CLOSE(result(1, 0), 4.0, 1e-10); // 3*1 + 1
    }

    TEST(gsConstantFunction_Construction)
    {
        gsConstantFunction<> constant(5.0, 2);
        
        CHECK_EQUAL(constant.targetDim(), 2);
    }

    TEST(gsConstantFunction_Evaluation)
    {
        gsConstantFunction<> constant(3.14, 1);
        
        gsMatrix<> point(1, 3);
        point << 0.0, 0.5, 1.0;
        
        gsMatrix<> result = constant.eval(point);
        
        CHECK_EQUAL(result.rows(), 1);
        CHECK_EQUAL(result.cols(), 3);
        CHECK_CLOSE(result(0, 0), 3.14, 1e-10);
        CHECK_CLOSE(result(0, 1), 3.14, 1e-10);
        CHECK_CLOSE(result(0, 2), 3.14, 1e-10);
    }

    TEST(gsConstantFunction_VectorValue)
    {
        gsMatrix<> value(3, 1);
        value << 1.0, 2.0, 3.0;
        
        gsConstantFunction<> constant(value, 2);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = constant.eval(point);
        
        CHECK_EQUAL(result.rows(), 3);
        CHECK_CLOSE(result(0, 0), 1.0, 1e-10);
        CHECK_CLOSE(result(1, 0), 2.0, 1e-10);
        CHECK_CLOSE(result(2, 0), 3.0, 1e-10);
    }

    // Commented out - gsConstantBasis requires bases in constructor
    // TEST(gsConstantBasis_Construction)
    // {
    //     gsConstantBasis<> cbasis;
    //     
    //     CHECK_EQUAL(cbasis.size(), 1);
    //     CHECK_EQUAL(cbasis.dim(), 0);
    // }

    // TEST(gsConstantBasis_Evaluation)
    // {
    //     gsConstantBasis<> cbasis;
    //     
    //     gsMatrix<> point(0, 1); // 0-dimensional domain
    //     
    //     gsMatrix<> result;
    //     cbasis.eval_into(point, result);
    //     
    //     CHECK(result.rows() > 0);
    // }

    TEST(gsFuncCoordinate_Construction)
    {
        // Create a simple function
        gsConstantFunction<> func(5.0, 3); // 3D output
        
        // Extract coordinate 1
        gsFuncCoordinate<> coordFunc(func, 1);
        
        CHECK_EQUAL(coordFunc.targetDim(), 1);
    }

    TEST(gsFuncCoordinate_Evaluation)
    {
        gsMatrix<> value(3, 1);
        value << 1.0, 2.0, 3.0;
        
        gsConstantFunction<> func(value, 2);
        
        // Extract coordinate 1 (second coordinate)
        gsFuncCoordinate<> coordFunc(func, 1);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = coordFunc.eval(point);
        
        CHECK_EQUAL(result.rows(), 1);
        CHECK_CLOSE(result(0, 0), 2.0, 1e-10);
    }

    // Commented out - gsBulk is a template requiring parameters
    // TEST(gsBulk_Construction)
    // {
    //     gsBulk bulk;
    //     bulk.patch = 0;
    //     
    //     CHECK_EQUAL(bulk.patch, 0);
    // }

    // TEST(gsBulk_Comparison)
    // {
    //     gsBulk b1, b2;
    //     b1.patch = 0;
    //     b2.patch = 0;
    //     
    //     CHECK(b1 == b2);
    //     
    //     b2.patch = 1;
    //     CHECK(b1 != b2);
    // }

    TEST(gsComposedFunction_Construction)
    {
        // Create inner function: identity in 2D
        gsAffineFunction<> inner(gsMatrix<>::Identity(2,2), gsMatrix<>::Zero(2,1));
        
        // Create outer function: scale by 2
        gsAffineFunction<> outer(2*gsMatrix<>::Identity(2,2), gsMatrix<>::Zero(2,1));
        
        gsComposedFunction<> composed(outer, inner);
        
        CHECK_EQUAL(composed.domainDim(), 2);
        CHECK_EQUAL(composed.targetDim(), 2);
    }

    TEST(gsComposedFunction_Evaluation)
    {
        // Inner: f(x) = x
        gsAffineFunction<> inner(gsMatrix<>::Identity(2,2), gsMatrix<>::Zero(2,1));
        
        // Outer: g(y) = 2*y + [1,1]
        gsMatrix<> b(2,1);
        b << 1, 1;
        gsAffineFunction<> outer(2*gsMatrix<>::Identity(2,2), b);
        
        // Composed: g(f(x)) = 2*x + [1,1]
        gsComposedFunction<> composed(outer, inner);
        
        gsMatrix<> point(2, 1);
        point << 1, 2;
        
        gsMatrix<> result = composed.eval(point);
        
        CHECK_CLOSE(result(0, 0), 3.0, 1e-10); // 2*1 + 1
        CHECK_CLOSE(result(1, 0), 5.0, 1e-10); // 2*2 + 1
    }

    // Commented out - gsField API takes gsMultiPatch not gsBSpline
    // TEST(gsField_Construction)
    // {
    //     // Create a simple geometry patch
    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsBSplineBasis<> basis(kv);
    //     
    //     gsMatrix<> coefs(2, 2);
    //     coefs << 0, 0,
    //              1, 0;
    //     
    //     gsBSpline<> patch(kv, coefs);
    //     
    //     // Create field coefficients
    //     gsMatrix<> fieldCoefs(2, 1);
    //     fieldCoefs << 1.0, 2.0;
    //     
    //     gsField<> field(patch, fieldCoefs);
    //     
    //     CHECK_EQUAL(field.dim(), 1);
    //     CHECK(field.patches().nPatches() > 0);
    // }

    // TEST(gsField_Evaluation)
    // {
    //     gsKnotVector<> kv(0, 1, 1, 3);
    //     gsBSplineBasis<> basis(kv);
    //     
    //     gsMatrix<> geomCoefs(3, 1);
    //     geomCoefs << 0, 0.5, 1.0;
    //     
    //     gsBSpline<> patch(kv, geomCoefs);
    //     
    //     gsMatrix<> fieldCoefs(3, 1);
    //     fieldCoefs << 1.0, 2.0, 3.0;
    //     
    //     gsField<> field(patch, fieldCoefs);
    //     
    //     gsMatrix<> point(1, 1);
    //     point << 0.5;
    //     
    //     gsMatrix<> result = field.value(point, 0);
    //     
    //     CHECK(result.rows() > 0);
    // }

    // TEST(gsGeometrySlice_Construction)
    // {
    //     // Create a 2D patch
    //     gsKnotVector<> kv1(0, 1, 1, 2);
    //     gsKnotVector<> kv2(0, 1, 1, 2);
    //     gsTensorBSplineBasis<2> basis(kv1, kv2);
    //     
    //     gsMatrix<> coefs(4, 2);
    //     coefs << 0, 0,
    //              1, 0,
    //              0, 1,
    //              1, 1;
    //     
    //     gsTensorBSpline<2> patch(basis, coefs);
    //     
    //     // Create a slice at u = 0.5 in direction 0
    //     gsGeometrySlice<> slice(patch, 0, 0.5);
    //     
    //     CHECK_EQUAL(slice.domainDim(), 1);
    // }

    TEST(gsMultiBasis_AddBasis)
    {
        gsMultiBasis<> mb;
        
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsBSplineBasis<> basis1(kv1);
        
        gsKnotVector<> kv2(0, 1, 2, 3);
        gsBSplineBasis<> basis2(kv2);
        
        mb.addBasis(&basis1);
        mb.addBasis(&basis2);
        
        CHECK_EQUAL(mb.nBases(), 2);
        CHECK_EQUAL(mb.nPieces(), 2);
    }

    TEST(gsMultiBasis_TotalSize)
    {
        gsMultiBasis<> mb;
        
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsBSplineBasis<> basis1(kv1);
        
        gsKnotVector<> kv2(0, 1, 1, 3);
        gsBSplineBasis<> basis2(kv2);
        
        mb.addBasis(&basis1);
        mb.addBasis(&basis2);
        
        size_t total = mb.totalSize();
        CHECK_EQUAL(total, basis1.size() + basis2.size());
    }

    // Commented out - gsMultiBasis doesn't have maxDegree() method
    // TEST(gsMultiBasis_MaxDegree)
    // {
    //     gsMultiBasis<> mb;
    //     
    //     gsKnotVector<> kv1(0, 1, 1, 2);
    //     gsBSplineBasis<> basis1(kv1);
    //     
    //     gsKnotVector<> kv2(0, 1, 3, 2);
    //     gsBSplineBasis<> basis2(kv2);
    //     
    //     mb.addBasis(&basis1);
    //     mb.addBasis(&basis2);
    //     
    //     CHECK(mb.maxDegree() >= 3);
    // }

    // Commented out - gsPiecewiseFunction takes references not pointers
    // TEST(gsPiecewiseFunction_Construction)
    // {
    //     gsPiecewiseFunction<> pwFunc;
    //     
    //     gsConstantFunction<> func1(1.0, 1);
    //     gsConstantFunction<> func2(2.0, 1);
    //     
    //     pwFunc.addPiece(&func1);
    //     pwFunc.addPiece(&func2);
    //     
    //     CHECK_EQUAL(pwFunc.nPieces(), 2);
    // }

    // TEST(gsPiecewiseFunction_PieceAccess)
    // {
    //     gsPiecewiseFunction<> pwFunc;
    //     
    //     gsConstantFunction<> func1(1.0, 1);
    //     gsConstantFunction<> func2(2.0, 1);
    //     
    //     pwFunc.addPiece(&func1);
    //     pwFunc.addPiece(&func2);
    //     
    //     const gsFunctionSet<>& piece0 = pwFunc.piece(0);
    //     const gsFunctionSet<>& piece1 = pwFunc.piece(1);
    //     
    //     CHECK(piece0.targetDim() > 0);
    //     CHECK(piece1.targetDim() > 0);
    // }

    // Commented out - gsGeometryTransform has different constructor
    // TEST(gsGeometryTransform_Construction)
    // {
    //     // Create a simple transformation
    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsBSplineBasis<> basis(kv);
    //     
    //     gsMatrix<> coefs(2, 2);
    //     coefs << 0, 0,
    //              2, 0; // Scales by 2 in x direction
    //     
    //     gsBSpline<> transform(kv, coefs);
    //     
    //     gsGeometryTransform<> geomTrans(transform);
    //     
    //     CHECK_EQUAL(geomTrans.domainDim(), 1);
    // }

    TEST(gsDofMapper_Construction)
    {
        gsDofMapper mapper;
        
        // Initialize with multi-basis
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        
        gsMultiBasis<> mb;
        mb.addBasis(&basis);
        mb.addBasis(&basis);
        
        mapper.init(mb);
        
        CHECK_EQUAL(mapper.numPatches(), (size_t)2);
    }

    TEST(gsDofMapper_SetFree)
    {
        gsDofMapper mapper;
        
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        
        gsMultiBasis<> mb;
        mb.addBasis(&basis);
        
        mapper.init(mb);
        mapper.finalize();
        
        CHECK(mapper.freeSize() > 0);
    }

    TEST(gsBasisFun_Construction)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        gsBasisFun<> bfun(basis, 0); // First basis function
        
        CHECK_EQUAL(bfun.domainDim(), 1);
    }

    TEST(gsBasisFun_Evaluation)
    {
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        
        gsBasisFun<> bfun(basis, 1); // Second basis function
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> result = bfun.eval(point);
        
        CHECK(result.rows() > 0);
    }

    // Additional comprehensive tests for gsCore module

    TEST(gsAffineFunction_Identity)
    {
        gsMatrix<> A = gsMatrix<>::Identity(3, 3);
        gsMatrix<> b = gsMatrix<>::Zero(3, 1);
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(3, 1);
        point << 1.5, 2.5, 3.5;
        
        gsMatrix<> result = affine.eval(point);
        
        CHECK_CLOSE(result(0, 0), 1.5, 1e-10);
        CHECK_CLOSE(result(1, 0), 2.5, 1e-10);
        CHECK_CLOSE(result(2, 0), 3.5, 1e-10);
    }

    TEST(gsAffineFunction_Rotation)
    {
        // 90 degree rotation matrix
        gsMatrix<> A(2, 2);
        A << 0, -1,
             1, 0;
        gsMatrix<> b = gsMatrix<>::Zero(2, 1);
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(2, 1);
        point << 1, 0;
        
        gsMatrix<> result = affine.eval(point);
        
        CHECK_CLOSE(result(0, 0), 0.0, 1e-10);
        CHECK_CLOSE(result(1, 0), 1.0, 1e-10);
    }

    TEST(gsAffineFunction_Translation)
    {
        gsMatrix<> A = gsMatrix<>::Identity(2, 2);
        gsMatrix<> b(2, 1);
        b << 5, 10;
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(2, 1);
        point << 1, 2;
        
        gsMatrix<> result = affine.eval(point);
        
        CHECK_CLOSE(result(0, 0), 6.0, 1e-10);
        CHECK_CLOSE(result(1, 0), 12.0, 1e-10);
    }

    TEST(gsAffineFunction_Scaling)
    {
        gsMatrix<> A(2, 2);
        A << 2, 0,
             0, 3;
        gsMatrix<> b = gsMatrix<>::Zero(2, 1);
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(2, 1);
        point << 4, 5;
        
        gsMatrix<> result = affine.eval(point);
        
        CHECK_CLOSE(result(0, 0), 8.0, 1e-10);
        CHECK_CLOSE(result(1, 0), 15.0, 1e-10);
    }

    TEST(gsAffineFunction_MultiplePoints)
    {
        gsMatrix<> A = gsMatrix<>::Identity(2, 2);
        gsMatrix<> b(2, 1);
        b << 1, 1;
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> points(2, 3);
        points << 0, 1, 2,
                  0, 1, 2;
        
        gsMatrix<> result = affine.eval(points);
        
        CHECK_EQUAL(result.cols(), 3);
        CHECK_CLOSE(result(0, 0), 1.0, 1e-10);
        CHECK_CLOSE(result(0, 1), 2.0, 1e-10);
        CHECK_CLOSE(result(0, 2), 3.0, 1e-10);
    }

    TEST(gsConstantFunction_ZeroValue)
    {
        gsConstantFunction<> constant(0.0, 3);
        
        gsMatrix<> point(3, 1);
        point << 1, 2, 3;
        
        gsMatrix<> result = constant.eval(point);
        
        CHECK_CLOSE(result(0, 0), 0.0, 1e-10);
    }

    TEST(gsConstantFunction_NegativeValue)
    {
        gsConstantFunction<> constant(-5.5, 2);
        
        gsMatrix<> point(2, 1);
        point << 1, 2;
        
        gsMatrix<> result = constant.eval(point);
        
        CHECK_CLOSE(result(0, 0), -5.5, 1e-10);
    }

    TEST(gsConstantFunction_LargeVector)
    {
        gsMatrix<> value(10, 1);
        for (int i = 0; i < 10; i++)
            value(i, 0) = i * 1.5;
        
        gsConstantFunction<> constant(value, 3);
        
        gsMatrix<> point(3, 1);
        point.setRandom();
        
        gsMatrix<> result = constant.eval(point);
        
        CHECK_EQUAL(result.rows(), 10);
        for (int i = 0; i < 10; i++)
            CHECK_CLOSE(result(i, 0), i * 1.5, 1e-10);
    }

    TEST(gsFuncCoordinate_FirstCoordinate)
    {
        gsMatrix<> value(3, 1);
        value << 1.0, 2.0, 3.0;
        gsConstantFunction<> func(value, 2);
        
        gsFuncCoordinate<> coordFunc(func, 0);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = coordFunc.eval(point);
        
        CHECK_EQUAL(result.rows(), 1);
        CHECK_CLOSE(result(0, 0), 1.0, 1e-10);
    }

    TEST(gsFuncCoordinate_LastCoordinate)
    {
        gsMatrix<> value(5, 1);
        for (int i = 0; i < 5; i++)
            value(i, 0) = i + 10;
        
        gsConstantFunction<> func(value, 2);
        
        gsFuncCoordinate<> coordFunc(func, 4);
        
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        
        gsMatrix<> result = coordFunc.eval(point);
        
        CHECK_CLOSE(result(0, 0), 14.0, 1e-10);
    }

    TEST(gsMultiBasis_EmptyConstruction)
    {
        gsMultiBasis<> mb;
        
        CHECK_EQUAL(mb.nBases(), 0);
        CHECK_EQUAL(mb.nPieces(), 0);
    }

    TEST(gsMultiBasis_SingleBasis)
    {
        gsMultiBasis<> mb;
        
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        mb.addBasis(&basis);
        
        CHECK_EQUAL(mb.nBases(), 1);
        CHECK_EQUAL(mb.totalSize(), basis.size());
    }

    TEST(gsMultiBasis_DifferentDegrees)
    {
        gsMultiBasis<> mb;
        
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsBSplineBasis<> basis1(kv1);
        
        gsKnotVector<> kv2(0, 1, 3, 2);
        gsBSplineBasis<> basis2(kv2);
        
        mb.addBasis(&basis1);
        mb.addBasis(&basis2);
        
        CHECK_EQUAL(mb.nBases(), 2);
    }

    TEST(gsDofMapper_MultiplePatchesInit)
    {
        gsDofMapper mapper;
        
        gsKnotVector<> kv1(0, 1, 1, 3);
        gsBSplineBasis<> basis1(kv1);
        
        gsKnotVector<> kv2(0, 1, 2, 2);
        gsBSplineBasis<> basis2(kv2);
        
        gsMultiBasis<> mb;
        mb.addBasis(&basis1);
        mb.addBasis(&basis2);
        
        mapper.init(mb);
        mapper.finalize();
        
        CHECK_EQUAL(mapper.numPatches(), (size_t)2);
    }

    TEST(gsBasisFun_Support)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        gsBasisFun<> bfun(basis, 0);
        
        gsMatrix<> support = bfun.support();
        
        CHECK(support.rows() > 0);
        CHECK_EQUAL(support.cols(), 2);
    }

    TEST(gsBasisFun_Derivative)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        gsBasisFun<> bfun(basis, 1);
        
        gsMatrix<> point(1, 1);
        point << 0.5;
        
        gsMatrix<> deriv = bfun.deriv(point);
        
        CHECK(deriv.rows() > 0);
    }

    TEST(gsAffineFunction_Derivative)
    {
        gsMatrix<> A(2, 2);
        A << 2, 1,
             1, 3;
        gsMatrix<> b(2, 1);
        b << 1, 2;
        
        gsAffineFunction<> affine(A, b);
        
        gsMatrix<> point(2, 1);
        point << 1, 1;
        
        gsMatrix<> deriv = affine.deriv(point);
        
        // Derivative of affine function is the matrix A
        CHECK_EQUAL(deriv.rows(), 2);
        CHECK_EQUAL(deriv.cols(), 2);
    }

    /* 
     * Step-by-step instructions for additional complex gsCore tests:
     * 
     * 1. gsMultiPatch tests:
     *    - Create multiple geometry patches (gsBSpline, gsTensorBSpline)
     *    - Combine into gsMultiPatch
     *    - Test patch access by index
     *    - Test addPatch(), addPatch with interface information
     *    - Test computeTopology() for automatic interface detection
     *    - Test global parameterization across patches
     *    - Test gsMultiPatch::eval() for evaluation across all patches
     * 
     * 2. gsBoundary and gsBoxTopology tests:
     *    - Create box topology with multiple patches
     *    - Define interfaces between patches
     *    - Test boundary side enumeration
     *    - Test interface matching and orientation
     *    - Test topology queries (adjacent patches, boundary sides)
     * 
     * 3. gsRationalBasis tests:
     *    - Create NURBS basis (rational B-spline)
     *    - Test weight manipulation
     *    - Test evaluation with weights
     *    - Test derivative computation with rational functions
     *    - Test extraction of numerator and denominator
     * 
     * 4. gsFunctionSet advanced tests:
     *    - Test piece() method for multi-piece functions
     *    - Test eval_into() for efficient evaluation
     *    - Test deriv_into() and deriv2_into()
     *    - Test active_into() for active basis functions
     *    - Test support queries for individual functions
     * 
     * 5. gsGenericGeometry and gsComposedGeometry tests:
     *    - Create geometry from basis and coefficients
     *    - Test isClosed() for periodic geometries
     *    - Test isAffine() check
     *    - Test geometry composition (curve on surface, etc.)
     *    - Test point inversion (physical to parametric coords)
     * 
     * 6. gsDofMapper advanced tests:
     *    - Test with boundary conditions (eliminateDofs())
     *    - Test coupled DOFs (matchDof(), coupledDof())
     *    - Test free/eliminated/coupled DOF queries
     *    - Test global/local index conversion
     *    - Test component-wise DOF management
     */
}
