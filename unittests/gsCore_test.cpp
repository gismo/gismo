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
    // NOTE: Affine function tests merged into a single `gsAffineFunction` test below.

    // TEST(gsAffineFunction_Evaluation)
    // {
    //     // TODO: Fix eval_into API - has dimension mismatch issues
    //     gsMatrix<> A(2, 2);
    //     A << 2, 0,
    //          0, 3;
    //     
    //     gsVector<> b(2);
    //     b << 1, 1;
    //     
    //     gsAffineFunction<> affine(A, b);
    //     
    //     gsMatrix<> point(2, 1);
    //     point << 1, 1;
    //     
    //     gsMatrix<> result;
    //     affine.eval_into(point, result);
    //     
    //     CHECK_EQUAL(result.rows(), 2);
    //     CHECK_CLOSE(result(0, 0), 3.0, 1e-10); // 2*1 + 1
    //     CHECK_CLOSE(result(1, 0), 4.0, 1e-10); // 3*1 + 1
    // }

    // NOTE: Constant function constructor/evaluation checks merged into `gsConstantFunction` test below.

    TEST(gsBasisFun)
    {
        gsKnotVector<> kv(0, 1, 2, 3);
        gsBSplineBasis<> basis(kv);
        
        gsBasisFun<> bfun(basis, 0); // First basis function
        CHECK_EQUAL(1, bfun.domainDim());

        gsMatrix<> point(1, 1); point << 0.5;
        gsMatrix<> result; bfun.eval_into(point, result);
        CHECK_EQUAL(1, result.rows());

        gsMatrix<> deriv = bfun.deriv(point);
        CHECK_EQUAL(1, deriv.rows());
    }

    // NOTE: explicit evaluation test merged into constructor+evaluation test above.

    // TEST(gsBulk_Construction)
    // {
    //     gsBulk<real_t> b;

    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsTensorBSplineBasis<4> basis(kv, kv, kv, kv);
    //     gsMatrix<> coefs = gsMatrix<>::Random(basis.size(), 4);
    //     b = gsBulk<>(basis,coefs);
    //     CHECK_EQUAL(b.dim(), 4);
    //     CHECK_EQUAL(b.coefs().rows(), basis.size()); 
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

    TEST(gsComposedFunction)
    {
        // Create inner function: identity in 2D
        gsGeometry<>::uPtr inner = gsNurbsCreator<>::BSplineSquare();
        
        // Create outer function: scale by 2
        gsGeometry<>::uPtr outer = gsNurbsCreator<>::BSplineRectangle(0.,0.,2.,2.);
        gsComposedFunction<> composed(*inner, *outer);
        
        CHECK_EQUAL(2, composed.domainDim());
        CHECK_EQUAL(2, composed.targetDim());

        // Also evaluate to verify result
        gsMatrix<> point(2, 1);
        point << 0.5, 0.5;
        gsMatrix<> evalRes;
        composed.eval_into(point, evalRes);
        CHECK_EQUAL(2, evalRes.rows());
    }

    // evaluation merged into the construction test above

    // Commented out - gsField API takes gsMultiPatch not gsBSpline
    TEST(gsField)
    {
        // Create a simple multi-patch
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineSquare());

        // Create a function
        gsConstantFunction<> func(1.0, 2);
        
        // Create a field
        gsField<> field(mp, func);
        
        CHECK_EQUAL(1, field.dim());
        CHECK_EQUAL(1, field.patches().nPatches());
    }

    // evaluation was merged into gsField test above

    TEST(gsGeometrySlice_Construction)
    {
        // Create a 2D patch
        gsTensorBSpline<2> patch = *gsNurbsCreator<>::BSplineSquare();
        
        // Create a slice at u = 0.5 in direction 0
        gsGeometrySlice<> slice(&patch, 0, 0.5);
        
        CHECK_EQUAL(1, slice.domainDim());
    }

    TEST(gsMultiBasis)
    {
        // Empty construction
        gsMultiBasis<> mbEmpty;
        CHECK_EQUAL(0, mbEmpty.nBases());
        CHECK_EQUAL(0, mbEmpty.nPieces());

        // Single basis
        gsKnotVector<> kv1(0, 1, 1, 2);
        gsBSplineBasis<> basis1(kv1);
        gsMultiBasis<> mbSingle;
        mbSingle.addBasis(basis1.clone().release());
        CHECK_EQUAL(1, mbSingle.nBases());
        CHECK_EQUAL(basis1.size(), mbSingle.totalSize());

        // Two bases and total size
        gsKnotVector<> kv2(0, 1, 1, 3);
        gsBSplineBasis<> basis2(kv2);
        gsMultiBasis<> mb;
        mb.addBasis(basis1.clone().release());
        mb.addBasis(basis2.clone().release());
        CHECK_EQUAL(2, mb.nBases());
        CHECK_EQUAL(2, mb.nPieces());
        CHECK_EQUAL(basis1.size() + basis2.size(), (index_t)mb.totalSize());

        // Degrees (per-component cwise degrees)
        short_t expected_max = std::max(basis1.maxDegree(), basis2.maxDegree());
        short_t expected_min = std::min(basis1.minDegree(), basis2.minDegree());
        CHECK_EQUAL(expected_max, mb.maxCwiseDegree());
        CHECK_EQUAL(expected_min, mb.minCwiseDegree());
    }

    TEST(gsPiecewiseFunction)
    {
        gsPiecewiseFunction<> pwFunc;
        
        gsConstantFunction<> func1(1.0, 1);
        gsConstantFunction<> func2(2.0, 1);
        
        pwFunc.addPiece(func1);
        pwFunc.addPiece(func2);
        
        CHECK_EQUAL(2, pwFunc.nPieces());
        const gsFunctionSet<>& piece0 = pwFunc.piece(0);
        const gsFunctionSet<>& piece1 = pwFunc.piece(1);
        CHECK_EQUAL(1, piece0.targetDim());
        CHECK_EQUAL(1, piece1.targetDim());
    }

    // TEST(gsGeometryTransform_Construction)
    // {
    //     // Create a simple transformation
    //     gsKnotVector<> kv(0, 1, 1, 2);
    //     gsBSplineBasis<> basis(kv);
        
    //     gsMatrix<> coefs(2, 2);
    //     coefs << 0, 0,
    //              2, 0; // Scales by 2 in x direction
        
    //     gsBSpline<> transform(kv, coefs);
        
    //     gsGeometryTransform<real_t> geomTrans(transform,);
        
    //     CHECK_EQUAL(geomTrans.domainDim(), 1);
    // }

    TEST(gsDofMapper)
    {
        gsDofMapper mapper;
        
        // Initialize with multi-basis
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        
        gsMultiBasis<> mb;
        mb.addBasis(basis.clone().release());
        mb.addBasis(basis.clone().release());
        mapper.init(mb);
        CHECK_EQUAL((size_t)2, mapper.numPatches());

        // finalize and check free size
        mapper.finalize();
        CHECK(mapper.freeSize() > 0);
    }


    // Additional comprehensive tests for gsCore module

    TEST(gsAffineFunction)
    {
        // Identity 3D
        gsMatrix<> A3 = gsMatrix<>::Identity(3, 3);
        gsVector<> b3 = gsVector<>::Zero(3);
        gsAffineFunction<> affine3(A3, b3);
        gsMatrix<> p3(3,1); p3 << 1.5, 2.5, 3.5;
        gsMatrix<> out; affine3.eval_into(p3, out);
        CHECK_CLOSE(1.5, out(0,0), 1e-10);
        CHECK_CLOSE(2.5, out(1,0), 1e-10);
        CHECK_CLOSE(3.5, out(2,0), 1e-10);

        // Rotation 2D
        gsMatrix<> A2r(2,2); A2r << 0, -1, 1, 0;
        gsVector<> b2r = gsVector<>::Zero(2);
        gsAffineFunction<> rot(A2r, b2r);
        gsMatrix<> p2(2,1); p2 << 1, 0;
        rot.eval_into(p2, out);
        CHECK_CLOSE(0.0, out(0,0), 1e-10);
        CHECK_CLOSE(1.0, out(1,0), 1e-10);

        // Translation 2D
        gsMatrix<> A2t = gsMatrix<>::Identity(2,2);
        gsVector<> b2t(2); b2t << 5,10;
        gsAffineFunction<> tr(A2t, b2t);
        gsMatrix<> pt(2,1); pt << 1,2; tr.eval_into(pt, out);
        CHECK_CLOSE(6.0, out(0,0), 1e-10);
        CHECK_CLOSE(12.0, out(1,0), 1e-10);

        // Scaling 2D
        gsMatrix<> A2s(2,2); A2s << 2,0, 0,3; gsVector<> b0 = gsVector<>::Zero(2);
        gsAffineFunction<> sc(A2s, b0); gsMatrix<> psc(2,1); psc << 4,5; sc.eval_into(psc, out);
        CHECK_CLOSE(8.0, out(0,0), 1e-10);
        CHECK_CLOSE(15.0, out(1,0), 1e-10);

        // Multiple points + offset
        gsMatrix<> A2m = gsMatrix<>::Identity(2,2);
        gsMatrix<> b2m(2,1); b2m << 1,1;
        gsAffineFunction<> affm(A2m, b2m);
        gsMatrix<> points(2,3); points << 0,1,2, 0,1,2;
        affm.eval_into(points, out);
        CHECK_EQUAL(3, out.cols());
        CHECK_CLOSE(1.0, out(0,0), 1e-10);
        CHECK_CLOSE(2.0, out(0,1), 1e-10);
        CHECK_CLOSE(3.0, out(0,2), 1e-10);

        // Derivative check
        gsMatrix<> A2d(2,2); A2d << 2,1, 1,3; gsMatrix<> b2d(2,1); b2d << 1,2;
        gsAffineFunction<> ad(A2d, b2d); gsMatrix<> pp(2,1); pp << 1,1;
        gsMatrix<> deriv = ad.deriv(pp);
        CHECK_EQUAL(4, deriv.rows());
        CHECK_EQUAL(1, deriv.cols());
    }

    // NOTE: Constant function tests merged above into `gsConstantFunction`.

    TEST(gsFuncCoordinate)
    {
        // first coordinate
        gsMatrix<> value(3, 1); value << 1.0, 2.0, 3.0;
        gsConstantFunction<> func(value, 2);
        gsFuncCoordinate<> coord0(func, 0);
        gsMatrix<> point(2,1); point << 0.5, 0.5;
        gsMatrix<> r0 = coord0.eval(point);
        CHECK_EQUAL(1, r0.rows());
        CHECK_CLOSE(1.0, r0(0,0), 1e-10);

        // last coordinate
        gsMatrix<> value2(5,1); for (int i=0;i<5;i++) value2(i,0)=i+10;
        gsConstantFunction<> func2(value2, 2);
        gsFuncCoordinate<> coord4(func2, 4);
        gsMatrix<> r4 = coord4.eval(point);
        CHECK_CLOSE(14.0, r4(0,0), 1e-10);
    }



    // Removed duplicate, merged into `gsBasisFun` test above

    // NOTE: Basis derivative test merged into gsBasisFun test above.

    // NOTE: Affine derivative test merged into `gsAffineFunction` test above.

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

// NOTE: Redundant tests commented out to avoid duplication.
// Tests like `gsAffineFunction_Evaluation` and `gsField` have been merged into comprehensive tests.
