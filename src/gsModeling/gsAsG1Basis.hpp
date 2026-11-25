/** @file gsAsG1Basis.hpp

    @brief Derive a C1 basis for an AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

namespace gismo {

// ====================================================================
// Helper utilities (matrix manipulations algebra)
// ====================================================================

template <typename T>
gsSparseMatrix<T> asEmbeddingMatrix(index_t rows, const gsMatrix<index_t> &c)
{
    index_t cols = c.rows();
    GISMO_ASSERT(c.cols() == 1, "");
    gsSparseEntries<T> se;
    se.reserve(cols);
    for (index_t i = 0; i < c.size(); ++i)
        se.add(c(i, 0), i, 1);
    gsSparseMatrix<T> result(rows, cols);
    result.setFrom(se);
    return result;
}

// ====================================================================
// Helper utilities (basis evaluation)
// ====================================================================

template <typename T>
bool isNested(const gsBSplineBasis<T> &coarseBasis,
              const gsBSplineBasis<T> &fineBasis)
{
    if (coarseBasis.dim() != 1 || fineBasis.dim() != 1)
        return false;
    const int coarseDegree = coarseBasis.degree();
    const int fineDegree = fineBasis.degree();
    if (coarseDegree > fineDegree)
        return false;
    const gsKnotVector<T> &coarseKnots = coarseBasis.knots();
    const gsKnotVector<T> &fineKnots = fineBasis.knots();
    gsKnotVector<T> elevated = coarseKnots;
    elevated.degreeElevate(fineDegree - coarseDegree);
    return fineKnots.includes(elevated);
}

template <typename T>
gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T> &basis,
                                    const gsMatrix<T> &pts,
                                    const index_t derivative = 0)
{
    gsMatrix<T> vals;
    if (derivative == 0)
        vals = basis.eval(pts);
    else if (derivative == 1)
        vals = basis.deriv(pts);
    else
        GISMO_ERROR("Invalid value.");
    gsMatrix<index_t> idx = basis.active(pts);
    gsSparseEntries<T> entries;
    entries.reserve(vals.rows() * vals.cols());
    for (int i = 0; i < vals.cols(); ++i)
        for (int j = 0; j < vals.rows(); ++j)
            entries.add(i, idx(j, i), vals(j, i));
    gsSparseMatrix<T> C(pts.cols(), basis.size());
    C.setFrom(entries);
    return C;
}

template <typename T>
gsSparseMatrix<T> evaluateBasisAtBoundary(const gsBSplineBasis<T> &basis,
                                          index_t sideParameter,
                                          index_t derivative,
                                          index_t layers)
{
    GISMO_ASSERT(sideParameter == 0 || sideParameter == 1, "Invalid side.");

    const gsMatrix<T> pt = basis.support().col(sideParameter);
    gsMatrix<index_t> indices = basis.active(pt);

    gsMatrix<T> values;
    if (derivative == 0)
        values = basis.eval(pt);
    else if (derivative == 1)
    {
        values = basis.deriv(pt);
        // On left hand side, normal vector is -1
        if (sideParameter == 0)
            values *= -1;
    }
    else
        GISMO_ERROR("Invalid value.");

    // Flip indices if right side is selected
    if (sideParameter == 1)
        indices.array() = basis.size() - 1 - indices.array();

    // Store in sparse matrix
    gsSparseMatrix<T> result(1, layers);
    for (index_t i = 0; i < indices.rows(); ++i)
        if (indices(i, 0) < layers)
            result(0, indices(i, 0)) = values(i, 0);

    return result;
}

// ====================================================================
// Evaluation of tensor basis
// ====================================================================

template <typename T>
gsSparseMatrix<T> collocateBoundaryValues(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    boxSide side,
    const gsMatrix<T> &pts)
{
    const gsBSplineBasis<T> normalBasis = tensorBasis.component(side.direction());
    const gsBSplineBasis<T> tangentialBasis = *tensorBasis.boundaryBasis(side);

    return evaluateBasisAtBoundary(normalBasis, side.parameter(), 0, 2)
        .kron(collocationMatrix(tangentialBasis, pts));
}

template <typename T>
gsSparseMatrix<T> collocateBoundaryCrossingDerivative(
      const gsTensorBSplineBasis<2, T> &tensorBasis,
      boxSide side,
      const gsMatrix<T> &gluingData,
      const gsMatrix<T> &pts)
{
    const gsBSplineBasis<T> normalBasis = tensorBasis.component(side.direction());
    const gsBSplineBasis<T> tangentialBasis = *tensorBasis.boundaryBasis(side);

    const gsMatrix<T> tangentialBasisSupport = tangentialBasis.support();
    const T a = tangentialBasisSupport(0, 0);
    const T b = tangentialBasisSupport(0, 1);
    const gsMatrix<T> ptsNormalized = (pts.transpose().array() - a) / (b - a);
    const gsMatrix<T> alpha = gluingData(0, 0) * (1 - ptsNormalized.array())
                              + gluingData(0, 1) * ptsNormalized.array();
    const gsMatrix<T> beta = gluingData(0, 2) * (1 - ptsNormalized.array())
                             + gluingData(0, 3) * ptsNormalized.array();

    const T sgn = side.parameter() == side.direction() ? 1 : -1; // This is as in the paper
    const gsMatrix<T> scaling0 = beta.array() / alpha.array() * sgn;
    const gsMatrix<T> scaling1 = 1 / alpha.array();

    return evaluateBasisAtBoundary(normalBasis, side.parameter(), 0, 2)
              .kron(scaling0.asDiagonal() * collocationMatrix(tangentialBasis, pts, 1))
          + evaluateBasisAtBoundary(normalBasis, side.parameter(), 1, 2)
              .kron(scaling1.asDiagonal() * collocationMatrix(tangentialBasis, pts, 0));
}

template <typename T>
gsMatrix<T> reducedHessAsMatrix(const gsMatrix<T>& reduced)
{
    GISMO_ENSURE (reduced.rows()==1 && reduced.cols()==3, "Wrong dimensions");
    gsMatrix<T> matrix(2,2);
    matrix(0,0) = reduced(0,0); //d_xx
    matrix(0,1) = reduced(0,2); //d_xy
    matrix(1,0) = reduced(0,2); //d_yx
    matrix(1,1) = reduced(0,1); //d_yy
    return matrix;
}

template <typename T>
gsMatrix<T> hessMatrixAsReduced(const gsMatrix<T>& matrix)
{
    GISMO_ENSURE (matrix.rows()==2 && matrix.cols()==2, "Wrong dimensions");
    gsMatrix<T> reduced(1,3);
    reduced(0,0) = matrix(0,0); //d_xx
    reduced(0,1) = matrix(1,1); //d_yy
    reduced(0,2) = (matrix(1,0) + matrix(0,1))/2; //d_xy=d_yx
    return reduced;
}


template <typename T>
gsSparseMatrix<T> collocateCorners(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    const gsGeometry<T> &geo,
    const gsMatrix<T> &normals)
{
    const gsMatrix<T> support = tensorBasis.support();

    // Coordinates of corners in parameter domain
    gsMatrix<T> corners(2, 4);
    for (index_t i=0; i<4; ++i)
    {
        corners(0, i) = support(0, i%2);
        corners(1, i) = support(1, i/2);
    }

    // Evaluate at corners
    gsMatrix<index_t> idx = tensorBasis.active(corners); //  n x 4
    gsMatrix<T> vals = tensorBasis.eval(corners);        //  n x 4
    gsMatrix<T> der1 = tensorBasis.deriv(corners);       // 2n x 4
    gsMatrix<T> der2 = tensorBasis.deriv2(corners);      // 3n x 4

    // Reserve matrix for result
    gsMatrix<T> result(4*6, tensorBasis.size());
    result.setZero();

    // For each corner
    const index_t stride = normals.cols() / 4;
    GISMO_ENSURE(normals.rows() == 1 &&
                 (normals.cols() == 8 || normals.cols() == 16 || normals.cols() == 20 ||
                  normals.cols() == 28 || normals.cols() == 40),
                 "Invalid corner directions matrix: expected 1x8, 1x16, 1x20, 1x28, or 1x40, got "
                 << normals.rows() << "x" << normals.cols());

    if (stride == 7 || stride == 10)
    {
        // 3D surface corner 6-jets in the shared 2D vertex tangent-plane frame:
        //   stride == 7 : [d1(3), d2(3), kappa_g]
        //   stride == 10: [d1(3), d2(3), k0_e1, k0_e2, k1_e1, k1_e2]
        GISMO_ENSURE(geo.targetDim() == 3,
                     "stride == 7 or 10 corner frames require a 3D surface geometry (targetDim == 3).");
        gsMatrix<T> gDeriv, gDeriv2;
        geo.deriv_into(corners, gDeriv);   // 6 x 4
        geo.deriv2_into(corners, gDeriv2); // 9 x 4

        for (index_t i = 0; i < 4; ++i)
        {
            gsVector<T> d1(3), d2(3);
            d1 << normals(0, stride*i + 0), normals(0, stride*i + 1), normals(0, stride*i + 2);
            d2 << normals(0, stride*i + 3), normals(0, stride*i + 4), normals(0, stride*i + 5);
            const T cKappa = (stride == 7) ? normals(0, stride*i + 6) : T(0);
            const T k0_e1  = (stride == 10) ? normals(0, stride*i + 6) : T(0);
            const T k0_e2  = (stride == 10) ? normals(0, stride*i + 7) : T(0);
            const T k1_e1  = (stride == 10) ? normals(0, stride*i + 8) : T(0);
            const T k1_e2  = (stride == 10) ? normals(0, stride*i + 9) : T(0);

            // Orthonormal basis E = [e1, e2] of span(d1, d2) = T_v S
            gsVector<T> e1 = d1.normalized();
            gsVector<T> e2 = d2 - d2.dot(e1) * e1;
            GISMO_ENSURE(e2.norm() > T(1e-12), "Degenerate 3D corner tangent frame at corner " << i);
            e2.normalize();

            gsMatrix<T> E(3, 2);
            E.col(0) = e1;
            E.col(1) = e2;

            gsMatrix<T> D(3, 2);
            D.col(0) = d1;
            D.col(1) = d2;
            const gsMatrix<T> Tmat = E.transpose() * D; // 2x2: columns are d1, d2 in E-basis

            // Geometry Jacobian J (3x2) and projected 2x2 Jacobian Jhat = E^T J
            gsMatrix<T> J(3, 2);
            for (short_t c = 0; c < 3; ++c)
            {
                J(c, 0) = gDeriv(2*c + 0, i);
                J(c, 1) = gDeriv(2*c + 1, i);
            }
            const gsMatrix<T> Jhat = E.transpose() * J;
            const gsMatrix<T> JhatInv = Jhat.inverse();
            const gsMatrix<T> JhatInvT = JhatInv.transpose();

            // Projected parametric second derivatives H_{y_0}, H_{y_1} (2x2 each)
            gsMatrix<T> Hy0(2, 2), Hy1(2, 2);
            gsVector<T> guu(3), gvv(3), guv(3);
            for (short_t c = 0; c < 3; ++c)
            {
                guu(c) = gDeriv2(3*c + 0, i);
                gvv(c) = gDeriv2(3*c + 1, i);
                guv(c) = gDeriv2(3*c + 2, i);
            }
            Hy0(0, 0) = e1.dot(guu); Hy0(0, 1) = e1.dot(guv);
            Hy0(1, 0) = Hy0(0, 1);   Hy0(1, 1) = e1.dot(gvv);
            Hy1(0, 0) = e2.dot(guu); Hy1(0, 1) = e2.dot(guv);
            Hy1(1, 0) = Hy1(0, 1);   Hy1(1, 1) = e2.dot(gvv);

            for (index_t j = 0; j < idx.rows(); ++j)
            {
                gsVector<T> gXi(2);
                gXi << der1(2*j + 0, i), der1(2*j + 1, i);
                gsMatrix<T> HXi(2, 2);
                HXi(0, 0) = der2(3*j + 0, i);
                HXi(1, 1) = der2(3*j + 1, i);
                HXi(0, 1) = der2(3*j + 2, i);
                HXi(1, 0) = HXi(0, 1);

                const gsVector<T> gy = JhatInvT * gXi;
                const gsMatrix<T> Hy = JhatInvT * (HXi - gy(0) * Hy0 - gy(1) * Hy1) * JhatInv;

                const gsVector<T> gd = Tmat.transpose() * gy;
                const gsMatrix<T> Hd = Tmat.transpose() * Hy * Tmat;

                result(6*i + 0, idx(j, i)) = vals(j, i);
                result(6*i + 1, idx(j, i)) = gd(0);
                result(6*i + 2, idx(j, i)) = gd(1);
                if (stride == 10)
                {
                    result(6*i + 3, idx(j, i)) = Hd(0, 0) + k0_e1 * gy(0) + k0_e2 * gy(1);
                    result(6*i + 4, idx(j, i)) = Hd(1, 1) + k1_e1 * gy(0) + k1_e2 * gy(1);
                }
                else
                {
                    result(6*i + 3, idx(j, i)) = Hd(0, 0);
                    result(6*i + 4, idx(j, i)) = Hd(1, 1) + cKappa * gd(0);
                }
                result(6*i + 5, idx(j, i)) = Hd(0, 1);
            }
        }
        return result.sparseView(1, T(1e-12));
    }

    // Transform first and second derivatives to physical domain using map data (2D path)
    gsMapData<T> md(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER);
    md.points = corners;
    geo.computeMap(md);

    for (index_t i=0; i<4; ++i)
    {
        // Transform gradient to physical domain
        gsMatrix<T> der1Phys;
        transformGradients(md, i, der1, der1Phys);
        der1Phys = der1Phys.transpose();

        // Transform Hessian to physical domain
        gsMatrix<T> der2Phys;
        transformDeriv2Hgrad(md, i, der1, der2, der2Phys);

        if (stride == 2)
        {
            // Transform x and y derivatives into derivatives in normal and tangential direction
            gsMatrix<T> transform(2,2);
            transform(0,0) = normals(0+2*i);
            transform(1,0) = normals(1+2*i);
            transform(0,1) = normals(1+2*i);
            transform(1,1) = -normals(0+2*i);
            der1Phys = der1Phys * transform;

            for (index_t j=0; j<idx.rows(); ++j)
                der2Phys.row(j) = hessMatrixAsReduced<T>(
                        transform * reducedHessAsMatrix<T>(der2Phys.row(j)) * transform
                );

            for (index_t j=0; j<idx.rows(); ++j)
            {
                result(6*i+0, idx(j, i)) = vals(j, i);       //      u
                result(6*i+1, idx(j, i)) = der1Phys(j, 0);   // d_n  u
                result(6*i+2, idx(j, i)) = der1Phys(j, 1);   // d_t  u
                result(6*i+3, idx(j, i)) = der2Phys(j, 0);   // d_nn u
                result(6*i+4, idx(j, i)) = der2Phys(j, 1);   // d_tt u
                result(6*i+5, idx(j, i)) = der2Phys(j, 2);   // d_nt u
            }
        }
        else
        {
            // General local frame transform Q = [d1, d2] (Case B: [n, t], Case A: [t1, t2])
            gsMatrix<T> transform(2,2);
            transform(0,0) = normals(0+stride*i);
            transform(1,0) = normals(1+stride*i);
            transform(0,1) = normals(2+stride*i);
            transform(1,1) = normals(3+stride*i);
            der1Phys = der1Phys * transform;

            for (index_t j=0; j<idx.rows(); ++j)
                der2Phys.row(j) = hessMatrixAsReduced<T>(
                        transform.transpose() * reducedHessAsMatrix<T>(der2Phys.row(j)) * transform
                );

            const T cKappa = (stride == 5) ? normals(4+stride*i) : T(0);
            for (index_t j=0; j<idx.rows(); ++j)
            {
                result(6*i+0, idx(j, i)) = vals(j, i);                              // u
                result(6*i+1, idx(j, i)) = der1Phys(j, 0);                          // d_{d1} u
                result(6*i+2, idx(j, i)) = der1Phys(j, 1);                          // d_{d2} u
                result(6*i+3, idx(j, i)) = der2Phys(j, 0);                          // d_{d1 d1} u
                result(6*i+4, idx(j, i)) = der2Phys(j, 1) + cKappa * der1Phys(j, 0);// d_{d2 d2} u + c_k d_{d1} u
                result(6*i+5, idx(j, i)) = der2Phys(j, 2);                          // d_{d1 d2} u
            }
        }
    }
    return result.sparseView(1, T(1e-12));
}

// ====================================================================
// Derive various embeddings
// ====================================================================

template <typename T> struct gsEdgeEmbedding {
    gsSparseMatrix<T> matrix;
    gsVector<index_t, 2> sizes;
};

template <typename T>
gsEdgeEmbedding<T> deriveEdgeEmbedding(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    const gsMatrix<T> &localGluingData,
    boxSide side,
    T eps = 1e-12)
{
    const gsBSplineBasis<T> sideBasis = *tensorBasis.boundaryBasis(side);
    const gsMatrix<T> greville = sideBasis.anchors();

    const gsSparseMatrix<T> collocateLeft = gsBlockSparseMatrix<T>(2, 1)
            .set(0, 0, collocateBoundaryValues(tensorBasis, side, greville))
            .set(1, 0, collocateBoundaryCrossingDerivative(tensorBasis, side,
                                                    localGluingData, greville));

    gsBSplineBasis<T> sideLowerDegreeBasis = sideBasis;
    sideLowerDegreeBasis.degreeReduce(1);

    gsBSplineBasis<T> sideSmootherBasis = sideBasis;
    sideSmootherBasis.elevateContinuity(1);

    const gsSparseMatrix<T> collocateRight = gsBlockSparseMatrix<T>(2, 2)
            .set(0, 0, collocationMatrix(sideSmootherBasis, greville))
            .set(1, 1, collocationMatrix(sideLowerDegreeBasis, greville));

    const gsLinearOperator<>::Ptr solver = makeSparseLUSolver(collocateLeft);
    gsMatrix<T> result;
    solver->apply(collocateRight, result);

    gsEdgeEmbedding<T> ee;
    ee.matrix = result.sparseView(1, eps);
    ee.sizes[0] = sideSmootherBasis.size();
    ee.sizes[1] = sideLowerDegreeBasis.size();
    return ee;
}

template <typename T>
gsSparseMatrix<T> deriveInnerPortionEdgeEmbedding(const gsEdgeEmbedding<T> &ee)
{
    const index_t rows = ee.matrix.cols();
    const index_t cols = ee.sizes[0] + ee.sizes[1] - 10;
    GISMO_ENSURE(cols > 0,
                "Inner portion edge embedding is not possible: " << cols);
    gsVector<index_t> indices(cols);
    for (index_t i = 3; i < ee.sizes[0] - 3; ++i)
        indices[i - 3] = i;

    for (index_t i = 2; i < ee.sizes[1] - 2; ++i)
        indices[ee.sizes[0] - 8 + i] = ee.sizes[0] + i;

    return ee.matrix * asEmbeddingMatrix<T>(rows, indices);
}

template <typename T>
gsSparseMatrix<T> deriveCornersPortionEdgeEmbedding(const gsEdgeEmbedding<T> &ee)
{
    const index_t rows = ee.matrix.cols();
    gsVector<index_t> indices(10);
    indices[0] = 0;
    indices[1] = 1;
    indices[2] = 2;
    indices[3] = ee.sizes[0] - 3;
    indices[4] = ee.sizes[0] - 2;
    indices[5] = ee.sizes[0] - 1;

    indices[6] = ee.sizes[0];
    indices[7] = ee.sizes[0] + 1;
    indices[8] = ee.sizes[0] + ee.sizes[1] - 2;
    indices[9] = ee.sizes[0] + ee.sizes[1] - 1;

    return ee.matrix * asEmbeddingMatrix<T>(rows, indices);
}

template <typename T>
std::vector<index_t> getInteriorDofs(const gsTensorBSplineBasis<2, T> &tensorBasis, boxSide side)
{
    std::vector<index_t> removedSet;
    for (index_t i = 0; i < 2; ++i)
    {
        gsMatrix<index_t> bdyDofs = tensorBasis.boundaryOffset(side, i);
        std::copy(bdyDofs.data(), bdyDofs.data() + bdyDofs.rows(),
                std::back_inserter(removedSet));
    }

    std::sort(removedSet.begin(), removedSet.end());

    std::vector<index_t> all(tensorBasis.size());
    std::iota(all.begin(), all.end(), 0);

    std::vector<index_t> interior;
    interior.reserve(all.size());
    std::set_difference(all.begin(), all.end(), removedSet.begin(),
                        removedSet.end(), std::back_inserter(interior));
    return interior;
}

template <typename T>
gsSparseMatrix<T> deriveInnerEmbedding(const gsTensorBSplineBasis<2, T> &tensorBasis, boxSide side)
{
    std::vector<index_t> interior = getInteriorDofs(tensorBasis, side);
    // Convert to gsVector
    gsVector<index_t> interiorVec;
    interiorVec.assign(interior.begin(), interior.end());
    // Create embedding matrix
    return asEmbeddingMatrix<T>(tensorBasis.size(), interiorVec);
}

template <typename T>
std::vector<index_t> getInteriorDofs(const gsTensorBSplineBasis<2, T> &tensorBasis)
{
    std::vector<index_t> removedSet;
    for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
    {
        for (index_t i = 0; i < 2; ++i)
        {
          gsMatrix<index_t> bdyDofs = tensorBasis.boundaryOffset(side, i);
          std::copy(bdyDofs.data(), bdyDofs.data() + bdyDofs.rows(),
                    std::back_inserter(removedSet));
        }
    }

    std::sort(removedSet.begin(), removedSet.end());

    std::vector<index_t> all(tensorBasis.size());
    std::iota(all.begin(), all.end(), 0);

    std::vector<index_t> interior;
    interior.reserve(all.size());
    std::set_difference(all.begin(), all.end(), removedSet.begin(),
                        removedSet.end(), std::back_inserter(interior));
    return interior;
}

template <typename T>
gsSparseMatrix<T> deriveInnerEmbedding(const gsTensorBSplineBasis<2, T> &tensorBasis)
{
    std::vector<index_t> interior = getInteriorDofs(tensorBasis);
    // Convert to gsVector
    gsVector<index_t> interiorVec;
    interiorVec.assign(interior.begin(), interior.end());
    // Create embedding matrix
    return asEmbeddingMatrix<T>(tensorBasis.size(), interiorVec);
}

// ====================================================================
// Assemble overall embedding matrix
// ====================================================================

template <typename T>
gsSparseMatrix<T> deriveCornerEmbedding(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    const gsGeometry<T> &geo,
    const gsMatrix<T> &localGluingData,
    const gsMatrix<T> &normals)
{
    index_t rows = tensorBasis.size();
    gsBlockSparseMatrix<T> collocation(6, 5);

    collocation.set(0, 0, collocateCorners(tensorBasis, geo, normals));
    collocation.set(1, 0, deriveInnerEmbedding(tensorBasis).transpose());

    for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
    {
        const gsSparseMatrix<T> simpleEdgeEmbedding = gsBlockSparseMatrix<T>(1, 2)
                .set(0, 0, asEmbeddingMatrix<T>(rows, tensorBasis.boundary(side)))
                .set(0, 1, asEmbeddingMatrix<T>(rows, tensorBasis.boundaryOffset(side, 1)));

        collocation.set(1 + side.m_index, 0, simpleEdgeEmbedding.transpose());
        collocation.set(1 + side.m_index, side.m_index,
            deriveCornersPortionEdgeEmbedding(deriveEdgeEmbedding(
                tensorBasis,
                gsMatrix<T>(localGluingData.middleCols(4 * (side.m_index - 1), 4)),
                side))
        );
    }

    gsSparseMatrix<T> collocationMatrix = collocation;

    gsMatrix<T> rhs(collocationMatrix.rows(), 24);
    rhs.setZero();
    for (index_t i = 0; i < 24; i++)
        rhs(i, i) = 1;

    gsMatrix<T> result;
    makeSparseLUSolver(collocationMatrix)->apply(rhs, result);

    return result.topRows(rows).sparseView(T(1e-12));
}

template <typename T> struct gsAsG1Embedding {
public:
    gsSparseMatrix<T> matrix;
    // Sizes of blocks:
    //   0...interior dofs
    //   1...edge dofs level 0 (function values)
    //   2...edge dofs level 1 (crossing derivatives)
    gsVector<index_t, 3> sizes;
};

template <typename T>
gsAsG1Embedding<T> deriveTwoPatchASG1Embedding(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    boxSide side,
    const gsMatrix<T> &localGluingData,
    T eps = 1e-12)
{
    gsAsG1Embedding<T> result;
    index_t rows = tensorBasis.size();

    const gsSparseMatrix<T> embeddingInterior =
        deriveInnerEmbedding(tensorBasis, side);
    result.sizes[0] = embeddingInterior.cols();

    const gsSparseMatrix<T> simpleEdgeEmbedding = gsBlockSparseMatrix<T>(1, 2)
            .set(0, 0, asEmbeddingMatrix<T>(rows, tensorBasis.boundary(side)))
            .set(0, 1, asEmbeddingMatrix<T>(rows, tensorBasis.boundaryOffset(side, 1)));

    const gsEdgeEmbedding<T> asG1edgeEmbedding =
        deriveEdgeEmbedding(tensorBasis, localGluingData, side);

    result.sizes[1] = asG1edgeEmbedding.sizes[0];
    result.sizes[2] = asG1edgeEmbedding.sizes[1];

    result.matrix = gsBlockSparseMatrix<T>(1, 2)
            .set(0, 0, embeddingInterior)
            .set(0, 1, simpleEdgeEmbedding * asG1edgeEmbedding.matrix);

    return result;
}

template <typename T> struct gsArgyrisEmbedding {
public:
    gsSparseMatrix<T> matrix;
    // Sizes of blocks:
    //   0...interior dofs
    //   1...edge dofs level 0 (function values)
    //   2...edge dofs level 1 (crossing derivatives)
    gsVector<index_t, 13> sizes;
};

template <typename T>
gsArgyrisEmbedding<T> deriveArgyrisBasisEmbedding(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    const gsMatrix<T> &localGluingData,
    const gsMatrix<T> &normals,
    gsGeometry<T> &geo,
    T eps = 1e-12)
{
    gsArgyrisEmbedding<T> result;
    gsBlockSparseMatrix<T> blockMatrix(1, 6);
    index_t rows = tensorBasis.size();

    gsSparseMatrix<T> embeddingInterior = deriveInnerEmbedding(tensorBasis);
    blockMatrix.set(0, 0, embeddingInterior);
    result.sizes[0] = embeddingInterior.cols();

    for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
    {
        const gsSparseMatrix<T> simpleEdgeEmbedding = gsBlockSparseMatrix<T>(1, 2)
                .set(0, 0, asEmbeddingMatrix<T>(rows, tensorBasis.boundary(side)))
                .set(0, 1, asEmbeddingMatrix<T>(rows, tensorBasis.boundaryOffset(side, 1)));

        const gsEdgeEmbedding<T> asG1edgeEmbedding = deriveEdgeEmbedding(
            tensorBasis,
            gsMatrix<T>(localGluingData.middleCols(4 * (side.m_index - 1), 4)),
            side
        );

        blockMatrix.set(0, side.m_index, simpleEdgeEmbedding * deriveInnerPortionEdgeEmbedding(asG1edgeEmbedding));

        result.sizes[1 + 2 * (side.m_index - 1)] = asG1edgeEmbedding.sizes[0] - 6;
        result.sizes[2 + 2 * (side.m_index - 1)] = asG1edgeEmbedding.sizes[1] - 4;
    }

    blockMatrix.set(0, 5, deriveCornerEmbedding(tensorBasis, geo, localGluingData, normals));
    result.sizes[9] = 6;
    result.sizes[10] = 6;
    result.sizes[11] = 6;
    result.sizes[12] = 6;

    result.matrix = blockMatrix;

    return result;
}

template <typename T>
gsArgyrisEmbedding<T> deriveArgyrisBasisEmbedding(
    const gsTensorBSplineBasis<2, T> &tensorBasis,
    const gsMatrix<T> &localGluingData,
    gsGeometry<T> &geo,
    T eps = 1e-12)
{
    gsMatrix<T> normals(1,4*2);
    normals.setZero();
    for (index_t j=0; j<4; ++j)
        normals(0,2*j)=1;
    return deriveArgyrisBasisEmbedding(tensorBasis, localGluingData, normals, geo, eps);
}

index_t sumUntil(const gsVector<index_t, 13> &vec, index_t until) {
    index_t sum = 0;
    for (index_t i = 0; i < until; ++i)
        sum += vec(i);
    return sum;
}

enum class cornerConditionType {
    all = 6,                   // eliminate {0,1,2,3,4,5}
    valuesNormals = 5,         // smooth clamped (or sharp clampedSide2): eliminate {0,1,2,4,5}, free {3}
    valuesSharp = 51,          // sharp both-values: eliminate {0,1,2,3,4}, free {5}
    valuesSmooth = 3,          // smooth values (or sharp valuesSide2): eliminate {0,2,4}, free {1,3,5}
    derivatives = 2,           // smooth derivatives (or sharp derivativesSide2): eliminate {1,5}
    valuesSide1 = 31,          // sharp side-0 values only: eliminate {0,1,3}, free {2,4,5}
    valuesSide2 = 32,          // sharp side-1 values only: eliminate {0,2,4}, free {1,3,5}
    clampedSide1 = 52,         // sharp side-0 clamped, side-1 free/orth-sym: eliminate {0,1,2,3,5}, free {4}
    clampedSide2 = 53,         // sharp side-1 clamped, side-0 free/orth-sym: eliminate {0,1,2,4,5}, free {3}
    valuesSide1DerSide2 = 41,  // sharp side-0 values, side-1 orth-sym: eliminate {0,1,3,5}, free {2,4}
    derSide1ValuesSide2 = 42,  // sharp side-0 orth-sym, side-1 values: eliminate {0,2,4,5}, free {1,3}
    derivativesSharp = 23,     // sharp both-sides orth-sym: eliminate {1,2,5}, free {0,3,4}
    derivativesSide1 = 21,     // sharp side-0 orth-sym only: eliminate {2,5}, free {0,1,3,4}
    derivativesSide2 = 22,     // sharp side-1 orth-sym only: eliminate {1,5}, free {0,2,3,4}
    none = 0
};

inline void eliminateCornerConditionDofs(
    gsDofMapper& mapper,
    index_t off_corner,
    index_t patch,
    cornerConditionType tp,
    index_t comp = 0)
{
    auto elim = [&](index_t j) { mapper.eliminateDof(off_corner + j, patch, comp); };
    switch (tp)
    {
        case cornerConditionType::all:
            for (index_t j = 0; j < 6; ++j) elim(j);
            break;
        case cornerConditionType::valuesNormals:
        case cornerConditionType::clampedSide2:
            elim(0); elim(1); elim(2); elim(4); elim(5);
            break;
        case cornerConditionType::clampedSide1:
            elim(0); elim(1); elim(2); elim(3); elim(5);
            break;
        case cornerConditionType::valuesSharp:
            for (index_t j = 0; j < 5; ++j) elim(j);
            break;
        case cornerConditionType::valuesSmooth:
        case cornerConditionType::valuesSide2:
            elim(0); elim(2); elim(4);
            break;
        case cornerConditionType::valuesSide1:
            elim(0); elim(1); elim(3);
            break;
        case cornerConditionType::valuesSide1DerSide2:
            elim(0); elim(1); elim(3); elim(5);
            break;
        case cornerConditionType::derSide1ValuesSide2:
            elim(0); elim(2); elim(4); elim(5);
            break;
        case cornerConditionType::derivativesSharp:
            elim(1); elim(2); elim(5);
            break;
        case cornerConditionType::derivativesSide1:
            elim(2); elim(5);
            break;
        case cornerConditionType::derivatives:
        case cornerConditionType::derivativesSide2:
            elim(1); elim(5);
            break;
        case cornerConditionType::none:
            break;
        default:
            GISMO_ENSURE(false, "Unknown cornerConditionType.");
            break;
    }
}

typedef std::pair<patchCorner, cornerConditionType> cornerCondition;
typedef std::vector<cornerCondition> cornerConditions;

template<typename T>
gsDofMapper makeMapperForArgyrisBasis(
    const gsMultiPatch<T>& mp,
    const std::vector<gsArgyrisEmbedding<T>>& argBasis,
    const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
    const cornerConditions& cc = cornerConditions())
{
    gsVector<index_t> patchDofSizes(argBasis.size());
    for (size_t i = 0; i < argBasis.size(); ++i)
        patchDofSizes[i] = argBasis[i].matrix.cols();
    gsDofMapper mapper(patchDofSizes);
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface &ifc = *it;
        const patchSide ps1 = ifc.first();
        const patchSide ps2 = ifc.second();
        const gsVector<index_t>& sz1 = argBasis[ps1.patch].sizes;
        const gsVector<index_t>& sz2 = argBasis[ps2.patch].sizes;

        const index_t nLvl0 = sz1[1 + 2 * (ps1.m_index - 1)];
        const index_t offLvl0_1 = sumUntil(sz1, 1 + 2 * (ps1.m_index - 1));
        const index_t offLvl0_2 = sumUntil(sz2, 1 + 2 * (ps2.m_index - 1));

        const index_t nLvl1 = sz1[2 + 2 * (ps1.m_index - 1)];
        const index_t offLvl1_1 = sumUntil(sz1, 2 + 2 * (ps1.m_index - 1));
        const index_t offLvl1_2 = sumUntil(sz2, 2 + 2 * (ps2.m_index - 1));

        GISMO_ASSERT (nLvl0 == sz2[1+2*(ps2.m_index-1)], "Dimension missmatch.");
        GISMO_ASSERT (nLvl1 == sz2[2+2*(ps2.m_index-1)], "Dimension missmatch.");

        // Check interface orientation: if the tangential directions
        // run in opposite directions, we need to reverse the DOF mapping
        // for patch 2's shared columns.
        const short_t tanDir1 = 1 - ps1.direction();
        const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

        for (index_t j1 = 0; j1 < nLvl0; ++j1)
        {
            // If flipped, DOF j1 on patch 1 corresponds to DOF (nLvl0-1-j1) on
            // patch 2.
            const index_t j2 = flipped ? nLvl0 - 1 - j1 : j1;
            mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2);
        }
        for (index_t j1 = 0; j1 < nLvl1; ++j1)
        {
            // If flipped, DOF j1 on patch 1 corresponds to DOF (nLvl1-1-j1) on
            // patch 2.
            const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
            mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2);
        }

        // Also match the corners
        std::vector<patchCorner> corners;
        ps1.getContainedCorners(2, corners);
        GISMO_ASSERT (corners.size() == 2, "Unexpected number of corners");
        for (index_t i=0; i<2; ++i)
        {
            const index_t c1 = corners[i].m_index - 1;
            const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;

            const index_t off_corner_1 = sumUntil(argBasis[ps1.patch].sizes, 9 + c1);
            const index_t off_corner_2 = sumUntil(argBasis[ps2.patch].sizes, 9 + c2);

            for (index_t j=0; j<6; ++j)
                mapper.matchDof(ps1.patch, off_corner_1+j, ps2.patch, off_corner_2+j);

        }
    }

    // Now, incorporate the boundary conditions
    auto eliminateEdgeBlock = [&](const patchSide &ps, index_t levelOffset) {
        const index_t nLvl = argBasis[ps.patch].sizes[levelOffset + 2 * (ps.m_index - 1)];
        const index_t offLvl = sumUntil(argBasis[ps.patch].sizes, levelOffset + 2 * (ps.m_index - 1));
        for (index_t j = 0; j < nLvl; ++j)
            mapper.eliminateDof(offLvl + j, ps.patch);
    };

    for (auto it = bc.begin("ValuesAndDerivatives"); it != bc.end("ValuesAndDerivatives"); ++it)
    {
        eliminateEdgeBlock(it->ps, 1);
        eliminateEdgeBlock(it->ps, 2);
    }

    for (auto it = bc.begin("Values"); it != bc.end("Values"); ++it)
        eliminateEdgeBlock(it->ps, 1);

    for (auto it = bc.begin("Derivatives"); it != bc.end("Derivatives"); ++it)
        eliminateEdgeBlock(it->ps, 2);

    for (size_t i=0; i<cc.size(); ++i)
    {
        const index_t c = cc[i].first.m_index - 1;
        const index_t off_corner = sumUntil(argBasis[cc[i].first.patch].sizes, 9 + c);
        eliminateCornerConditionDofs(mapper, off_corner, cc[i].first.patch, cc[i].second, 0);
    }

    mapper.finalize();
    return mapper;
}

// ====================================================================
// Non-isoparametric (p, k) polynomial multi-basis, 3-component vector
// mapper, derivative scaling, and Hasanova-Takacs-Takacs (2026)
// Section 3.4 boundary lifting & patch-local projector \Pi_{AS}.
// ====================================================================

/**
 * @brief Build a polynomial tensor-product B-spline gsMultiBasis in S_{p, k, Z}
 *        (1 <= k <= p - 2, interior knot multiplicity m = p - k) on a multi-patch
 *        domain mp (which may be either polynomial B-spline or rational NURBS,
 *        realizing Section 3.1 of Hasanova, S. Takacs, T. Takacs 2026).
 *
 * Ensures every patch has at least enough interior knots so that the inner edge
 * spaces S_{p, k+1, Z}^{inner} and S_{p-1, k, Z}^{inner} are non-empty, then
 * applies numRef dyadic refinements preserving regularity k (multiplicity p - k).
 */
template <typename T>
gsMultiBasis<T> buildAsG1PolynomialMultiBasis(
    const gsMultiPatch<T>& mp,
    short_t degreeP,
    short_t regularityK = -1,
    index_t numRef = 0,
    index_t minInnerKnots = -1)
{
    GISMO_ENSURE(degreeP >= 3, "AS-G1 Argyris space requires degree p >= 3, got " << degreeP);
    if (regularityK < 0)
        regularityK = degreeP - 2; // Maximal AS-G1 smoothness k = p - 2
    GISMO_ENSURE(regularityK >= 1 && regularityK <= degreeP - 2,
                 "Regularity k must satisfy 1 <= k <= p - 2, got k=" << regularityK << ", p=" << degreeP);

    const short_t mult = degreeP - regularityK; // interior knot multiplicity in S_{p, k, Z}
    // Inner edge embedding requires dim(S_{p,k+1,Z}) >= 7 or dim(S_{p,k+1,Z}) >= 6 && dim(S_{p-1,k,Z}) >= 5
    // so cols = dim(S_{p,k+1,Z}) + dim(S_{p-1,k,Z}) - 10 > 0 and both inner loops are valid.
    // dim(S_{p,k+1,Z}) = p + 1 + nInner * (mult - 1) >= 6 => nInner >= ceil((5 - p) / (mult - 1))
    // dim(S_{p-1,k,Z}) = p + nInner * (mult - 1) >= 4
    // And for cols > 0: 2p + 1 + 2*nInner*(mult - 1) > 10 => nInner >= ceil((10 - 2p) / (2*(mult - 1)))
    index_t reqInner = 1;
    while ((degreeP + 1 + reqInner * (mult - 1) < 6) ||
           (degreeP + reqInner * (mult - 1) < 4) ||
           (2 * degreeP + 1 + 2 * reqInner * (mult - 1) <= 10))
    {
        ++reqInner;
    }
    if (minInnerKnots > reqInner)
        reqInner = minInnerKnots;

    gsMultiBasis<T> dbasis;
    dbasis.setTopology(mp.topology());
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        auto makeAlignedKV = [&](short_t dir) -> gsKnotVector<T> {
            const gsBasis<T>* bPtr = &mp.patch(p).basis();
            if (const auto* rBasis = dynamic_cast<const gsRationalBasis<gsTensorBSplineBasis<2, T>>*>(bPtr))
                bPtr = &rBasis->source();
            const gsBSplineBasis<T>& gbb =
                dynamic_cast<const gsBSplineBasis<T>&>(bPtr->component(dir));
            const gsKnotVector<T>& gkv = gbb.knots();
            const short_t gDeg = gbb.degree();

            // Collect geometry breakpoints and their parametric continuity r_geo = gDeg - gMult
            std::vector<T> brks = gkv.breaks();
            std::vector<short_t> brkMult(brks.size(), mult);
            for (size_t i = 1; i + 1 < brks.size(); ++i)
            {
                const index_t gMult = gkv.multiplicity(brks[i]);
                const short_t gCont = static_cast<short_t>(std::max<index_t>(0, gDeg - gMult));
                // Preserve geometry continuity across original interior knots so pullback u o G
                // is never over-constrained across C^{gCont} geometry knots (always >= mult >= 2).
                const short_t neededMult = static_cast<short_t>(std::max<int>(mult, degreeP - gCont));
                brkMult[i] = std::min<short_t>(neededMult, static_cast<short_t>(degreeP - 1));
            }

            // Refine breakpoints by midpoint insertion until at least reqInner interior breaks exist
            while (static_cast<index_t>(brks.size()) - 2 < reqInner)
            {
                std::vector<T> nextB;
                std::vector<short_t> nextM;
                nextB.reserve(2 * brks.size() - 1);
                nextM.reserve(2 * brks.size() - 1);
                for (size_t i = 0; i + 1 < brks.size(); ++i)
                {
                    nextB.push_back(brks[i]);
                    nextM.push_back(brkMult[i]);
                    nextB.push_back(T(0.5) * (brks[i] + brks[i + 1]));
                    nextM.push_back(mult);
                }
                nextB.push_back(brks.back());
                nextM.push_back(brkMult.back());
                brks.swap(nextB);
                brkMult.swap(nextM);
            }

            std::vector<T> knots;
            for (short_t m = 0; m < degreeP + 1; ++m) knots.push_back(brks.front());
            for (size_t i = 1; i + 1 < brks.size(); ++i)
                for (short_t m = 0; m < brkMult[i]; ++m)
                    knots.push_back(brks[i]);
            for (short_t m = 0; m < degreeP + 1; ++m) knots.push_back(brks.back());
            return gsKnotVector<T>(degreeP, knots.begin(), knots.end());
        };
        dbasis.addBasis(new gsTensorBSplineBasis<2, T>(makeAlignedKV(0), makeAlignedKV(1)));
    }
    if (numRef > 0)
        dbasis.uniformRefine(numRef, mult);
    return dbasis;
}

/**
 * @brief Scale derivative DOF columns of gsArgyrisEmbedding in place:
 *        - Edge level-1 (crossing derivative) columns by invH = 1 / hScale
 *        - Corner 1st-derivative columns (j = 1, 2) by invH = 1 / hScale
 *        - Corner 2nd-derivative columns (j = 3, 4, 5) by invH^2 = 1 / hScale^2
 *        Returns the per-patch column scaling vector so unscaled physical jets
 *        can be converted via d_scaled(k) = d_unscaled(k) / colScale[p](k).
 */
template <typename T>
std::vector<gsVector<T>> scaleArgyrisDerivativeBlocks(
    std::vector<gsArgyrisEmbedding<T>>& argBasis,
    T hScale)
{
    const T invH = (hScale > T(0)) ? (T(1) / hScale) : T(1);
    const T invH2 = invH * invH;
    std::vector<gsVector<T>> scales(argBasis.size());

    for (size_t p = 0; p < argBasis.size(); ++p)
    {
        const index_t nCols = argBasis[p].matrix.cols();
        scales[p].setOnes(nCols);
        const gsVector<index_t, 13>& sz = argBasis[p].sizes;

        // Edge level-1 blocks: indices 2, 4, 6, 8
        for (index_t s = 0; s < 4; ++s)
        {
            const index_t blk = 2 + 2 * s;
            const index_t off = sumUntil(sz, blk);
            for (index_t j = 0; j < sz[blk]; ++j)
                scales[p](off + j) = invH;
        }
        // Corner blocks: indices 9, 10, 11, 12 (6 DOFs each: u, d1, d2, d11, d22, d12)
        for (index_t c = 0; c < 4; ++c)
        {
            const index_t off = sumUntil(sz, 9 + c);
            scales[p](off + 1) = invH;
            scales[p](off + 2) = invH;
            scales[p](off + 3) = invH2;
            scales[p](off + 4) = invH2;
            scales[p](off + 5) = invH2;
        }

        argBasis[p].matrix = argBasis[p].matrix * scales[p].asDiagonal();
    }
    return scales;
}

/**
 * @brief Build a multi-component (e.g., nComp = 3 for Kirchhoff-Love shells)
 *        gsDofMapper for the conforming AS-G1 Argyris space [V_A^1]^{nComp}.
 *
 * Supports component-specific boundary conditions via bc[comp] (unknown index = comp)
 * using labels:
 *   - "ValuesAndDerivatives": strong value + crossing-derivative elimination (C^1 clamped per component)
 *   - "Values" / "Dirichlet" / "Clamped": strong value elimination on edge level-0
 *       (for "Clamped" on curved KL shells, physical normal rotation theta_nu = A_3 . d_nu u = 0
 *        is enforced by gsAsG1KLShellAssembler so in-plane membrane normal strains remain free)
 *   - "Derivatives": strong crossing-derivative elimination on edge level-1
 * plus per-component corner conditions ccPerComp[comp].
 */
template <typename T>
gsDofMapper makeVectorMapperForArgyrisBasis(
    const gsMultiPatch<T>& mp,
    const std::vector<gsArgyrisEmbedding<T>>& argBasis,
    index_t nComp,
    const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
    const std::vector<cornerConditions>& ccPerComp = std::vector<cornerConditions>())
{
    gsVector<index_t> patchDofSizes(argBasis.size());
    for (size_t i = 0; i < argBasis.size(); ++i)
        patchDofSizes[i] = argBasis[i].matrix.cols();

    gsDofMapper mapper(patchDofSizes, nComp);

    // 1. Match shared interface and corner DOFs for each component independently
    for (index_t comp = 0; comp < nComp; ++comp)
    {
        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            const boundaryInterface& ifc = *it;
            const patchSide ps1 = ifc.first();
            const patchSide ps2 = ifc.second();
            const gsVector<index_t, 13>& sz1 = argBasis[ps1.patch].sizes;
            const gsVector<index_t, 13>& sz2 = argBasis[ps2.patch].sizes;

            const index_t nLvl0 = sz1[1 + 2 * (ps1.m_index - 1)];
            const index_t offLvl0_1 = sumUntil(sz1, 1 + 2 * (ps1.m_index - 1));
            const index_t offLvl0_2 = sumUntil(sz2, 1 + 2 * (ps2.m_index - 1));

            const index_t nLvl1 = sz1[2 + 2 * (ps1.m_index - 1)];
            const index_t offLvl1_1 = sumUntil(sz1, 2 + 2 * (ps1.m_index - 1));
            const index_t offLvl1_2 = sumUntil(sz2, 2 + 2 * (ps2.m_index - 1));

            GISMO_ASSERT(nLvl0 == sz2[1 + 2 * (ps2.m_index - 1)], "Level-0 dimension mismatch.");
            GISMO_ASSERT(nLvl1 == sz2[2 + 2 * (ps2.m_index - 1)], "Level-1 dimension mismatch.");

            const short_t tanDir1 = 1 - ps1.direction();
            const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

            for (index_t j1 = 0; j1 < nLvl0; ++j1)
            {
                const index_t j2 = flipped ? (nLvl0 - 1 - j1) : j1;
                mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2, comp);
            }
            for (index_t j1 = 0; j1 < nLvl1; ++j1)
            {
                const index_t j2 = flipped ? (nLvl1 - 1 - j1) : j1;
                mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2, comp);
            }

            std::vector<patchCorner> corners;
            ps1.getContainedCorners(2, corners);
            GISMO_ASSERT(corners.size() == 2, "Unexpected number of corners.");
            for (index_t i = 0; i < 2; ++i)
            {
                const index_t c1 = corners[i].m_index - 1;
                const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;
                const index_t off_c1 = sumUntil(sz1, 9 + c1);
                const index_t off_c2 = sumUntil(sz2, 9 + c2);
                for (index_t j = 0; j < 6; ++j)
                    mapper.matchDof(ps1.patch, off_c1 + j, ps2.patch, off_c2 + j, comp);
            }
        }
    }

    // 2. Eliminate boundary edge DOFs per component
    auto eliminateEdgeBlock = [&](const patchSide& ps, index_t levelOffset, index_t comp) {
        if (comp < 0 || comp >= nComp) return;
        const index_t nLvl = argBasis[ps.patch].sizes[levelOffset + 2 * (ps.m_index - 1)];
        const index_t offLvl = sumUntil(argBasis[ps.patch].sizes, levelOffset + 2 * (ps.m_index - 1));
        for (index_t j = 0; j < nLvl; ++j)
            mapper.eliminateDof(offLvl + j, ps.patch, comp);
    };

    for (auto it = bc.begin("ValuesAndDerivatives"); it != bc.end("ValuesAndDerivatives"); ++it)
    {
        eliminateEdgeBlock(it->ps, 1, it->unknown());
        eliminateEdgeBlock(it->ps, 2, it->unknown());
    }
    for (auto it = bc.begin("Clamped"); it != bc.end("Clamped"); ++it)
        eliminateEdgeBlock(it->ps, 1, it->unknown());
    for (auto it = bc.begin("Values"); it != bc.end("Values"); ++it)
        eliminateEdgeBlock(it->ps, 1, it->unknown());
    for (auto it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
        eliminateEdgeBlock(it->ps, 1, it->unknown());
    for (auto it = bc.begin("Derivatives"); it != bc.end("Derivatives"); ++it)
        eliminateEdgeBlock(it->ps, 2, it->unknown());

    // 3. Eliminate boundary corner DOFs per component
    for (index_t comp = 0; comp < nComp; ++comp)
    {
        if (static_cast<size_t>(comp) >= ccPerComp.size())
            continue;
        const cornerConditions& cc = ccPerComp[comp];
        for (size_t i = 0; i < cc.size(); ++i)
        {
            const index_t c = cc[i].first.m_index - 1;
            const index_t p = cc[i].first.patch;
            const index_t off_corner = sumUntil(argBasis[p].sizes, 9 + c);
            eliminateCornerConditionDofs(mapper, off_corner, p, cc[i].second, comp);
        }
    }

    mapper.finalize();
    return mapper;
}

/**
 * @brief Build sparse global-to-discontinuous transfer matrices T_free and T_bnd
 *        for an nComp-component vector Argyris space.
 *
 * Row layout of T_free / T_bnd:
 *   row = comp * nDiscontScalar + patchRowOffset[patch] + localBsplineDof
 * Column layout:
 *   T_free: 0 .. mapper.freeSize() - 1
 *   T_bnd:  0 .. mapper.boundarySize() - 1
 */
template <typename T>
void buildVectorArgyrisTransferMatrices(
    const std::vector<gsArgyrisEmbedding<T>>& argBasis,
    const gsDofMapper& mapper,
    index_t nComp,
    gsSparseMatrix<T>& T_free,
    gsSparseMatrix<T>& T_bnd)
{
    const size_t nP = argBasis.size();
    std::vector<index_t> patchRowOff(nP + 1, 0);
    for (size_t p = 0; p < nP; ++p)
        patchRowOff[p + 1] = patchRowOff[p] + argBasis[p].matrix.rows();
    const index_t nDiscontScalar = patchRowOff[nP];
    const index_t nDiscontTotal = nComp * nDiscontScalar;

    gsSparseEntries<T> entFree, entBnd;
    index_t nnzEst = 0;
    for (size_t p = 0; p < nP; ++p)
        nnzEst += argBasis[p].matrix.nonZeros();
    entFree.reserve(nComp * nnzEst);
    entBnd.reserve(nComp * nnzEst / 4 + 16);

    for (index_t comp = 0; comp < nComp; ++comp)
    {
        const index_t compRowBase = comp * nDiscontScalar;
        for (size_t p = 0; p < nP; ++p)
        {
            const gsSparseMatrix<T>& E = argBasis[p].matrix;
            const index_t rowOff = compRowBase + patchRowOff[p];
            for (index_t k = 0; k < E.outerSize(); ++k)
            {
                for (typename gsSparseMatrix<T>::InnerIterator it(E, k); it; ++it)
                {
                    const index_t localCol = it.col();
                    const index_t row = rowOff + it.row();
                    const T val = it.value();
                    if (mapper.is_free(localCol, static_cast<index_t>(p), comp))
                    {
                        const index_t gIdx = mapper.index(localCol, static_cast<index_t>(p), comp);
                        entFree.add(row, gIdx, val);
                    }
                    else
                    {
                        const index_t bIdx = mapper.bindex(localCol, static_cast<index_t>(p), comp);
                        entBnd.add(row, bIdx, val);
                    }
                }
            }
        }
    }

    T_free.resize(nDiscontTotal, mapper.freeSize());
    T_free.setFrom(entFree);
    T_bnd.resize(nDiscontTotal, mapper.boundarySize());
    T_bnd.setFrom(entBnd);
}

/**
 * @brief Exact 6-jet of a smooth physical scalar function U(x,y,z) on a 3D surface
 *        patch corner in the shared tangent-plane frame [d1, d2, kappa_g].
 *
 * Uses the exact parametric chain rule matching collocateCorners(stride == 7),
 * automatically capturing the Weingarten shape-operator curvature contribution.
 */
template <typename T>
gsVector<T> evalSurfaceCornerJet3D(
    const gsGeometry<T>& geo,
    index_t cornerIdx,
    const gsMatrix<T>& cornerFramesRow,
    T val,
    const gsVector<T>& gradPhys,
    const gsMatrix<T>& hessPhys)
{
    const gsMatrix<T> supp = geo.support();
    gsMatrix<T> pt(2, 1);
    pt(0, 0) = supp(0, cornerIdx % 2);
    pt(1, 0) = supp(1, cornerIdx / 2);

    gsMatrix<T> gDeriv, gDeriv2;
    geo.deriv_into(pt, gDeriv);
    geo.deriv2_into(pt, gDeriv2);

    const index_t stride = cornerFramesRow.cols() / 4;
    GISMO_ENSURE(stride == 7 || stride == 10,
                 "evalSurfaceCornerJet3D expects cornerFramesRow with 28 or 40 columns, got "
                 << cornerFramesRow.cols());

    gsVector<T> d1(3), d2(3);
    d1 << cornerFramesRow(0, stride * cornerIdx + 0),
          cornerFramesRow(0, stride * cornerIdx + 1),
          cornerFramesRow(0, stride * cornerIdx + 2);
    d2 << cornerFramesRow(0, stride * cornerIdx + 3),
          cornerFramesRow(0, stride * cornerIdx + 4),
          cornerFramesRow(0, stride * cornerIdx + 5);
    const T cKappa = (stride == 7) ? cornerFramesRow(0, stride * cornerIdx + 6) : T(0);
    const T k0_e1  = (stride == 10) ? cornerFramesRow(0, stride * cornerIdx + 6) : T(0);
    const T k0_e2  = (stride == 10) ? cornerFramesRow(0, stride * cornerIdx + 7) : T(0);
    const T k1_e1  = (stride == 10) ? cornerFramesRow(0, stride * cornerIdx + 8) : T(0);
    const T k1_e2  = (stride == 10) ? cornerFramesRow(0, stride * cornerIdx + 9) : T(0);

    gsVector<T> e1 = d1.normalized();
    gsVector<T> e2 = d2 - d2.dot(e1) * e1;
    e2.normalize();

    gsMatrix<T> E(3, 2);
    E.col(0) = e1; E.col(1) = e2;
    gsMatrix<T> D(3, 2);
    D.col(0) = d1; D.col(1) = d2;
    const gsMatrix<T> Tmat = E.transpose() * D;

    gsMatrix<T> J(3, 2);
    gsVector<T> A1(3), A2(3), Auu(3), Avv(3), Auv(3);
    for (short_t c = 0; c < 3; ++c)
    {
        A1(c)  = gDeriv(2 * c + 0, 0);
        A2(c)  = gDeriv(2 * c + 1, 0);
        J(c, 0) = A1(c);
        J(c, 1) = A2(c);
        Auu(c) = gDeriv2(3 * c + 0, 0);
        Avv(c) = gDeriv2(3 * c + 1, 0);
        Auv(c) = gDeriv2(3 * c + 2, 0);
    }

    const gsMatrix<T> Jhat = E.transpose() * J;
    const gsMatrix<T> JhatInv = Jhat.inverse();
    const gsMatrix<T> JhatInvT = JhatInv.transpose();

    gsMatrix<T> Hy0(2, 2), Hy1(2, 2);
    Hy0(0, 0) = e1.dot(Auu); Hy0(0, 1) = e1.dot(Auv);
    Hy0(1, 0) = Hy0(0, 1);   Hy0(1, 1) = e1.dot(Avv);
    Hy1(0, 0) = e2.dot(Auu); Hy1(0, 1) = e2.dot(Auv);
    Hy1(1, 0) = Hy1(0, 1);   Hy1(1, 1) = e2.dot(Avv);

    // Parametric gradient and Hessian of u(xi) = U(G(xi))
    gsVector<T> gXi(2);
    gXi(0) = gradPhys.dot(A1);
    gXi(1) = gradPhys.dot(A2);

    gsMatrix<T> HXi(2, 2);
    HXi(0, 0) = A1.dot(hessPhys * A1) + gradPhys.dot(Auu);
    HXi(1, 1) = A2.dot(hessPhys * A2) + gradPhys.dot(Avv);
    HXi(0, 1) = A1.dot(hessPhys * A2) + gradPhys.dot(Auv);
    HXi(1, 0) = HXi(0, 1);

    const gsVector<T> gy = JhatInvT * gXi;
    const gsMatrix<T> Hy = JhatInvT * (HXi - gy(0) * Hy0 - gy(1) * Hy1) * JhatInv;

    const gsVector<T> gd = Tmat.transpose() * gy;
    const gsMatrix<T> Hd = Tmat.transpose() * Hy * Tmat;

    gsVector<T> jet(6);
    jet(0) = val;
    jet(1) = gd(0);
    jet(2) = gd(1);
    if (stride == 10)
    {
        jet(3) = Hd(0, 0) + k0_e1 * gy(0) + k0_e2 * gy(1);
        jet(4) = Hd(1, 1) + k1_e1 * gy(0) + k1_e2 * gy(1);
    }
    else
    {
        jet(3) = Hd(0, 0);
        jet(4) = Hd(1, 1) + cKappa * gd(0);
    }
    jet(5) = Hd(0, 1);
    return jet;
}

/**
 * @brief Compute patch-local AS-G1 degrees of freedom via the Hasanova, S. Takacs,
 *        T. Takacs (2026) Section 3.4 projector \Pi_{AS} (Eq. (41)-(43), (61)-(63)):
 *          - Exact 6-jets at all 4 patch corners (C^2 vertex interpolation),
 *          - 1D H^2-Ritz trace projection \Pi_j^{(0)} in S_{p, k+1, Z}^{inner},
 *          - 1D H^1-Ritz crossing-derivative projection \Pi_j^{(1)} in S_{p-1, k, Z}^{inner},
 *          - (Optional) 2D patch-interior H^2-Ritz projection Q in S_{p, k, Z}^0.
 *
 * ExactFn signature:
 *   void evalExact(const gsVector<T>& xPhys, T& val, gsVector<T>& gradPhys, gsMatrix<T>& hessPhys)
 *
 * If includeInterior == false, only corner and edge DOFs are populated (sufficient
 * for exact p-robust boundary lifting u_bnd!).
 */
template <typename T, typename ExactFn>
gsVector<T> computePatchAsG1ProjectionDofs3D(
    const gsTensorBSplineBasis<2, T>& tb,
    const gsGeometry<T>& geo,
    const gsMatrix<T>& localGluingRow,
    const gsMatrix<T>& cornerFramesRow,
    const gsArgyrisEmbedding<T>& unscaledEmb,
    const gsVector<T>& colScale,
    const ExactFn& evalExact,
    bool includeInterior = true)
{
    const gsVector<index_t, 13>& sz = unscaledEmb.sizes;
    const index_t nLoc = unscaledEmb.matrix.cols();
    gsVector<T> dUnscaled = gsVector<T>::Zero(nLoc);

    const gsMatrix<T> supp = tb.support();

    // 1. Exact 6-jets at all 4 corners
    for (index_t c = 0; c < 4; ++c)
    {
        gsMatrix<T> pt(2, 1);
        pt(0, 0) = supp(0, c % 2);
        pt(1, 0) = supp(1, c / 2);
        gsMatrix<T> xMat = geo.eval(pt);
        gsVector<T> xPhys = xMat.col(0);

        T val = T(0);
        gsVector<T> gP(3);
        gsMatrix<T> hP(3, 3);
        evalExact(xPhys, val, gP, hP);

        const gsVector<T> jet = evalSurfaceCornerJet3D(geo, c, cornerFramesRow, val, gP, hP);
        const index_t off = sumUntil(sz, 9 + c);
        dUnscaled.segment(off, 6) = jet;
    }

    // Reconstruct the corner-only B-spline control points c_corner = E_corners * d_corners
    const index_t offC0 = sumUntil(sz, 9);
    const gsVector<T> cCorner = unscaledEmb.matrix.middleCols(offC0, 24) * dUnscaled.segment(offC0, 24);

    // Helper: evaluate a 2D B-spline with coefficients c at parameter points pts(2, N),
    // returning val(N), du(N), dv(N), duu(N), dvv(N), duv(N).
    auto evalSpline2D = [&](const gsVector<T>& cVec, const gsMatrix<T>& pts,
                            gsVector<T>& v, gsVector<T>& du, gsVector<T>& dv,
                            gsVector<T>& duu, gsVector<T>& dvv, gsVector<T>& duv)
    {
        const index_t N = pts.cols();
        v.setZero(N); du.setZero(N); dv.setZero(N);
        duu.setZero(N); dvv.setZero(N); duv.setZero(N);
        std::vector<gsMatrix<T>> ders;
        tb.evalAllDers_into(pts, 2, ders);
        gsMatrix<index_t> act = tb.active(pts);
        for (index_t q = 0; q < N; ++q)
        {
            for (index_t a = 0; a < act.rows(); ++a)
            {
                const T ca = cVec(act(a, q));
                v(q)   += ca * ders[0](a, q);
                du(q)  += ca * ders[1](2 * a + 0, q);
                dv(q)  += ca * ders[1](2 * a + 1, q);
                duu(q) += ca * ders[2](3 * a + 0, q);
                dvv(q) += ca * ders[2](3 * a + 1, q);
                duv(q) += ca * ders[2](3 * a + 2, q);
            }
        }
    };

    // 2. For each of the 4 sides: 1D H^2-Ritz trace projection & 1D H^1-Ritz crossing-deriv projection
    for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
    {
        const index_t sIdx = side.m_index - 1;
        const short_t nDir = side.direction();
        const short_t tDir = 1 - nDir;
        const bool par = side.parameter();
        const T nOutSign = par ? T(1) : T(-1); // outward parametric normal sign (matching evaluateBasisAtBoundary)
        const T sgn = (par == (nDir != 0)) ? T(1) : T(-1); // matching collocateBoundaryCrossingDerivative

        const gsBSplineBasis<T> sideBasis = *tb.boundaryBasis(side);
        gsBSplineBasis<T> b0 = sideBasis; b0.elevateContinuity(1); // S_{p, k+1, Z}
        gsBSplineBasis<T> b1 = sideBasis; b1.degreeReduce(1);      // S_{p-1, k, Z}

        const index_t n0 = b0.size();
        const index_t n1 = b1.size();
        const index_t nIn0 = n0 - 6;
        const index_t nIn1 = n1 - 4;

        const T ta = sideBasis.support()(0, 0);
        const T tbEnd = sideBasis.support()(0, 1);
        const T invLen = T(1) / (tbEnd - ta);
        const T a0 = localGluingRow(0, 4 * sIdx + 0);
        const T a1 = localGluingRow(0, 4 * sIdx + 1);
        const T b0g = localGluingRow(0, 4 * sIdx + 2);
        const T b1g = localGluingRow(0, 4 * sIdx + 3);
        const T dAlphaDt = (a1 - a0) * invLen;
        const T dBetaDt  = (b1g - b0g) * invLen;

        gsGaussRule<T> rule(sideBasis.degree() + 3);
        gsMatrix<T> qNodes;
        gsVector<T> qWeights;
        rule.mapToAll(sideBasis.knots().breaks(), qNodes, qWeights);
        const index_t nQ = qNodes.cols();

        gsMatrix<T> pts2D(2, nQ);
        pts2D.row(nDir).setConstant(supp(nDir, par ? 1 : 0));
        pts2D.row(tDir) = qNodes.row(0);

        gsMatrix<T> Xq = geo.eval(pts2D);
        gsMatrix<T> gD1 = geo.deriv(pts2D);
        gsMatrix<T> gD2 = geo.deriv2(pts2D);

        gsVector<T> cV, cDu, cDv, cDuu, cDvv, cDuv;
        evalSpline2D(cCorner, pts2D, cV, cDu, cDv, cDuu, cDvv, cDuv);

        gsVector<T> rVal(nQ), rDt(nQ), rDtt(nQ), gVal(nQ), gDt(nQ);
        for (index_t q = 0; q < nQ; ++q)
        {
            gsVector<T> xPhys = Xq.col(q);
            T uVal = T(0);
            gsVector<T> gP(3);
            gsMatrix<T> hP(3, 3);
            evalExact(xPhys, uVal, gP, hP);

            gsVector<T> A1(3), A2(3), Auu(3), Avv(3), Auv(3);
            for (short_t c = 0; c < 3; ++c)
            {
                A1(c)  = gD1(2 * c + 0, q);
                A2(c)  = gD1(2 * c + 1, q);
                Auu(c) = gD2(3 * c + 0, q);
                Avv(c) = gD2(3 * c + 1, q);
                Auv(c) = gD2(3 * c + 2, q);
            }
            const T uDu  = gP.dot(A1);
            const T uDv  = gP.dot(A2);
            const T uDuu = A1.dot(hP * A1) + gP.dot(Auu);
            const T uDvv = A2.dot(hP * A2) + gP.dot(Avv);
            const T uDuv = A1.dot(hP * A2) + gP.dot(Auv);

            // Residual r = u_ex - u_corner and its parametric derivatives
            const T rv   = uVal - cV(q);
            const T rdu  = uDu  - cDu(q);
            const T rdv  = uDv  - cDv(q);
            const T rduu = uDuu - cDuu(q);
            const T rdvv = uDvv - cDvv(q);
            const T rduv = uDuv - cDuv(q);

            const T rt   = (tDir == 0) ? rdu  : rdv;
            const T rn   = nOutSign * ((nDir == 0) ? rdu : rdv);
            const T rtt  = (tDir == 0) ? rduu : rdvv;
            const T rnt  = nOutSign * rduv;

            rVal(q) = rv;
            rDt(q)  = rt;
            rDtt(q) = rtt;

            const T sNorm = (qNodes(0, q) - ta) * invLen;
            const T alpha = a0 * (T(1) - sNorm) + a1 * sNorm;
            const T beta  = b0g * (T(1) - sNorm) + b1g * sNorm;
            const T num   = rn + sgn * beta * rt;
            const T dNumDt = rnt + sgn * (dBetaDt * rt + beta * rtt);

            gVal(q) = num / alpha;
            gDt(q)  = (dNumDt * alpha - num * dAlphaDt) / (alpha * alpha);
        }

        // 2a. Level-0 trace projection in S_{p, k+1, Z}^{inner} (indices 3 .. n0-4)
        if (nIn0 > 0)
        {
            std::vector<gsMatrix<T>> ders0;
            b0.evalAllDers_into(qNodes, 2, ders0);
            gsMatrix<index_t> act0 = b0.active(qNodes);

            gsMatrix<T> A0 = gsMatrix<T>::Zero(nIn0, nIn0);
            gsVector<T> f0 = gsVector<T>::Zero(nIn0);
            for (index_t q = 0; q < nQ; ++q)
            {
                const T w = qWeights(q);
                for (index_t a = 0; a < act0.rows(); ++a)
                {
                    const index_t ia = act0(a, q) - 3;
                    if (ia < 0 || ia >= nIn0) continue;
                    const T va = ders0[0](a, q), da = ders0[1](a, q), dda = ders0[2](a, q);
                    f0(ia) += w * (va * rVal(q) + da * rDt(q) + dda * rDtt(q));
                    for (index_t b = 0; b < act0.rows(); ++b)
                    {
                        const index_t ib = act0(b, q) - 3;
                        if (ib < 0 || ib >= nIn0) continue;
                        A0(ia, ib) += w * (va * ders0[0](b, q) + da * ders0[1](b, q) + dda * ders0[2](b, q));
                    }
                }
            }
            const gsVector<T> sol0 = A0.ldlt().solve(f0);
            const index_t off0 = sumUntil(sz, 1 + 2 * sIdx);
            dUnscaled.segment(off0, nIn0) = sol0;
        }

        // 2b. Level-1 crossing-derivative projection in S_{p-1, k, Z}^{inner} (indices 2 .. n1-3)
        if (nIn1 > 0)
        {
            std::vector<gsMatrix<T>> ders1;
            b1.evalAllDers_into(qNodes, 1, ders1);
            gsMatrix<index_t> act1 = b1.active(qNodes);

            gsMatrix<T> A1 = gsMatrix<T>::Zero(nIn1, nIn1);
            gsVector<T> f1 = gsVector<T>::Zero(nIn1);
            for (index_t q = 0; q < nQ; ++q)
            {
                const T w = qWeights(q);
                for (index_t a = 0; a < act1.rows(); ++a)
                {
                    const index_t ia = act1(a, q) - 2;
                    if (ia < 0 || ia >= nIn1) continue;
                    const T va = ders1[0](a, q), da = ders1[1](a, q);
                    f1(ia) += w * (va * gVal(q) + da * gDt(q));
                    for (index_t b = 0; b < act1.rows(); ++b)
                    {
                        const index_t ib = act1(b, q) - 2;
                        if (ib < 0 || ib >= nIn1) continue;
                        A1(ia, ib) += w * (va * ders1[0](b, q) + da * ders1[1](b, q));
                    }
                }
            }
            const gsVector<T> sol1 = A1.ldlt().solve(f1);
            const index_t off1 = sumUntil(sz, 2 + 2 * sIdx);
            dUnscaled.segment(off1, nIn1) = sol1;
        }
    }

    // 3. Optional patch-interior H^2-Ritz projection Q (Paper Eq. (41)-(43))
    if (includeInterior && sz[0] > 0)
    {
        const index_t nInt = sz[0];
        const gsVector<T> cSkel = unscaledEmb.matrix.rightCols(nLoc - nInt) * dUnscaled.tail(nLoc - nInt);

        const gsBSplineBasis<T>& bu = tb.component(0);
        const gsBSplineBasis<T>& bv = tb.component(1);
        gsGaussRule<T> ruleU(bu.degree() + 2), ruleV(bv.degree() + 2);
        gsMatrix<T> qu, qv;
        gsVector<T> wu, wv;
        ruleU.mapToAll(bu.knots().breaks(), qu, wu);
        ruleV.mapToAll(bv.knots().breaks(), qv, wv);
        const index_t nQu = qu.cols(), nQv = qv.cols();
        const index_t nQ2 = nQu * nQv;

        gsMatrix<T> pts2D(2, nQ2);
        gsVector<T> w2D(nQ2);
        for (index_t j = 0; j < nQv; ++j)
            for (index_t i = 0; i < nQu; ++i)
            {
                const index_t q = i + j * nQu;
                pts2D(0, q) = qu(0, i);
                pts2D(1, q) = qv(0, j);
                w2D(q) = wu(i) * wv(j);
            }

        gsMatrix<T> Xq = geo.eval(pts2D);
        gsMatrix<T> gD1 = geo.deriv(pts2D);
        gsMatrix<T> gD2 = geo.deriv2(pts2D);
        gsVector<T> sV, sDu, sDv, sDuu, sDvv, sDuv;
        evalSpline2D(cSkel, pts2D, sV, sDu, sDv, sDuu, sDvv, sDuv);

        // Map interior B-spline index -> local interior column (0 .. nInt-1)
        const std::vector<index_t> intIds = getInteriorDofs(tb);
        std::vector<index_t> bspToInt(tb.size(), -1);
        for (index_t k = 0; k < nInt; ++k)
            bspToInt[intIds[k]] = k;

        std::vector<gsMatrix<T>> ders;
        tb.evalAllDers_into(pts2D, 2, ders);
        gsMatrix<index_t> act = tb.active(pts2D);

        gsSparseEntries<T> entries;
        entries.reserve(nQ2 * 16);
        gsVector<T> rhsInt = gsVector<T>::Zero(nInt);

        for (index_t q = 0; q < nQ2; ++q)
        {
            const T w = w2D(q);
            gsVector<T> xPhys = Xq.col(q);
            T uVal = T(0);
            gsVector<T> gP(3);
            gsMatrix<T> hP(3, 3);
            evalExact(xPhys, uVal, gP, hP);

            gsVector<T> A1(3), A2(3), Auu(3), Avv(3), Auv(3);
            for (short_t c = 0; c < 3; ++c)
            {
                A1(c)  = gD1(2 * c + 0, q);
                A2(c)  = gD1(2 * c + 1, q);
                Auu(c) = gD2(3 * c + 0, q);
                Avv(c) = gD2(3 * c + 1, q);
                Auv(c) = gD2(3 * c + 2, q);
            }
            const T rv   = uVal - sV(q);
            const T rdu  = gP.dot(A1) - sDu(q);
            const T rdv  = gP.dot(A2) - sDv(q);
            const T rduu = (A1.dot(hP * A1) + gP.dot(Auu)) - sDuu(q);
            const T rdvv = (A2.dot(hP * A2) + gP.dot(Avv)) - sDvv(q);
            const T rduv = (A1.dot(hP * A2) + gP.dot(Auv)) - sDuv(q);

            for (index_t a = 0; a < act.rows(); ++a)
            {
                const index_t ia = bspToInt[act(a, q)];
                if (ia < 0) continue;
                const T va   = ders[0](a, q);
                const T dua  = ders[1](2 * a + 0, q), dva  = ders[1](2 * a + 1, q);
                const T duua = ders[2](3 * a + 0, q), dvva = ders[2](3 * a + 1, q), duva = ders[2](3 * a + 2, q);

                rhsInt(ia) += w * (va * rv + dua * rdu + dva * rdv
                                   + duua * rduu + dvva * rdvv + T(2) * duva * rduv);
                for (index_t b = 0; b < act.rows(); ++b)
                {
                    const index_t ib = bspToInt[act(b, q)];
                    if (ib < 0) continue;
                    const T vb   = ders[0](b, q);
                    const T dub  = ders[1](2 * b + 0, q), dvb  = ders[1](2 * b + 1, q);
                    const T duub = ders[2](3 * b + 0, q), dvvb = ders[2](3 * b + 1, q), duvb = ders[2](3 * b + 2, q);
                    const T val = w * (va * vb + dua * dub + dva * dvb
                                       + duua * duub + dvva * dvvb + T(2) * duva * duvb);
                    entries.add(ia, ib, val);
                }
            }
        }

        gsSparseMatrix<T> Aint(nInt, nInt);
        Aint.setFrom(entries);
        gsEigen::SimplicialLDLT<gsSparseMatrix<T>> ldlt;
        ldlt.compute(Aint);
        dUnscaled.head(nInt) = ldlt.solve(rhsInt);
    }

    // Convert unscaled DOFs to scaled DOFs (since E_scaled = E_unscaled * diag(colScale))
    gsVector<T> dScaled(nLoc);
    for (index_t k = 0; k < nLoc; ++k)
        dScaled(k) = (colScale.size() == nLoc && colScale(k) != T(0))
                   ? (dUnscaled(k) / colScale(k))
                   : dUnscaled(k);
    return dScaled;
}

/**
 * @brief Populate the eliminated boundary DOF vector u_bnd (size mapper.boundarySize())
 *        and optionally the full standalone paper quasi-interpolant u_qi_free (size mapper.freeSize())
 *        for a 3-component vector displacement field on a 3D AS-G1 surface.
 *
 * ExactVecFn signature:
 *   void evalExactComp(index_t comp, const gsVector<T>& xPhys,
 *                      T& val, gsVector<T>& gradPhys, gsMatrix<T>& hessPhys)
 */
template <typename T, typename ExactVecFn>
void projectVectorAsG1Exact3D(
    const gsMultiPatch<T>& mp,
    const gsMultiBasis<T>& dbasis,
    const gsMatrix<T>& gluingData,
    const gsMatrix<T>& cornerFrames,
    const std::vector<gsArgyrisEmbedding<T>>& unscaledArgBasis,
    const std::vector<gsVector<T>>& colScales,
    const gsDofMapper& mapper,
    index_t nComp,
    const ExactVecFn& evalExactComp,
    gsVector<T>& u_bnd,
    gsVector<T>* u_free_qi = nullptr)
{
    const index_t nBnd = mapper.boundarySize();
    const index_t nFree = mapper.freeSize();
    u_bnd.setZero(nBnd);
    gsVector<index_t> bndCount = gsVector<index_t>::Zero(nBnd);

    const bool doInterior = (u_free_qi != nullptr);
    gsVector<index_t> freeCount;
    if (doInterior)
    {
        u_free_qi->setZero(nFree);
        freeCount.setZero(nFree);
    }

    const size_t nP = mp.nPatches();
    for (index_t comp = 0; comp < nComp; ++comp)
    {
        auto evalComp = [&](const gsVector<T>& xP, T& v, gsVector<T>& gP, gsMatrix<T>& hP) {
            evalExactComp(comp, xP, v, gP, hP);
        };

        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(dbasis[p]);
            const gsVector<T> dLoc = computePatchAsG1ProjectionDofs3D(
                tb, mp.patch(p), gsMatrix<T>(gluingData.row(p)), gsMatrix<T>(cornerFrames.row(p)),
                unscaledArgBasis[p], colScales[p], evalComp, doInterior);

            for (index_t k = 0; k < dLoc.size(); ++k)
            {
                if (mapper.is_boundary(k, static_cast<index_t>(p), comp))
                {
                    const index_t bIdx = mapper.bindex(k, static_cast<index_t>(p), comp);
                    u_bnd(bIdx) += dLoc(k);
                    bndCount(bIdx) += 1;
                }
                else if (doInterior)
                {
                    const index_t fIdx = mapper.index(k, static_cast<index_t>(p), comp);
                    (*u_free_qi)(fIdx) += dLoc(k);
                    freeCount(fIdx) += 1;
                }
            }
        }
    }

    for (index_t i = 0; i < nBnd; ++i)
        if (bndCount(i) > 1) u_bnd(i) /= T(bndCount(i));
    if (doInterior)
        for (index_t i = 0; i < nFree; ++i)
            if (freeCount(i) > 1) (*u_free_qi)(i) /= T(freeCount(i));
}

} // namespace gismo
