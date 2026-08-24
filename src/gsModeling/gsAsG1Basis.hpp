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

    // Transform first and second derivatives to physical domain using map data
    gsMapData<T> md(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER);
    md.points = corners;
    geo.computeMap(md);

    // Reserve matrix for result
    gsMatrix<T> result(4*6, tensorBasis.size());
    result.setZero();

    // For each corner
    for (index_t i=0; i<4; ++i)
    {
        // Transform gradient to phyiscal domain
        gsMatrix<T> der1Phys;
        transformGradients(md, i, der1, der1Phys);
        // Transform x and y derivatives into derivatives in normal and tangenial direction
        gsMatrix<T> transform(2,2);
        transform(0,0) = normals(0+2*i);
        transform(1,0) = normals(1+2*i);
        transform(0,1) = normals(1+2*i);
        transform(1,1) = -normals(0+2*i);
        der1Phys = transform * der1Phys;
        // Transform into uniform shape
        der1Phys.resize(der1.rows(), 1);

        // Transform Hessian to physical domain
        gsMatrix<T> der2Phys;
        transformDeriv2Hgrad(md, i, der1, der2, der2Phys);
        // Transform x and y derivatives into derivatives in normal and tangenial direction
        for (index_t j=0; j<idx.rows(); ++j)
            der2Phys.row(j) = hessMatrixAsReduced<T>(
                    transform * reducedHessAsMatrix<T>(der2Phys.row(j)) * transform
            );
        // Transform into uniform shape
        der2Phys = der2Phys.transpose();
        der2Phys.resize(der2.rows(), 1);

        for (index_t j=0; j<idx.rows(); ++j)
        {
            result(6*i+0, idx(j, i)) = vals(j, i);            //      u
            result(6*i+1, idx(j, i)) = der1Phys(2 * j);       // d_n  u
            result(6*i+2, idx(j, i)) = der1Phys(2 * j + 1);   // d_t  u
            result(6*i+3, idx(j, i)) = der2Phys(3 * j);       // d_nn u
            result(6*i+4, idx(j, i)) = der2Phys(3 * j + 1);   // d_tt u
            result(6*i+5, idx(j, i)) = der2Phys(3 * j + 2);   // d_nt u
        }
    }
    return result.sparseView(1, 1e-4);
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
    T eps = 1e-4)
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

    return result.topRows(rows).sparseView(1e-4);
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
    all = 6,
    valuesNormals = 5,
    none = 0
};

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

    // Now, incorporate the boundary conditons
    // TODO: Implement separate conditions for "Values" and "Derivatives"
    for (auto it = bc.begin("ValuesAndDerivatives"); it != bc.end("ValuesAndDerivatives"); ++it)
    {
        const patchSide &ps = it->ps;

        // External boundary side
        const index_t nLvl0 = argBasis[ps.patch].sizes[1 + 2 * (ps.m_index - 1)];
        const index_t offLvl0 = sumUntil(argBasis[ps.patch].sizes, 1 + 2 * (ps.m_index - 1));
        for (index_t j = 0; j < nLvl0; ++j)
          mapper.eliminateDof(offLvl0 + j, ps.patch);

        const index_t nLvl1 = argBasis[ps.patch].sizes[2 + 2 * (ps.m_index - 1)];
        const index_t offLvl1 =
            sumUntil(argBasis[ps.patch].sizes, 2 + 2 * (ps.m_index - 1));
        for (index_t j = 0; j < nLvl1; ++j)
          mapper.eliminateDof(offLvl1 + j, ps.patch);
    }

    for (size_t i=0; i<cc.size(); ++i)
    {
        const index_t c = cc[i].first.m_index - 1;
        const index_t off_corner = sumUntil(argBasis[cc[i].first.patch].sizes, 9 + c);
        if (cc[i].second == cornerConditionType::all)
        {
            for (index_t j=0; j<6; ++j)
                mapper.eliminateDof(off_corner+j, cc[i].first.patch);
        }
        else if (cc[i].second == cornerConditionType::valuesNormals)
        {
            for (index_t j=0; j<6; ++j)
                if (j!=3)
                    mapper.eliminateDof(off_corner+j, cc[i].first.patch);
        }
        else if (cc[i].second == cornerConditionType::none)
        {
            // Do nothing
        }
        else
        {
            GISMO_ENSURE(false, "Not implemented.");
        }
    }

    mapper.finalize();
    return mapper;
}


} // namespace gismo
