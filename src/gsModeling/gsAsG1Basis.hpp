/** @file gsAsG1Basis.hpp

    @brief Derive a C1 basis for an AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

namespace gismo {

// ====================================================================
// Helper utilities (matrix manipulations algebra)
// ====================================================================

template<typename T>
gsSparseMatrix<T> asEmbeddingMatrix(index_t rows, const gsMatrix<index_t>& c)
{
    index_t cols = c.rows();
    GISMO_ASSERT(c.cols() == 1, "");
    gsSparseEntries<T> se;
    se.reserve(cols);
    for (index_t i=0; i<c.size(); ++i)
        se.add(c(i,0),i,1);
    gsSparseMatrix<T>result(rows, cols);
    result.setFrom(se);
    return result;
}


// ====================================================================
// Helper utilities (basis evaluation)
// ====================================================================

template<typename T>
bool isNested(const gsBSplineBasis<T>& coarseBasis, const gsBSplineBasis<T>& fineBasis)
{
    if (coarseBasis.dim() != 1 || fineBasis.dim() != 1)
        return false;
    const int coarseDegree = coarseBasis.degree();
    const int fineDegree   = fineBasis.degree();
    if (coarseDegree > fineDegree)
        return false;
    const gsKnotVector<T>& coarseKnots = coarseBasis.knots();
    const gsKnotVector<T>& fineKnots   = fineBasis.knots();
    gsKnotVector<T> elevated = coarseKnots;
    elevated.degreeElevate(fineDegree - coarseDegree);
    return fineKnots.includes(elevated);
}

template<typename T>
gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>& basis,
                                    const gsMatrix<T>& pts,
                                    const index_t derivative = 0)
{
    gsMatrix<T> vals;
    if (derivative == 0)
        vals = basis.eval(pts);
    else if (derivative == 1)
        vals = basis.deriv(pts);
    else
        GISMO_ERROR ("Invalid value.");
    gsMatrix<index_t> idx  = basis.active(pts);
    gsSparseEntries<T> entries;
    entries.reserve(vals.rows() * vals.cols());
    for (int i = 0; i < vals.cols(); ++i)
        for (int j = 0; j < vals.rows(); ++j)
            entries.add(i, idx(j, i), vals(j, i));
    gsSparseMatrix<T> C(pts.cols(), basis.size());
    C.setFrom(entries);
    return C;
}

template<typename T>
gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarse,
                                  const gsBSplineBasis<T>& fine)
{
    gsMatrix<T> greville = fine.anchors();
    gsSparseMatrix<T> Cc = collocationMatrix(coarse, greville);
    gsSparseMatrix<T> Cf = collocationMatrix(fine,   greville);
    gsMatrix<T> result;
    makeSparseLUSolver(Cf)->apply(Cc, result);
    return result.sparseView(1, 1e-10);
}

template<typename T>
gsSparseMatrix<T> evaluateBasisAtBoundary(
        const gsBSplineBasis<T>& basis,
        index_t sideParameter,
        index_t derivative,
        index_t layers
)
{
    GISMO_ASSERT (sideParameter==0 || sideParameter==1, "Invalid side.");

    const gsMatrix<T> pt = basis.support().col(sideParameter);
    gsMatrix<index_t> indices = basis.active(pt);

    gsMatrix<T> values;
    if (derivative==0)
        values = basis.eval(pt);
    else if (derivative==1)
    {
        values = basis.deriv(pt);
        // On left hand side, normal vector is -1
        if (sideParameter==0)
            values *= -1;
    }
    else
        GISMO_ERROR ("Invalid value.");

    // Flip indices if right side is selected
    if (sideParameter==1)
        indices.array() = basis.size() - 1 - indices.array();

    // Store in sparse matrix
    gsSparseMatrix<T> result(1, layers);
    for (index_t i=0; i<indices.rows(); ++i)
        if (indices(i,0)<layers)
            result(0, indices(i,0)) = values(i,0);


    return result;
}


// ====================================================================
// Evaluation of tensor basis
// ====================================================================

template<typename T>
gsSparseMatrix<T> collocateBoundaryValues(
        const gsTensorBSplineBasis<2,T>& tensorBasis,
        boxSide side,
        const gsMatrix<T>& pts
)
{
    const gsBSplineBasis<T> normalBasis     =  tensorBasis.component(side.direction());
    const gsBSplineBasis<T> tangentialBasis = *tensorBasis.boundaryBasis(side);

    return evaluateBasisAtBoundary(normalBasis, side.parameter(), 0, 2).kron(collocationMatrix(tangentialBasis, pts));
}

template<typename T>
gsSparseMatrix<T> collocateBoundaryCrossingDerivative(
        const gsTensorBSplineBasis<2,T>& tensorBasis,
        boxSide side,
        const gsMatrix<T>& gluingData,
        const gsMatrix<T>& pts
)
{
    const gsBSplineBasis<T> normalBasis      = tensorBasis.component(side.direction());
    const gsBSplineBasis<T> tangentialBasis = *tensorBasis.boundaryBasis(side);

    const gsMatrix<T> tangentialBasisSupport = tangentialBasis.support();
    const T a = tangentialBasisSupport(0,0);
    const T b = tangentialBasisSupport(0,1);
    const gsMatrix<T> ptsNormalized = (pts.transpose().array()-a)/(b-a);
    const gsMatrix<T> alpha = gluingData(0,0)*(1-ptsNormalized.array()) + gluingData(0,1)*ptsNormalized.array();
    const gsMatrix<T> beta  = gluingData(0,2)*(1-ptsNormalized.array()) + gluingData(0,3)*ptsNormalized.array();
    const gsMatrix<T> scaling0 = beta.array()/alpha.array() * (side.parameter()==side.direction() ? 1 : -1); //This is as in the paper
    const gsMatrix<T> scaling1 = 1/alpha.array();

    return evaluateBasisAtBoundary(normalBasis, side.parameter(), 0, 2)
              .kron(scaling0.asDiagonal() * collocationMatrix(tangentialBasis, pts, 1))
            + evaluateBasisAtBoundary(normalBasis, side.parameter(), 1, 2)
              .kron(scaling1.asDiagonal() * collocationMatrix(tangentialBasis, pts, 0));
}


// ====================================================================
// Derive various embeddings
// ====================================================================

template<typename T>
struct gsEdgeEmbedding {
    gsSparseMatrix<T> matrix;
    gsVector<index_t,2> sizes;
};

template<typename T>
gsEdgeEmbedding<T> deriveEdgeEmbedding(
        const gsTensorBSplineBasis<2,T>& tensorBasis,
        const gsMatrix<T>& localGluingData,
        boxSide side,
        T eps = 1e-4
)
{
    const gsBSplineBasis<T> sideBasis = *tensorBasis.boundaryBasis(side);
    const gsMatrix<T> greville = sideBasis.anchors();

    const gsSparseMatrix<T> collocateLeft = gsBlockSparseMatrix<T>(2,1)
        .set(0,0,collocateBoundaryValues(tensorBasis, side, greville))
        .set(1,0,collocateBoundaryCrossingDerivative(tensorBasis, side, localGluingData, greville));

    gsBSplineBasis<T> sideLowerDegreeBasis = sideBasis;
    sideLowerDegreeBasis.degreeReduce(1);

    gsBSplineBasis<T> sideSmootherBasis = sideBasis;
    sideSmootherBasis.elevateContinuity(1);

    const gsSparseMatrix<T> collocateRight = gsBlockSparseMatrix<T>(2,2)
            .set(0,0,collocationMatrix(sideSmootherBasis, greville))
            .set(1,1,collocationMatrix(sideLowerDegreeBasis, greville));

    const gsLinearOperator<>::Ptr solver = makeSparseLUSolver(collocateLeft);
    gsMatrix<T> result;
    solver->apply(collocateRight, result);

    gsEdgeEmbedding<T> ee;
    ee.matrix = result.sparseView(1, eps);
    ee.sizes[0] = sideSmootherBasis.size();
    ee.sizes[1] = sideLowerDegreeBasis.size();
    return ee;
}

template<typename T>
std::vector<index_t> getInteriorDofs(const gsTensorBSplineBasis<2,T>& tensorBasis, boxSide side)
{
    std::vector<index_t> removedSet;
    for (index_t i=0; i<2; ++i)
    {
        gsMatrix<index_t> bdyDofs = tensorBasis.boundaryOffset(side, i);
        std::copy(bdyDofs.data(), bdyDofs.data()+bdyDofs.rows(), std::back_inserter(removedSet));
    }

    std::sort(removedSet.begin(), removedSet.end());

    std::vector<index_t> all(tensorBasis.size());
    std::iota(all.begin(), all.end(), 0);

    std::vector<index_t> interior;
    interior.reserve(all.size());
    std::set_difference(all.begin(), all.end(),
                        removedSet.begin(), removedSet.end(),
                        std::back_inserter(interior));
    return interior;
}

template<typename T>
gsSparseMatrix<T> deriveInnerEmbedding(const gsTensorBSplineBasis<2,T>& tensorBasis, boxSide side)
{
    std::vector<index_t> interior = getInteriorDofs(tensorBasis, side);
    // Convert to gsVector
    gsVector<index_t> interiorVec;
    interiorVec.assign(interior.begin(), interior.end());
    // Create embedding matrix
    return asEmbeddingMatrix<T>(tensorBasis.size(), interiorVec);
}

// ====================================================================
// Assemble overall embedding matrix
// ====================================================================

template<typename T>
struct gsAsG1Embedding {
public:
    gsSparseMatrix<T> matrix;
    // Sizes of blocks:
    //   0...interior dofs
    //   1...edge dofs level 0 (function values)
    //   2...edge dofs level 1 (crossing derivatives)
    gsVector<index_t,3> sizes;
};

template<typename T>
gsAsG1Embedding<T> createGluingDataArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    const gsMatrix<T> & localGluingData,
    T eps = 1e-12)
{
    gsAsG1Embedding<T> result;
    index_t rows = tensorBasis.size();

    const gsSparseMatrix<T> embeddingInterior = deriveInnerEmbedding(tensorBasis, side);
    result.sizes[0] = embeddingInterior.cols();

    const gsSparseMatrix<T> simpleEdgeEmbedding = gsBlockSparseMatrix<T>(1,2)
        .set(0,0,asEmbeddingMatrix<T>(rows, tensorBasis.boundary(side)))
        .set(0,1,asEmbeddingMatrix<T>(rows, tensorBasis.boundaryOffset(side,1)));

    const gsEdgeEmbedding<T> asG1edgeEmbedding = deriveEdgeEmbedding(tensorBasis, localGluingData, side);

    result.sizes[1] = asG1edgeEmbedding.sizes[0];
    result.sizes[2] = asG1edgeEmbedding.sizes[1];

    result.matrix = gsBlockSparseMatrix<T>(1,2)
         .set(0,0,embeddingInterior)
         .set(0,1,simpleEdgeEmbedding * asG1edgeEmbedding.matrix);

    return result;
}


} // namespace gismo
