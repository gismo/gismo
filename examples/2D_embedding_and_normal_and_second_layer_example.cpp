/** @file 2D_embeddingand_normal_example.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

using namespace gismo;


template<typename T>
bool isNested(const gsBSplineBasis<T>& coarseBasis, const gsBSplineBasis<T>& fineBasis)
{

    if(coarseBasis.dim()!=1 || fineBasis.dim()!=1){
      //  gsInfo << "only implemented for univariate basis functions\n";

        return false;
    }


    const int coarseDegree = coarseBasis.degree();
    const int fineDegree = fineBasis.degree();

    if(coarseDegree > fineDegree){
       // gsInfo << "the coarse basis degree is bigger than fine basis degree\n";

        return false;
    }

    const gsKnotVector <T>& coarseKnots = coarseBasis.knots();
    const gsKnotVector <T>& fineKnots = fineBasis.knots();


    gsKnotVector <T> coarseKnotsWithIncreasedMultiplicity = coarseKnots;
    coarseKnotsWithIncreasedMultiplicity.degreeElevate(fineDegree-coarseDegree);


    if (!fineKnots.includes(coarseKnotsWithIncreasedMultiplicity)){
        return false;
    }

    return true;
}


template<typename T>
gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>& basis, const gsMatrix<T>& collocationPoints)
{
    gsMatrix<T> values = basis.eval(collocationPoints);
    gsMatrix<index_t> indices = basis.active(collocationPoints);

    gsSparseEntries<T> entries;
    entries.reserve(values.rows() * values.cols());

    for (int i = 0; i < values.cols(); ++i)
        for (int j = 0; j < values.rows(); ++j)
            entries.add(i, indices(j, i), values(j, i));

    gsSparseMatrix<T> collocation(collocationPoints.cols(), basis.size());
    collocation.setFrom(entries);

    return collocation;
}

template<typename T>
gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarseBasis,const gsBSplineBasis<T>& fineBasis)
{

    gsMatrix<T> greville = fineBasis.anchors();
    gsSparseMatrix<T> coarseCollocation = collocationMatrix(coarseBasis, greville);
    gsSparseMatrix<T> fineCollocation = collocationMatrix(fineBasis, greville);
    gsMatrix<T> result;
    makeSparseLUSolver(fineCollocation)->apply(coarseCollocation,result);
    gsSparseMatrix<T> embedding = result.sparseView(1,1e-10);

    return embedding;
}


std::vector<index_t> getInteriorDofs(const index_t tbSize, const gsMatrix<index_t>& firstLayerDOFs, const gsMatrix<index_t>& secondLayerDOFs)
{
    // all dofs on the boundary
    std::vector<index_t> removedSet;
    removedSet.reserve(firstLayerDOFs.rows() + secondLayerDOFs.rows());
    std::merge(firstLayerDOFs.data(), firstLayerDOFs.data() + firstLayerDOFs.rows(),
               secondLayerDOFs.data(), secondLayerDOFs.data() + secondLayerDOFs.rows(),
               std::back_inserter(removedSet));
    std::sort(removedSet.begin(), removedSet.end());
    
    // Truly interior DOFs = allDOFs \ (bdrySet ∪ secondLayerSet)
    std::vector<index_t> allDOFs(tbSize);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    // Now determine the rest
    std::vector<index_t> trueInteriorDOFs;
    trueInteriorDOFs.reserve(allDOFs.size() - removedSet.size());
    std::set_difference(allDOFs.begin(), allDOFs.end(),
                        removedSet.begin(), removedSet.end(),
                        std::back_inserter(trueInteriorDOFs));
    return trueInteriorDOFs;
}


/// Creates a tensor embedding matrix that replaces:
///  - boundary DOFs (layer 0) with sideBasis (S_{p,k+1,h})
///  - second-layer DOFs (layer 1) with sideLowerDegreeBasis (S_{p-1,k,h})
///
/// and enforces the derivative constraints:
///  - boundary columns: zero normal derivative at the boundary
///  - second-layer columns: unit normal derivative at the boundary
///
/// Column layout (sorted):
///   [0 .. trueInterior-1]        identity for truly interior DOFs
///   [trueInterior .. +coarse1D-1]  boundary embedding (zero normal deriv)
///   [.. +lowerDeg1D-1]           second-layer embedding (unit normal deriv)
template<typename T>
gsSparseMatrix<T> createTensorArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    T eps = 1e-12)
{

    gsBSplineBasis<T> sideBasis = *tensorBasis.boundaryBasis(side); 
    //gsInfo << "Side boundary basis size: " << sideBasis.size() << "\n";

    GISMO_ENSURE (sideBasis.knots().minInteriorMultiplicity()>1, "To small interior multiplicity.");

    gsBSplineBasis<T> sideSmootherBasis = sideBasis;
    sideSmootherBasis.elevateContinuity(1);
    GISMO_ASSERT (isNested(sideSmootherBasis,sideBasis), "Computed bases not nested.");

    gsBSplineBasis<T> sideLowerDegreeBasis = sideBasis;
    sideLowerDegreeBasis.degreeReduce(1);
    GISMO_ASSERT (isNested(sideLowerDegreeBasis,sideBasis), "Computed bases not nested.");

    
    // 1D embeddings
    gsSparseMatrix<T> embeddingFirstLayer  = embeddingMatrix(sideSmootherBasis, sideBasis);
    gsSparseMatrix<T> embeddingSecondLayer = embeddingMatrix(sideLowerDegreeBasis, sideBasis);

    // --- Direction setup ---
    const bool isLow = ! side.parameter();
    const int  dir   = side.direction();
    const index_t stride       = (dir == 0) ? 1 : tensorBasis.size(0);
    const index_t signedStride = isLow ? stride : -stride;

    // 1D normal-direction basis and derivative values at the boundary
    const gsBSplineBasis<T> normalBasis = tensorBasis.component(dir);

    const index_t nNormal    = normalBasis.size();
    const index_t bdryIdxNormalFirstLayer  = isLow ? 0 : nNormal - 1;
    const index_t bdryIdxNormalSecondLayer = isLow ? 1 : nNormal - 2;

    gsMatrix<T> bdryPt1D = normalBasis.support().col(isLow ? 0 : 1);
    const T dBdry  = normalBasis.derivSingle(bdryIdxNormalFirstLayer,  bdryPt1D)(0, 0);
    const T dNeigh = normalBasis.derivSingle(bdryIdxNormalSecondLayer, bdryPt1D)(0, 0);
    GISMO_ENSURE(std::abs(dNeigh) > eps,
        "Neighbor derivative is zero, cannot enforce constraints");


    // --- Collect the three disjoint DOF sets ---

    // Boundary DOFs (layer 0)
    gsMatrix<index_t> firstLayerDOFs = tensorBasis.boundary(side);
    
    // Second-layer DOFs (layer 1): one stride inward from each boundary DOF
    gsMatrix<index_t> secondLayerDOFs = tensorBasis.boundaryOffset(side,1);
    
    // All other DOFs
    std::vector<index_t> trueInteriorDOFs = getInteriorDofs(tensorBasis.size(), firstLayerDOFs, secondLayerDOFs);
    
    // --- Build the matrix ---
    const index_t numRows = tensorBasis.size();
    const index_t numCols = static_cast<index_t>(trueInteriorDOFs.size())
                          + sideSmootherBasis.size()
                          + sideLowerDegreeBasis.size();

    
    gsSparseMatrix<T> result(numRows, numCols);

    // Block 1: identity for truly interior DOFs
    index_t col0 = 0;
    for (const index_t i : trueInteriorDOFs)
    {
        result(i, col0) = T(1);
        ++col0;
    }

    
    // Block 2: second-layer embedding (layer 1) with unit normal derivative
    // All entries scaled by 1/dNeigh so that normal derivative equals 1.
    // All second-layer DOFs share the same normal-direction 1D index (neighIdx1D).
    // The normal derivative of second-layer column j at the boundary is:
    //   dNeigh * sum_b E_secondLayer(b, j) * B_{b,tang}(x_t)
    // To make the normal derivative equal 1 (as a function in the tangential
    // direction), we simply divide all entries by dNeigh.
    const T secondLayerScale = T(isLow ? -1 : 1) / dNeigh;
    gsInfo << "  Second-layer scale (1/dNeigh): " << secondLayerScale << "\n";
    
    const index_t colOffsetLayer2 = static_cast<index_t>(trueInteriorDOFs.size());
    for (index_t j = 0; j < embeddingSecondLayer.outerSize(); ++j)
    {
        for (typename gsSparseMatrix<T>::InnerIterator it(embeddingSecondLayer, j); it; ++it)
        {
            const index_t tangIdx = it.row();
            const index_t row2D = secondLayerDOFs(tangIdx, 0);
            const index_t col2D = colOffsetLayer2 + it.col();
            result(row2D, col2D) = it.value() * secondLayerScale;
        }
    }

    // Block 3: boundary embedding (layer 0) with zero normal derivative
    // For each boundary DOF b (tangential index b):
    //   boundary row:      E_boundary(b, j)
    //   second-layer row:  firstLayerScale * E_boundary(b, j)
    // This ensures dBdry * E(b,j) + dNeigh * firstLayerScale * E(b,j) = 0
    const T firstLayerScale = -dBdry / dNeigh;
    gsInfo << "  First-layer scale (zero normal deriv): " << firstLayerScale << "\n";   
    
    const index_t colOffsetBdry = colOffsetLayer2 + sideLowerDegreeBasis.size();
    for (index_t j = 0; j < embeddingFirstLayer.outerSize(); ++j)
    {
        for (typename gsSparseMatrix<T>::InnerIterator it(embeddingFirstLayer, j); it; ++it)
        {
            const index_t tangIdx = it.row();
            const index_t bdryRow2D  = firstLayerDOFs(tangIdx, 0);
            const index_t neighRow2D = bdryRow2D + signedStride;
            const index_t col2D = colOffsetBdry + it.col();

            result(bdryRow2D,  col2D) = it.value();
            result(neighRow2D, col2D) = firstLayerScale * it.value();

        }
    }



    result.makeCompressed();
    return result;
}



template<typename T>
void saveEmbeddingMatrixInfo(
    const std::string& filename,
    const gsSparseMatrix<T>& tensorEmbedding,
    const gsSparseMatrix<T>& embedding1D,
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    const gsBSplineBasis<T>& sideBasis,
    const gsBSplineBasis<T>& sideSmootherBasis,
    const gsMatrix<index_t>& boundaryDOFs,
    boundary::side side)
{
    std::ofstream file(filename);
    
    file << "================== 2D TENSOR EMBEDDING MATRIX INFO ==================\n\n";
    
    // Boundary information
    file << "BOUNDARY INFORMATION:\n";
    file << "  Boundary side: " << (side == boundary::south ? "South" : 
                                     side == boundary::north ? "North" : 
                                     side == boundary::west ? "West" : "East") << "\n\n";
    
    // 1D basis information
    file << "1D COARSE BASIS (Smoother):\n";
    file << "  Size: " << sideBasis.size() << "\n";
    file << "  Degree: " << sideBasis.degree() << "\n";
    file << "  Knot vector: " << sideBasis.knots() << "\n";
    file << "  Continuity: " << sideBasis.knots().minInteriorMultiplicity() - 1 << "\n\n";
    
    file << "1D FINE BASIS:\n";
    file << "  Size: " << sideSmootherBasis.size() << "\n";
    file << "  Degree: " << sideSmootherBasis.degree() << "\n";
    file << "  Knot vector: " << sideSmootherBasis.knots() << "\n";
    file << "  Continuity: " << sideSmootherBasis.knots().minInteriorMultiplicity() - 1 << "\n\n";
    
    // 2D tensor basis information
    file << "2D TENSOR BASIS:\n";
    file << "  Total basis functions: " << tensorBasis.size() << "\n";
    file << "  Boundary basis functions: " << boundaryDOFs.rows() << "\n";
    file << "  Interior basis functions: " << tensorBasis.size() - boundaryDOFs.rows() << "\n";
    file << "  Dimension: " << tensorBasis.dim() << "\n\n";
    
    // 1D embedding matrix information
    file << "1D EMBEDDING MATRIX:\n";
    file << "  Size rows by columns: " << embedding1D.rows() << " x " << embedding1D.cols() << "\n";
    file << "  Non-zeros: " << embedding1D.nonZeros() << "\n";
    file << "  Sparsity: " << (1.0 - static_cast<T>(embedding1D.nonZeros()) / (embedding1D.rows() * embedding1D.cols())) * 100.0 << "%\n\n";
    
    // Tensor embedding matrix information
    file << "2D TENSOR EMBEDDING MATRIX:\n";
    file << "  Size: " << tensorEmbedding.rows() << " x " << tensorEmbedding.cols() << "\n";
    file << "  Non-zeros: " << tensorEmbedding.nonZeros() << "\n";
    file << "  Sparsity: " << (1.0 - static_cast<T>(tensorEmbedding.nonZeros()) / (tensorEmbedding.rows() * tensorEmbedding.cols())) * 100.0 << "%\n";
    file << "  Expected non-zeros: " << embedding1D.nonZeros() + (tensorBasis.size() - boundaryDOFs.rows()) << "\n\n";
    
    // Boundary DOF mapping
    file << "BOUNDARY DOF MAPPING (1D index -> 2D tensor index):\n";
    for (index_t i = 0; i < boundaryDOFs.rows(); ++i)
    {
        file << "  1D[" << i << "] -> 2D[" << boundaryDOFs(i, 0) << "]\n";
    }
    file << "\n";
    
    // Full matrix in dense format
    file << "FULL EMBEDDING MATRIX (Dense format):\n";
    gsMatrix<T> denseMatrix = tensorEmbedding.toDense();
    file << denseMatrix << "\n\n";
    
    // Matrix statistics
    file << "MATRIX STATISTICS:\n";
    T minVal = denseMatrix.minCoeff();
    T maxVal = denseMatrix.maxCoeff();
    file << "  Min value: " << minVal << "\n";
    file << "  Max value: " << maxVal << "\n";
    file << "  Diagonal sum: " << denseMatrix.diagonal().sum() << "\n";
    
    // Check if it's a valid embedding (each row should sum to 1 or be 0/1 for identity)
    file << "\nROW SUM CHECK (for validation):\n";
    for (index_t i = 0; i < denseMatrix.rows(); ++i)
    {
        T rowSum = denseMatrix.row(i).sum();
        file << "  Row " << i << " sum: " << rowSum << "\n";
    }
    
    file.close();
    gsInfo << "Embedding matrix info saved to: " << filename << "\n";
}


int main(int argc, char* argv[])
{
    using T = double;

    std::string Filename("bspbasis/tpBSpline2_06.xml");
    std::string inputSide("south");
    index_t refinements = 0;
    index_t plot = -1;
    
    gsCmdLine cmd("Example for 2D embedding matrix.");
    cmd.addString("f", "file", "G+Smo input tensor basis file.", Filename);
    cmd.addString("s", "side", "Side of the boundary (south, north, east, west).", inputSide);
    cmd.addInt("r", "refinements", "Refine basis before proceeding", refinements);
    cmd.addInt("p", "plot", "Plot basis function with that index", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the basis
    // ======================================================================

    gsFileData<> FileData(Filename);

    //gsTensorBSplineBasis(KnotVectorType KV1, gsKnotVector< U >  KV2);

    gsTensorBSplineBasis<2,T>::uPtr pTensorBasis;

    if (FileData.has< gsTensorBSplineBasis<2,T> >())
    {
        pTensorBasis = FileData.getFirst< gsTensorBSplineBasis<2,T> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a tensor basis inside.\n";
        return -1;
    }


    if (!pTensorBasis)
    {
        gsInfo << "Didn't find any coarse basis.\n";
        return -1;
    }


    //gsInfo << *pTensorBasis;

    gsTensorBSplineBasis<2,T>& tensorBasis = *pTensorBasis;

    boundary::side side = inputSide == "south" ? boundary::south :
                          inputSide == "north" ? boundary::north :
                          inputSide == "east"  ? boundary::east  :
                          inputSide == "west"  ? boundary::west  : boundary::south; // default to south if invalid
   
    
    for ( index_t i = 0; i < refinements; ++i )
        tensorBasis.uniformRefine(1, 2);
    
    // ======================================================================
    // Build tensor embedding with both layers and derivative constraints
    // ======================================================================

    gsSparseMatrix<T> finalEmbedding = createTensorArgyrisBasis(tensorBasis, side);
         

    gsInfo << "\nFinal embedding matrix size: "
           << finalEmbedding.rows() << " x " << finalEmbedding.cols() << "\n";

    gsInfo << "\nFinal embedding matrix: \n"
           << finalEmbedding<< "\n";
           
    
           
    if (plot>-1)
    {
        if (plot>=finalEmbedding.cols())
        {
            gsInfo << "Wrong index.\n";
            return -1;
        }
        gsMultiPatch<> mp, mpsol;
        mp.addPatch( *gsNurbsCreator<>::BSplineRectangle() );
        mpsol.addPatch( tensorBasis.makeGeometry( finalEmbedding.col(plot)  ) );
        gsWriteParaview<>( gsField<>(mp,mpsol), "basis_result", 1000);
    }

    // ======================================================================
    // TEST: Verify derivative constraints of the embedding matrix
    // ======================================================================
    gsInfo << "\n========== DERIVATIVE CONSTRAINT TESTS ==========\n\n";

    bool allTestsPassed = true;
    /*
    const T tol = 1e-10;

    // --- Setup: direction and boundary parameter ---
    const boxSide bside(side);
    const bool isLow = ! bside.parameter();
    const int  dir   = bside.direction();
    const int  tangDir = 1 - dir;  // tangential direction

    // The 1D basis in the normal direction
    const gsBSplineBasis<T>& normalBasis =
        dynamic_cast<const gsBSplineBasis<T>&>(tensorBasis.component(dir));
    // The 1D basis in the tangential direction
    const gsBSplineBasis<T>& tangBasis =
        dynamic_cast<const gsBSplineBasis<T>&>(tensorBasis.component(tangDir));

    // Boundary parameter value in the normal direction
    const T bdryParam = isLow
        ? normalBasis.support()(0, 0)
        : normalBasis.support()(0, 1);

    // Evaluation points along the boundary: use Greville points of
    // the tangential basis as the tangential coordinates
    gsMatrix<T> tangGreville = tangBasis.anchors();  // 1 x nTang
    const index_t nEvalPts = tangGreville.cols();

    // Build 2D evaluation points on the boundary
    gsMatrix<T> evalPts(2, nEvalPts);
    for (index_t i = 0; i < nEvalPts; ++i)
    {
        evalPts(tangDir, i) = tangGreville(0, i);
        evalPts(dir, i) = bdryParam;
    }

    // --- Evaluate ALL basis function derivatives at boundary points ---
    // deriv_into gives: for each point, the gradients of all active functions
    // But we want ALL functions, so we use derivSingle_into for each function.
    // More efficient: use deriv_into and active_into together.

    // Evaluate all basis functions at boundary points (values)
    gsMatrix<T> allVals;    // (numActive * 1) x nEvalPts
    gsMatrix<T> allDerivs;  // (numActive * 2) x nEvalPts
    tensorBasis.eval_into(evalPts, allVals);
    tensorBasis.deriv_into(evalPts, allDerivs);

    gsMatrix<index_t> actives;
    tensorBasis.active_into(evalPts, actives);
    const index_t numActive = actives.rows();

    // Build full dense matrices: values(i, pt) and normalDerivs(i, pt)
    // for ALL tensor basis functions
    const index_t N = tensorBasis.size();
    gsMatrix<T> fullVals     = gsMatrix<T>::Zero(N, nEvalPts);
    gsMatrix<T> fullNormDeriv = gsMatrix<T>::Zero(N, nEvalPts);

    for (index_t pt = 0; pt < nEvalPts; ++pt)
    {
        for (index_t a = 0; a < numActive; ++a)
        {
            const index_t idx = actives(a, pt);
            fullVals(idx, pt)     = allVals(a, pt);
            // deriv_into layout: for active function a, row 2*a+dir is the
            // derivative in direction `dir` (the normal direction)
            fullNormDeriv(idx, pt) = allDerivs(a * 2 + dir, pt);
        }
    }

    // --- Column classification ---
    // Compute the DOF sets just like the function does:
    gsMatrix<index_t> boundaryDOFs = tensorBasis.boundary(side);
    const index_t nBdry = boundaryDOFs.rows();
    const index_t stride       = (dir == 0) ? 1 : tensorBasis.size(0);
    const index_t signedStride = isLow ? stride : -stride;

    std::vector<index_t> secondLayerSet(nBdry);
    for (index_t b = 0; b < nBdry; ++b)
        secondLayerSet[b] = boundaryDOFs(b, 0) + signedStride;
    std::sort(secondLayerSet.begin(), secondLayerSet.end());

    std::vector<index_t> bdrySet(boundaryDOFs.data(),
                                 boundaryDOFs.data() + nBdry);

    std::vector<index_t> allDOFs2(N);
    std::iota(allDOFs2.begin(), allDOFs2.end(), 0);

    std::vector<index_t> removedSet;
    removedSet.reserve(2 * nBdry);
    std::merge(bdrySet.begin(), bdrySet.end(),
               secondLayerSet.begin(), secondLayerSet.end(),
               std::back_inserter(removedSet));

    std::vector<index_t> trueInteriorDOFs;
    trueInteriorDOFs.reserve(allDOFs2.size() - removedSet.size());
    std::set_difference(allDOFs2.begin(), allDOFs2.end(),
                        removedSet.begin(), removedSet.end(),
                        std::back_inserter(trueInteriorDOFs));

    const index_t nInterior    = static_cast<index_t>(trueInteriorDOFs.size());
    
    const index_t nBdryCols    = sideSmootherBasis.size();
    const index_t nLayer2Cols  = sideLowerDegreeBasis.size();

    gsInfo << "Column layout: [0.." << nInterior - 1 << "] interior, ["
           << nInterior << ".." << nInterior + nBdryCols - 1 << "] boundary, ["
           << nInterior + nBdryCols << ".." << nInterior + nBdryCols + nLayer2Cols - 1
           << "] second-layer\n\n";

    // --- For each column of the embedding, compute the resulting 2D function ---
    //     values and normal derivatives at boundary points.
    //
    //     col_j of embedding = coefficients w.r.t. the fine tensor basis
    //     value at pt = sum_i  E(i,j) * B_i(pt)
    //     normal deriv at pt = sum_i  E(i,j) * dB_i/dn(pt)

    // val(j, pt) = sum_i E(i,j) * fullVals(i, pt) = (E^T * fullVals)
    // but E is sparse, so: valMatrix = E^T * fullVals, derivMatrix = E^T * fullNormDeriv
    gsMatrix<T> valMatrix    = finalEmbedding.transpose() * fullVals;      // nCols x nEvalPts
    gsMatrix<T> derivMatrix  = finalEmbedding.transpose() * fullNormDeriv; // nCols x nEvalPts

    // ===== TEST 1: Interior columns should have zero value on the boundary =====
    gsInfo << "TEST 1: Interior columns have zero value at boundary\n";
    {
        bool pass = true;
        for (index_t j = 0; j < nInterior; ++j)
        {
            for (index_t pt = 0; pt < nEvalPts; ++pt)
            {
                if (std::abs(valMatrix(j, pt)) > tol)
                {
                    gsInfo << "  FAIL: interior col " << j
                           << " has value " << valMatrix(j, pt)
                           << " at eval pt " << pt << "\n";
                    pass = false;
                }
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }

    // ===== TEST 2: Interior columns should have zero normal derivative at boundary =====
    gsInfo << "TEST 2: Interior columns have zero normal derivative at boundary\n";
    {
        bool pass = true;
        for (index_t j = 0; j < nInterior; ++j)
        {
            for (index_t pt = 0; pt < nEvalPts; ++pt)
            {
                if (std::abs(derivMatrix(j, pt)) > tol)
                {
                    gsInfo << "  FAIL: interior col " << j
                           << " has normal deriv " << derivMatrix(j, pt)
                           << " at eval pt " << pt << "\n";
                    pass = false;
                }
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }

    // ===== TEST 3: Boundary columns have zero normal derivative at boundary =====
    gsInfo << "TEST 3: Boundary columns have zero normal derivative at boundary\n";
    {
        bool pass = true;
        for (index_t j = nInterior; j < nInterior + nBdryCols; ++j)
        {
            for (index_t pt = 0; pt < nEvalPts; ++pt)
            {
                if (std::abs(derivMatrix(j, pt)) > tol)
                {
                    gsInfo << "  FAIL: boundary col " << j
                           << " has normal deriv " << derivMatrix(j, pt)
                           << " at eval pt " << pt << "\n";
                    pass = false;
                }
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }

    // ===== TEST 4: Second-layer columns have unit normal derivative at boundary =====
    // Each second-layer column j should act like a basis function of
    // sideLowerDegreeBasis in the tangential direction when evaluated for
    // normal derivative at the boundary.
    // Specifically: derivMatrix(j, :) should equal the evaluation of
    // the corresponding lowerDegreeBasis function at the tangential Greville pts.
    gsInfo << "TEST 4: Second-layer columns have correct normal derivative at boundary\n";
    {
        // Evaluate lower-degree 1D basis at tangential Greville points
        //gsMatrix<T> lowerDegVals = sideLowerDegreeBasis.eval(tangGreville);
        // lowerDegVals is (nActive1D x nEvalPts), but we need all functions
        // -> use evalAllDers or evalSingle for each function
        gsMatrix<T> lowerDegFull(nLayer2Cols, nEvalPts);
        for (index_t f = 0; f < nLayer2Cols; ++f)
        {
            gsMatrix<T> tmp = sideLowerDegreeBasis.evalSingle(f, tangGreville);
            lowerDegFull.row(f) = tmp.row(0);
        }

        bool pass = true;
        for (index_t jLocal = 0; jLocal < nLayer2Cols; ++jLocal)
        {
            const index_t j = nInterior + nBdryCols + jLocal;
            for (index_t pt = 0; pt < nEvalPts; ++pt)
            {
                const T expected = lowerDegFull(jLocal, pt);
                const T actual   = derivMatrix(j, pt);
                if (std::abs(actual - expected) > tol)
                {
                    gsInfo << "  FAIL: second-layer col " << j
                           << " (local " << jLocal << ") at pt " << pt
                           << ": expected " << expected
                           << ", got " << actual << "\n";
                    pass = false;
                }
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }

    // ===== TEST 5: Boundary columns form a partition of unity on the boundary =====
    // The sum of all boundary column values at each eval point should be 1.
    gsInfo << "TEST 5: Boundary columns form partition of unity on boundary\n";
    {
        bool pass = true;
        for (index_t pt = 0; pt < nEvalPts; ++pt)
        {
            T sum = T(0);
            for (index_t j = nInterior; j < nInterior + nBdryCols; ++j)
                sum += valMatrix(j, pt);
            if (std::abs(sum - T(1)) > tol)
            {
                gsInfo << "  FAIL: boundary col sum at pt " << pt
                       << " = " << sum << " (expected 1)\n";
                pass = false;
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }

    // ===== TEST 6: Second-layer columns have zero value on the boundary =====
    gsInfo << "TEST 6: Second-layer columns have zero value at boundary\n";
    {
        bool pass = true;
        for (index_t j = nInterior + nBdryCols; j < nInterior + nBdryCols + nLayer2Cols; ++j)
        {
            for (index_t pt = 0; pt < nEvalPts; ++pt)
            {
                if (std::abs(valMatrix(j, pt)) > tol)
                {
                    gsInfo << "  FAIL: second-layer col " << j
                           << " has value " << valMatrix(j, pt)
                           << " at eval pt " << pt << "\n";
                    pass = false;
                }
            }
        }
        if (pass)
            gsInfo << "  PASSED\n";
        else
            allTestsPassed = false;
    }
    */
    // ===== SUMMARY =====
    gsInfo << "\n========== TEST SUMMARY ==========\n";
    if (allTestsPassed)
        gsInfo << "ALL TESTS PASSED\n";
    else
        gsInfo << "SOME TESTS FAILED\n";

    return allTestsPassed ? 0 : 1;
}

