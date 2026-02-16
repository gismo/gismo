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


// Returns the string with the size of a matrix.
template <typename T>
std::string size(const gsMatrix<T>& matrix)
{
    std::string result = "(" + util::to_string(matrix.rows()) + " x " +
        util::to_string(matrix.cols()) + ")";

    return result;
}

template<typename T>
void print(const std::vector<T>& v)
{
    std::cout << "[ ";
    for (std::size_t i=0; i<v.size(); ++i)
        std::cout << v[i] << " ";
    std::cout << "]";
}

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




template<typename T>
gsSparseMatrix<T> createTensorEmbedding(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    const gsBSplineBasis<T>& coarse1DBasis,
    const gsBSplineBasis<T>& fine1DBasis,
    boundary::side side)
{
    gsSparseMatrix<T> embedding1D = embeddingMatrix(coarse1DBasis, fine1DBasis);

    gsMatrix<index_t> boundaryDOFs = tensorBasis.boundary(side);

    std::vector<index_t> bdrySet(boundaryDOFs.data(),
                              boundaryDOFs.data() + boundaryDOFs.rows());

    // Compute interior DOFs via set difference: allDOFs \ bdrySet
    std::vector<index_t> allDOFs(tensorBasis.size());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    std::vector<index_t> interiorDOFs;
    interiorDOFs.reserve(allDOFs.size() - bdrySet.size());
    std::set_difference(allDOFs.begin(), allDOFs.end(),
                        bdrySet.begin(), bdrySet.end(),
                        std::back_inserter(interiorDOFs));

    const index_t numFineTensorFuncs = tensorBasis.size();
    const index_t numCoarseTensorFuncs =
        static_cast<index_t>(interiorDOFs.size()) + coarse1DBasis.size();

    gsSparseMatrix<T> result(numFineTensorFuncs, numCoarseTensorFuncs);

    // Identity block for interior basis functions
    index_t coarseIdx = 0;
    for (const index_t i : interiorDOFs)
    {
        result(i, coarseIdx) = 1.0;
        ++coarseIdx;
    }
    
    // Insert the 1D embedding for boundary functions
    // Iterate only over nonzero entries of embedding1D (column-major)
    for (index_t j = 0; j < embedding1D.outerSize(); ++j)
    {
        for (typename gsSparseMatrix<T>::InnerIterator it(embedding1D, j); it; ++it)
        {
            // it.row() = fine 1D index, j = coarse 1D index
            const index_t row2D = boundaryDOFs(it.row(), 0);
            const index_t col2D = coarseIdx + it.col();
            result(row2D, col2D) = it.value();
        }
    }

    return result;
}


template<typename T>
gsSparseMatrix<T> enforceZeroNormalDerivative(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    const gsSparseMatrix<T>& tensorEmbedding,
    const gsMatrix<index_t>& boundaryDOFs,
    boxSide side,
    T eps = 1e-12)
{
    // Determine direction properties once
    const bool isLow = ! side.parameter();
    GISMO_ASSERT ( isLow == (side == boundary::south || side == boundary::west), "");
    
    const int dir = side.direction();
    GISMO_ASSERT ( dir == ((side == boundary::south || side == boundary::north) ? 1 : 0), "");
    
    
    // Get the 1D basis in the normal direction
    const gsBSplineBasis<T>& normalBasis =
        dynamic_cast<const gsBSplineBasis<T>&>(tensorBasis.component(dir));

    const index_t nNormal = normalBasis.size();
    const index_t bdryIdx1D  = isLow ? 0 : nNormal - 1;
    const index_t neighIdx1D = isLow ? 1 : nNormal - 2;

    // Evaluate 1D derivatives at the boundary parameter
    
    gsMatrix<T> bdryPt = normalBasis.support().col(isLow ? 0 : 1);
    //  gsMatrix<T> bdryPt(1, 1);
    //  bdryPt(0, 0) = isLow ? T(0) : T(1);

    const T dBdry  = normalBasis.derivSingle(bdryIdx1D, bdryPt)(0, 0);
    const T dNeigh = normalBasis.derivSingle(neighIdx1D, bdryPt)(0, 0);

    GISMO_ENSURE (std::abs(dNeigh) > eps, "Neighbor derivative is zero, cannot enforce zero normal derivative");

    const T ratio = -dBdry / dNeigh;
    gsInfo << "  Correction ratio: " << ratio << "\n";

    // Stride: offset between a boundary DOF and its neighbor in the normal direction
    const index_t stride = (dir == 0) ? 1 : tensorBasis.size(0);
    const index_t signedStride = isLow ? stride : -stride;

    gsSparseMatrix<T> result = tensorEmbedding;

    // For each boundary DOF, set the neighbor row = ratio * boundary row
    for (index_t b = 0; b < boundaryDOFs.rows(); ++b)
    {
        const index_t fineRow  = boundaryDOFs(b, 0);
        const index_t neighbor = fineRow + signedStride;

        // Iterate over columns that have a nonzero at fineRow
        for (index_t col = 0; col < result.outerSize(); ++col)
        {
            const T val = result.coeff(fineRow, col);
            if (std::abs(val) < eps) continue;

            result(neighbor, col) = ratio * val;
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
    const gsBSplineBasis<T>& coarse1DBasis,
    const gsBSplineBasis<T>& fine1DBasis,
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
    file << "  Size: " << coarse1DBasis.size() << "\n";
    file << "  Degree: " << coarse1DBasis.degree() << "\n";
    file << "  Knot vector: " << coarse1DBasis.knots() << "\n";
    file << "  Continuity: " << coarse1DBasis.knots().minInteriorMultiplicity() - 1 << "\n\n";
    
    file << "1D FINE BASIS:\n";
    file << "  Size: " << fine1DBasis.size() << "\n";
    file << "  Degree: " << fine1DBasis.degree() << "\n";
    file << "  Knot vector: " << fine1DBasis.knots() << "\n";
    file << "  Continuity: " << fine1DBasis.knots().minInteriorMultiplicity() - 1 << "\n\n";
    
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

   // std::string output("");

    std::string inputSide="south";
    gsCmdLine cmd("Example for 2D embedding matrix.");
    cmd.addString("f", "file", "G+Smo input tensor basis file.", Filename);
    cmd.addString("s", "side", "Side of the boundary (south, north, east, west).", inputSide);
  //  cmd.addString("o", "output", "Name of the output file.", output);

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
   
    gsBSplineBasis<T> sideBasis = *tensorBasis.boundaryBasis(side); 
   
    gsInfo << "Side boundary basis size: " << sideBasis.size() << "\n";
   

    if( sideBasis.knots().minInteriorMultiplicity()<2){
        gsInfo<<"incorrect multiplicity\n";
        return -1;
    }
    
   
    gsBSplineBasis<T> sideSmootherBasis = sideBasis;
    sideSmootherBasis.elevateContinuity(1);


    // ======================================================================
    // checking if the fine basis is a refined coarse basis
    // ======================================================================

    if(!isNested(sideSmootherBasis,sideBasis)){

        gsInfo << "spaces are not nested\n";
        return -1;
    };

    gsInfo << "spaces are nested\n";


    // ======================================================================
    // embedding matrix computation
    // ======================================================================


    using T = double;
    gsSparseMatrix<T> embedding = embeddingMatrix(sideSmootherBasis, sideBasis);
    gsInfo << "embedding\n"<< embedding << "\n";

    gsMatrix<index_t> sideDOFs = tensorBasis.boundary(side);

    gsInfo << "Side boundary DOFs: " << sideDOFs.rows() << " functions\n";

    // Print the indices
    for (index_t i = 0; i < sideDOFs.rows(); ++i) {
        gsInfo << "DOF " << i << ": " << sideDOFs(i, 0) << "\n";
    }

    gsSparseMatrix<T> tensorEmbedding = createTensorEmbedding(tensorBasis, sideSmootherBasis, sideBasis, side);
    gsInfo << "\nTensor embedding matrix:\n" << tensorEmbedding << "\n";
    gsInfo << "Modified embedding matrix size: " << tensorEmbedding.rows() << " x " << tensorEmbedding.cols() << "\n";
   
    // ======================================================================
    // Debug: Check tensor embedding matrix
    // ======================================================================
   /* 
    index_t numTensorBasisFuncs = tensorBasis.size();
    
    gsInfo << "\n=== Tensor Embedding Verification ===\n";
    gsInfo << "1D embedding matrix size: " << embedding.rows() << " x " << embedding.cols() << "\n";
    gsInfo << "1D embedding non-zeros: " << embedding.nonZeros() << "\n";
    gsInfo << "Boundary DOFs count: " << sideDOFs.rows() << "\n";
    gsInfo << "Total tensor basis functions: " << numTensorBasisFuncs << "\n";
    gsInfo << "Tensor embedding matrix size: " << tensorEmbedding.rows() << " x " << tensorEmbedding.cols() << "\n";
    gsInfo << "Tensor embedding non-zeros: " << tensorEmbedding.nonZeros() << "\n";
    gsInfo << "Expected non-zeros: " << embedding.nonZeros() + (numTensorBasisFuncs - sideDOFs.rows()) << "\n";
    
    gsInfo << "\nTensor embedding matrix:\n" << tensorEmbedding << "\n";

    saveEmbeddingMatrixInfo(
        "../examples/embedding_matrix_info.txt",
        tensorEmbedding,
        embedding,
        tensorBasis,
        sideSmootherBasis,
        sideBasis,
        sideDOFs,
        side);

*/


    // ======================================================================
    // Enforce zero normal derivative at boundary
    // ======================================================================
    
    gsInfo << "\n=== Enforcing Zero Normal Derivative ===\n";
    
    gsSparseMatrix<T> modifiedEmbedding = enforceZeroNormalDerivative(
        tensorBasis,
        tensorEmbedding,
        sideDOFs,
        side);
    
    gsInfo << "Modified embedding matrix size: " << modifiedEmbedding.rows() << " x " << modifiedEmbedding.cols() << "\n";
    //gsInfo << "Modified embedding non-zeros: " << modifiedEmbedding.nonZeros() << "\n";
   // gsInfo << "Original embedding non-zeros: " << tensorEmbedding.nonZeros() << "\n";
    //gsInfo << "Added non-zeros: " << modifiedEmbedding.nonZeros() - tensorEmbedding.nonZeros() << "\n";
    gsInfo << "\nModified embedding matrix:\n" << modifiedEmbedding << "\n";

   
    gsInfo << "\n=== Verifying Zero Normal Derivative ===\n";
    
    int dir = (side == boundary::south || side == boundary::north) ? 1 : 0;
    
    const index_t nBoundary = sideDOFs.rows();

    // Get the boundary parameter value in the normal direction
    T bdryParam = (side == boundary::south || side == boundary::west) ? T(0) : T(1);
    
    // For verification, we evaluate derivatives of each basis function at multiple
    // boundary points and check the normal derivative of each column of the modified embedding
    
    // Use the Greville abscissae of the tangential basis as test points
    int tangDir = 1 - dir;
    const gsBSplineBasis<T>& tangBasis =
        dynamic_cast<const gsBSplineBasis<T>&>(tensorBasis.component(tangDir));
    gsMatrix<T> tangAnchors = tangBasis.anchors();
    
    T maxOverall = 0;
    
    for (index_t col = 0; col < modifiedEmbedding.cols(); ++col)
    {
        // Check if this column involves any boundary DOF
        bool hasBoundary = false;
        for (index_t b = 0; b < nBoundary; ++b)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(modifiedEmbedding, col); it; ++it)
            {
                if (it.row() == sideDOFs(b, 0))
                {
                    hasBoundary = true;
                    break;
                }
            }
            if (hasBoundary) break;
        }
        
        if (!hasBoundary) continue;
        
        // Check the normal derivative at several points along the boundary
        T maxForCol = 0;
        for (index_t t = 0; t < tangAnchors.cols(); ++t)
        {
            gsMatrix<T> pt(2, 1);
            pt(tangDir, 0) = tangAnchors(0, t);
            pt(dir, 0) = bdryParam;
            
            // Compute normal derivative: sum over all basis functions i with nonzero coefficient
            T normalDeriv = 0;
            for (typename gsSparseMatrix<T>::InnerIterator it(modifiedEmbedding, col); it; ++it)
            {
                gsMatrix<T> d = tensorBasis.derivSingle(it.row(), pt);
                // derivSingle returns (dim x nPoints), so d(dir, 0) is the normal derivative
                normalDeriv += it.value() * d(dir, 0);
            }
            maxForCol = std::max(maxForCol, std::abs(normalDeriv));
        }
        
        if (maxForCol > 1e-10)
        {
            gsInfo << "  Column " << col << ": max |normal deriv| = " << maxForCol << " *** NOT ZERO\n";
        }
        maxOverall = std::max(maxOverall, maxForCol);
    }
    gsInfo << "  Max normal derivative across all boundary columns: " << maxOverall << "\n";
    if (maxOverall < 1e-10)
        gsInfo << "  SUCCESS: All normal derivatives are zero at the boundary.\n";
    else
        gsInfo << "  FAILURE: Some normal derivatives are NOT zero.\n";

    return 0;
}



