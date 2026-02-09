/** @file 2Dembedding_example.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
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

    gsMatrix<size_t> indices = basis.active(collocationPoints);

    //TODO: use gsSparseEntries

    gsSparseMatrix<T> collocation(collocationPoints.cols(), basis.size());
    for(int i = 0; i < values.cols(); i++){
        for(int j = 0; j < values.rows(); j++){
            collocation(i,indices(j,i)) = values(j,i);
        }
    }

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


/*
// Create a transformation matrix for the full 2D tensor basis
template<typename T>
gsSparseMatrix<T> createTensorEmbedding(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    const gsBSplineBasis<T>& coarse1DBasis,
    const gsBSplineBasis<T>& fine1DBasis,
    boundary::side side)
{
    // Get the embedding matrix for the 1D boundary
    gsSparseMatrix<T> embedding1D = embeddingMatrix(coarse1DBasis, fine1DBasis);
    
    // Get boundary DOF mapping
    gsMatrix<index_t> boundaryDOFs = tensorBasis.boundary(side);
    index_t numTensorBasisFuncs = tensorBasis.size();
    
    // Create the full transformation matrix
    gsSparseMatrix<T> result(numTensorBasisFuncs, numTensorBasisFuncs);
    result.setIdentity(); // Initialize with identity
    
    // Insert the 1D embedding into the appropriate block
    for (index_t i = 0; i < embedding1D.rows(); ++i)
    {
        for (index_t j = 0; j < embedding1D.cols(); ++j)
        {
            if (embedding1D.coeff(i, j) != 0)
            {
                // Map 1D indices to 2D tensor indices
                index_t row2D = boundaryDOFs(i, 0);
                index_t col2D = boundaryDOFs(j, 0);
                result.coeffRef(row2D, col2D) = embedding1D.coeff(i, j);
            }
        }
    }
    
    return result;
}

*/




template<typename T>
gsSparseMatrix<T> createTensorEmbedding(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    const gsBSplineBasis<T>& coarse1DBasis,
    const gsBSplineBasis<T>& fine1DBasis,
    boundary::side side)
{
    gsSparseMatrix<T> embedding1D = embeddingMatrix(coarse1DBasis, fine1DBasis);
    
    gsMatrix<index_t> boundaryDOFs = tensorBasis.boundary(side);
    
    index_t numFineTensorFuncs = tensorBasis.size();  // rows (fine basis)
    index_t numCoarseTensorFuncs = tensorBasis.size() - boundaryDOFs.rows() + coarse1DBasis.size();  // cols (coarse basis)
    
    gsSparseMatrix<T> result(numFineTensorFuncs, numCoarseTensorFuncs);
    
    //  identity matrix for interior basis functions
    index_t coarseIdx = 0;
    for (index_t i = 0; i < numFineTensorFuncs; ++i)
    {
        // Check if this is a boundary DOF
        bool isBoundary = false;
        for (index_t b = 0; b < boundaryDOFs.rows(); ++b)
        {
            if (boundaryDOFs(b, 0) == i)
            {
                isBoundary = true;
                break;
            }
        }
        
        if (!isBoundary)
        {
            result.coeffRef(i, coarseIdx) = 1.0;
            coarseIdx++;
        }
    }
    
    // Insert the 1D embedding for boundary functions
    for (index_t i = 0; i < embedding1D.rows(); ++i)
    {
        for (index_t j = 0; j < embedding1D.cols(); ++j)
        {
            if (embedding1D.coeff(i, j) != 0)
            {
                // Map 1D indices to 2D tensor indices
                index_t row2D = boundaryDOFs(i, 0);
                index_t col2D = coarseIdx + j;  // Offset by interior functions
                result.coeffRef(row2D, col2D) = embedding1D.coeff(i, j);
            }
        }
    }
    
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

    
   
    // ======================================================================
    // Debug: Check tensor embedding matrix
    // ======================================================================
    
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

    return 0;
}




