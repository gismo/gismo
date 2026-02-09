/** @file embedding_example.cpp

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


int main(int argc, char* argv[])
{

    std::string coarseFilename("bspbasis/tpBSpline2_03.xml");
    std::string fineFilename("bspbasis/tpBSpline2_04.xml");
   // std::string output("");

    gsCmdLine cmd("Example for embedding matrix.");
    cmd.addString("c", "coarse", "G+Smo input coarse basis file.", coarseFilename);
    cmd.addString("f", "fine", "G+Smo input fine basis file.", fineFilename);
  //  cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the basis
    // ======================================================================

    gsFileData<> coarseFileData(coarseFilename);
    gsFileData<> fineFileData(fineFilename);

    gsBSplineBasis<>::uPtr pCoarseBasis;

    if (coarseFileData.has< gsBSplineBasis<> >())
    {
        pCoarseBasis = coarseFileData.getFirst< gsBSplineBasis<> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a coarse basis inside.\n";
        return -1;
    }


    if (!pCoarseBasis)
    {
        gsInfo << "Didn't find any coarse basis.\n";
        return -1;
    }

    gsBSplineBasis<>::uPtr pFineBasis;
    if (fineFileData.has< gsBSplineBasis<> >())
    {
        pFineBasis = fineFileData.getFirst< gsBSplineBasis<> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a fine basis inside.\n";
        return -1;
    }


    if (!pFineBasis)
    {
        gsInfo << "Didn't find any fine basis.\n";
        return -1;
    }


    // ======================================================================
    // checking if the fine basis is a refined coarse basis
    // ======================================================================
    gsBSplineBasis<>& coarseBasis = *pCoarseBasis;
    gsBSplineBasis<>& fineBasis = *pFineBasis;


    if(!isNested(coarseBasis,fineBasis)){

        gsInfo << "spaces are not nested\n";
        return -1;
    };

    gsInfo << "spaces are nested\n";


    // ======================================================================
    // embedding matrix computation
    // ======================================================================

    using T = double;
    gsSparseMatrix<T> embedding = embeddingMatrix(coarseBasis,fineBasis);
    gsInfo << "embedding\n"<< embedding << "\n";


    return 0;
}




