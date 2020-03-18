//
// Created by afarahat on 2/25/20.
//

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include <gsG1Basis/gsG1AuxiliaryPatch.h>
# include <gsG1Basis/gsG1BasisEdge.h>
# include <gsG1Basis/gsG1BasisVertex.h>


namespace gismo
{

class gsG1AuxiliaryVertexMultiplePatches
{

public:

// Constructor for n patches around a common vertex
    gsG1AuxiliaryVertexMultiplePatches(const gsMultiPatch<> & mp, const std::vector<size_t> patchesAroundVertex, const std::vector<size_t> vertexIndices)
    {
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(patchesAroundVertex[i]), patchesAroundVertex[i]));
            auxVertexIndices.push_back(vertexIndices[i]);
            checkBoundary(mp, patchesAroundVertex[i], i);
        }
        sigma = 0.0;
        gsInfo << "\n";
    }

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    gsMultiPatch<> computeAuxTopology(){
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++)
        {
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        auxTop.computeTopology();
        return auxTop;
    }

    void reparametrizeG1Vertex()
    {
        for(size_t i = 0; i < auxGeom.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            switch (auxVertexIndices[i])
            {
                case 1:
                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " not rotated\n";
                    break;
                case 4:
                    auxGeom[i].rotateParamAntiClockTwice();
                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex()
                           << " rotated twice anticlockwise\n";
                    break;
                case 2:
                    auxGeom[i].rotateParamAntiClock();
                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " rotated anticlockwise\n";
                    break;
                case 3:
                    auxGeom[i].rotateParamClock();
                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " rotated clockwise\n";
                    break;
            }
        }
        gsInfo << "-----------------------------------------------------------------\n";
    }


    index_t kindOfVertex()
    {
        if(auxGeom.size() == 1)
            return -1; // Boundary vertex

        gsMultiPatch<> top(computeAuxTopology());
        size_t nInt = top.interfaces().size();
        if(auxGeom.size() == nInt)
            return 0; // Internal vertex
        else
            return 1; // Interface-Boundary vertex
    }



    void checkOrientation(size_t i)
    {
        if (auxGeom[i].getPatch().orientation() == -1)
        {
            auxGeom[i].swapAxis();
            gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";

            this->swapBdy(i); //Swap boundary edge bool-value

            // Swap vertices index after swapping axis
            if(auxVertexIndices[i] == 2)
                auxVertexIndices[i] = 3;
            else
            if(auxVertexIndices[i] == 3)
                auxVertexIndices[i] = 2;

        }
    }


    void computeSigma()
    {
        real_t p = 0;
        real_t h_geo = 0;
        for(size_t i = 0; i < auxGeom.size(); i++)
        {
            gsTensorBSplineBasis<2, real_t> & bsp_temp = dynamic_cast<gsTensorBSplineBasis<2, real_t> & >(auxGeom[i].getPatch().basis());
            real_t p_temp = bsp_temp.maxDegree();
            p = (p < p_temp ? p_temp : p);

            for(index_t j = 0; j < auxGeom[i].getPatch().parDim(); j++)
            {
               real_t h_geo_temp = bsp_temp.component(j).knots().at(p + 2);
               h_geo = (h_geo < h_geo_temp ? h_geo_temp : h_geo);
            }

        }
        real_t val = auxGeom.size();

        gsMatrix<> zero;
        zero.setZero(2,1);
        for (index_t i = 0; i < val; i++)
            sigma += auxGeom[i].getPatch().deriv(zero).lpNorm<Eigen::Infinity>();
        sigma *= h_geo/(val*p);
        sigma = 1 / sigma;
    }


    void checkBoundary(const gsMultiPatch<> & mpTmp, size_t  patchInd, size_t Ind)
    {
        std::vector<bool> tmp;
        switch (auxVertexIndices[Ind])
        {
            case 1: tmp.push_back(mpTmp.isBoundary(patchInd,3));
                    tmp.push_back(mpTmp.isBoundary(patchInd,1));
                    gsInfo << "Edge 3: " << mpTmp.isBoundary(patchInd, 3) << "\t Edge 1: " << mpTmp.isBoundary(patchInd, 1) << "\n";
                break;
            case 2: tmp.push_back(mpTmp.isBoundary(patchInd, 2));
                    tmp.push_back(mpTmp.isBoundary(patchInd, 3));
                    gsInfo << "Edge 2: " << mpTmp.isBoundary(patchInd, 2) << "\t Edge 3: " << mpTmp.isBoundary(patchInd, 3) << "\n";

                break;
            case 3: tmp.push_back(mpTmp.isBoundary(patchInd, 1));
                    tmp.push_back(mpTmp.isBoundary(patchInd, 4));
                    gsInfo << "Edge 1: " << mpTmp.isBoundary(patchInd, 1) << "\t Edge 4: " << mpTmp.isBoundary(patchInd, 4) << "\n";

                break;
            case 4: tmp.push_back(mpTmp.isBoundary(patchInd, 4));
                    tmp.push_back(mpTmp.isBoundary(patchInd, 2));
                    gsInfo << "Edge 4: " << mpTmp.isBoundary(patchInd, 4) << "\t Edge 2: " << mpTmp.isBoundary(patchInd, 2) << "\n";
                break;
            default:
                break;
        }
        isBdy.push_back(tmp);
    }


    void swapBdy(size_t i)
    {
        bool tmp = isBdy[i][0];
        isBdy[i][0] = isBdy[i][1];
        isBdy[i][1] = tmp;
    }


    gsMatrix<> computeBigSystemMatrix( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);
        size_t dimV = temp_L.size(1);

        gsMatrix<> BigMatrix;
        BigMatrix.setZero( 2 * (dimU + dimV - 2),6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < 2 * dimU; i++)
            {
                BigMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (size_t i = 1; i < dimV - 1; i++)
            {
                for(size_t j = i; j < i + 2; j++)
                {
                    BigMatrix(i + j + (2 * dimU ) - 2, bf) = auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i);
                }
            }
        }
        return BigMatrix;
    }


    gsMatrix<> computeSmallSystemMatrix( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);
        size_t dimV = temp_L.size(1);

        gsMatrix<> SmallMatrix;
        SmallMatrix.setZero((dimU + dimV - 1),6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < dimU; i++)
            {
                SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (size_t i = 1; i < dimV; i++)
            {
                SmallMatrix(i + dimU -1, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i * dimU);
            }
        }
        return SmallMatrix;
    }


    gsMatrix<> leftBoundaryBigSystem(index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);
        size_t dimV = temp_L.size(1);

        gsMatrix<> BigMatrix;
        BigMatrix.setZero( 2 * dimV,6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < dimV ; i++)
            {
                for(size_t j = i; j < i + 2; j++)
                {
                    BigMatrix(i + j, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i  * dimU + j - i);
                }
            }
        }
        return BigMatrix;
    }


    gsMatrix<> leftBoundarySmallSystem( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);
        size_t dimV = temp_L.size(1);

        gsMatrix<> SmallMatrix;
        SmallMatrix.setZero(dimV,6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < dimV; i++)
            {
                SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i * dimU);
            }
        }
        return SmallMatrix;
    }


    gsMatrix<> rightBoundaryBigSystem( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);

        gsMatrix<> BigMatrix;
        BigMatrix.setZero( 2 * dimU ,6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < 2 * dimU; i++)
            {
                BigMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }
        }
        return BigMatrix;
    }


    gsMatrix<> rightBoundarySmallSystem( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);

        gsMatrix<> SmallMatrix;
        SmallMatrix.setZero( dimU, 6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            for (size_t i = 0; i < dimU; i++)
            {
                SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }
        }
        return SmallMatrix;
    }


    gsMatrix<> bigInternalBoundaryPatchSystem( index_t np)
    {
        gsMultiBasis<> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(bas.basis(0));
        size_t dimU = temp_L.size(0);

        gsMatrix<> Matrix;
        Matrix.setZero( 3 ,6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            Matrix(0, bf) = auxGeom[np].getG1BasisCoefs(bf).at(0);
            Matrix(1, bf) = auxGeom[np].getG1BasisCoefs(bf).at(1);
            Matrix(2, bf) = auxGeom[np].getG1BasisCoefs(bf).at(dimU);

        }
        return Matrix;
    }

    gsMatrix<> smallInternalBoundaryPatchSystem( index_t np)
    {
        gsMatrix<> Matrix;
        Matrix.setZero( 1 ,6);

        for(size_t bf = 0; bf < 6; bf++)
        {
            Matrix(0, bf) = auxGeom[np].getG1BasisCoefs(bf).at(0);
        }
        return Matrix;
    }


    std::pair<gsMatrix<>, gsMatrix<>> createSinglePatchSystem(index_t np)
    {
        if(isBdy[np][1] == 1)
            return std::make_pair(leftBoundaryBigSystem(np), leftBoundarySmallSystem(np));
        else
        {
        if (isBdy[np][0] == 1)
            return std::make_pair(rightBoundaryBigSystem(np), rightBoundarySmallSystem(np));
        else
            return std::make_pair(bigInternalBoundaryPatchSystem(np), smallInternalBoundaryPatchSystem(np));
        }
    }


    std::pair<gsMatrix<>, std::vector<index_t>> selectVertexBoundaryBasisFunction(gsMatrix<> bigKernel, index_t bigKerDim, gsMatrix<> smallKernel, index_t smallKerDim)
    {
        gsMatrix<> basisVect;
        std::vector<index_t> numberPerType;

        numberPerType.push_back(bigKerDim); // Number of basis which has to be moved to the internal
        numberPerType.push_back(smallKerDim - bigKerDim); // Number of basis which are boundary function of FIRST TYPE
        numberPerType.push_back(6 - smallKerDim); // Number of basis which are boundary function of SECOND TYPE

        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);


        if(bigKerDim != 0)
        {
            for(index_t bk = 0; bk < bigKerDim; bk++ )
            {
                for(index_t r=0; r < 6; r++)
                {
                    if( abs(bigKernel(r, bk) )  < (10e-10) )
                        bigKernel(r, bk) = 0;
                }
            }

            for(index_t bk = 0; bk < smallKerDim; bk++ )
            {
                for(index_t r=0; r < 6; r++)
                {
                    if( abs( smallKernel(r, bk)) < (10e-10))
                        smallKernel(r, bk) = 0;
                }
            }

            basisVect = bigKernel;

            for(index_t i=0; i < smallKerDim; i++)
            {
                basisVect.conservativeResize(basisVect.rows(), basisVect.cols() + 1);
                basisVect.col(basisVect.cols()-1) = smallKernel.col(i);

                Eigen::FullPivLU<gsMatrix<>> ker(basisVect);
                if(ker.dimensionOfKernel() != 0)
                {
                    basisVect = basisVect.block(0, 0, basisVect.rows(), basisVect.cols()-1);
                }
            }

            size_t count=0;
            while (basisVect.cols() < 6)
            {
                basisVect.conservativeResize(basisVect.rows(), basisVect.cols() + 1);
                basisVect.col(basisVect.cols()-1) = vertBas.col(count);

                Eigen::FullPivLU<gsMatrix<>> ker(basisVect);
                if(ker.dimensionOfKernel() != 0)
                {
                    basisVect = basisVect.block(0, 0, basisVect.rows(), basisVect.cols()-1);
                }
                count++;
            }
        }
        else
        {
            for(index_t bk = 0; bk < smallKerDim; bk++ )
            {
                for(index_t r=0; r < 6; r++)
                {
                    if( ( abs(smallKernel(r, bk) )) < (10e-10))
                        smallKernel(r, bk) = 0;

                }

            }

            basisVect = smallKernel;

            size_t count=0;
            while (basisVect.cols() < 6)
            {
                basisVect.conservativeResize(basisVect.rows(), basisVect.cols() + 1);
                basisVect.col(basisVect.cols()-1) = vertBas.col(count);

                Eigen::FullPivLU<gsMatrix<>> ker(basisVect);
                if(ker.dimensionOfKernel() != 0)
                {
                    basisVect = basisVect.block(0, 0, basisVect.rows(), basisVect.cols()-1);
                }
                count++;
            }
        }


        gsInfo << "Big kernel:\n";
        gsInfo << bigKernel << "\n ";

        gsInfo << "Small kernel:\n";
        gsInfo << smallKernel << "\n ";

        gsInfo << "Basis:\n";
        gsInfo << basisVect << "\n";



        return std::make_pair(basisVect, numberPerType);
    }


    void computeG1InternalVertexBasis(gsOptionList optionList)
    {
        gsMultiPatch<> test_mp(this->computeAuxTopology());
        gsMultiBasis<> test_mb(test_mp);

        this->reparametrizeG1Vertex();

        this->computeSigma();

        std::string fileName;
        std::string basename = "singelFunktions";

        gsParaviewCollection collection(basename);

        std::vector<gsMultiPatch<>> g1BasisVector;

        std::pair<gsMatrix<>, std::vector<index_t>> vertexBoundaryBasis;


        for(size_t i = 0; i < auxGeom.size(); i++)
        {
            gsG1BasisVertex<real_t> g1BasisVertex_0(auxGeom[i].getPatch(),auxGeom[i].getPatch().basis(), isBdy[i], sigma, optionList);
            gsMultiPatch<> g1Basis;
            g1BasisVertex_0.constructSolution(g1Basis);
            g1BasisVector.push_back(g1Basis);
            auxGeom[i].setG1Basis(g1Basis);
        }


        if (this->kindOfVertex() == 1)
        {
            gsMatrix<> bigMatrix;
            gsMatrix<> smallMatrix;
            for (size_t i = 0; i < auxGeom.size(); i++)
            {
                std::pair<gsMatrix<>, gsMatrix<>> tmp(createSinglePatchSystem(i));
                bigMatrix.conservativeResize(bigMatrix.rows() + tmp.first.rows(), 6);
                smallMatrix.conservativeResize(smallMatrix.rows() + tmp.second.rows(), 6);

                bigMatrix.block(bigMatrix.rows(), 0, tmp.first.rows(), 6) = tmp.first;
                smallMatrix.block(smallMatrix.rows(), 0, tmp.second.rows(), 6) = tmp.second;
            }
            Eigen::FullPivLU<gsMatrix<>> BigLU(bigMatrix);
            Eigen::FullPivLU<gsMatrix<>> SmallLU(smallMatrix);

            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());

        }

        else
            if(this->kindOfVertex() == -1)
        {
            Eigen::FullPivLU<gsMatrix<>> BigLU(computeBigSystemMatrix(0));
            Eigen::FullPivLU<gsMatrix<>> SmallLU(computeSmallSystemMatrix(0));
            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());
        }




        for (size_t i = 0; i < auxGeom.size(); i++)
        {
            auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);
        }


        collection.save();
    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i){
        return auxGeom[i];
    }


protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;
    std::vector<size_t> auxVertexIndices;
    std::vector< std::vector<bool>> isBdy;
    real_t sigma;


};

}