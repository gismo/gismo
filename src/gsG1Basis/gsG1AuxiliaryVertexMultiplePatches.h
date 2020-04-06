//
// Created by afarahat on 2/25/20.
//

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include <gsG1Basis/gsG1AuxiliaryPatch.h>

# include <gsG1Basis/gsG1BasisVertex.h>

# include <gsG1Basis/gsG1OptionList.h>

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
            checkBoundary(mp, patchesAroundVertex[i], vertexIndices[i]);
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


    void checkBoundary(const gsMultiPatch<> & mpTmp, size_t  patchInd, size_t sideInd)
    {
        std::vector<bool> tmp;
        switch (sideInd)
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
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
                    BigMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (size_t i = 1; i < dimV - 1; i++)
            {
                for(size_t j = i; j < i + 2; j++)
                {
                    if (auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i) * auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i) > m_zero*m_zero)
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
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
                    SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (size_t i = 1; i < dimV; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) > m_zero*m_zero)
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
                    if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU + j - i) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU + j - i) > m_zero*m_zero)
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
                if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) > m_zero*m_zero)
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
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
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
                if (auxGeom[np].getG1BasisCoefs(bf).at(i ) * auxGeom[np].getG1BasisCoefs(bf).at(i ) > m_zero*m_zero)
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


    void checkValues(gsMatrix<> & mat)
    {
        for(index_t bk = 0; bk < mat.cols(); bk++ )
        {
            for(index_t r=0; r < mat.rows(); r++)
            {
                if( abs( mat(r, bk) * mat(r, bk) )   < (m_zero*m_zero) )
                    mat(r, bk) = 0;
            }
        }
    }


    void addVertexBasis(gsMatrix<> & basisV)
    {
        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);
        size_t count = 0;
        while (basisV.cols() < 6)
        {
            basisV.conservativeResize(basisV.rows(), basisV.cols() + 1);
            basisV.col(basisV.cols() - 1) = vertBas.col(count);

            Eigen::FullPivLU<gsMatrix<>> ker(basisV);
            if (ker.dimensionOfKernel() != 0)
            {
                basisV = basisV.block(0, 0, basisV.rows(), basisV.cols() - 1);
            }
            count++;
        }
    }


    void addSmallKerBasis(gsMatrix<> & basisV, gsMatrix<> & smallK, index_t smallKDim)
    {
        for(index_t i=0; i < smallKDim; i++)
        {
            basisV.conservativeResize(basisV.rows(), basisV.cols() + 1);
            basisV.col(basisV.cols()-1) = smallK.col(i);

            Eigen::FullPivLU<gsMatrix<>> ker(basisV);
            if(ker.dimensionOfKernel() != 0)
            {
                basisV = basisV.block(0, 0, basisV.rows(), basisV.cols()-1);
            }
        }
    }


    std::pair<gsMatrix<>, std::vector<index_t>> selectVertexBoundaryBasisFunction(gsMatrix<> bigKernel, index_t bigKerDim, gsMatrix<> smallKernel, index_t smallKerDim)
    {
        gsMatrix<> basisVect;
        std::vector<index_t> numberPerType;

        numberPerType.push_back(bigKerDim); // Number of basis which has to be moved to the internal
        numberPerType.push_back(smallKerDim - bigKerDim); // Number of basis which are boundary function of FIRST TYPE
        numberPerType.push_back(6 - smallKerDim); // Number of basis which are boundary function of SECOND TYPE

        if(bigKerDim != 0)
        {
            checkValues(bigKernel);
            checkValues(smallKernel);

            basisVect = bigKernel;
            addSmallKerBasis(basisVect, smallKernel, smallKerDim);
            addVertexBasis(basisVect);
        }
        else
        {
            if (smallKerDim != 0)
            {
                checkValues(smallKernel);

                basisVect = smallKernel;
                addVertexBasis(basisVect);
            }
            else
            {
                gsMatrix<> vertBas;
                vertBas.setIdentity(6, 6);
                basisVect = vertBas;
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


    void computeG1InternalVertexBasis(gsG1OptionList g1OptionList)
    {
        //gsMultiPatch<> test_mp(this->computeAuxTopology());
        //gsMultiBasis<> test_mb(test_mp);
        m_zero = g1OptionList.getReal("zero");

        this->reparametrizeG1Vertex();
        this->computeSigma();

        //g1OptionList.setInt("gluingData",gluingData::global);

        std::vector<gsMultiPatch<>> g1BasisVector;
        std::pair<gsMatrix<>, std::vector<index_t>> vertexBoundaryBasis;
        std::vector<gsBSpline<>> alpha, beta_S;
        std::vector<gsG1BasisVertex<real_t>> g1BasisVertexVector;
        for(size_t i = 0; i < auxGeom.size(); i++)
        {
            gsInfo << "Index " << auxVertexIndices[i] << " Patch " << auxGeom[i].getGlobalPatchIndex() <<  "\n";

            gsG1BasisVertex<real_t> g1BasisVertex_0(auxGeom[i].getPatch(),auxGeom[i].getPatch().basis(), isBdy[i], sigma, g1OptionList);
            g1BasisVertexVector.push_back(g1BasisVertex_0);


            if (g1OptionList.getInt("gluingData")==gluingData::global)
            {
                alpha.push_back(g1BasisVertex_0.get_alpha_tilde(0));
                beta_S.push_back(g1BasisVertex_0.get_beta_tilde(0));

                alpha.push_back(g1BasisVertex_0.get_alpha_tilde(1));
                beta_S.push_back(g1BasisVertex_0.get_beta_tilde(1));
            }
            else if (g1OptionList.getInt("gluingData")==gluingData::local)
            {
                alpha.push_back(g1BasisVertex_0.get_local_alpha_tilde(0));
                beta_S.push_back(g1BasisVertex_0.get_local_beta_tilde(0));

                alpha.push_back(g1BasisVertex_0.get_local_alpha_tilde(1));
                beta_S.push_back(g1BasisVertex_0.get_local_beta_tilde(1));

            }
        }
        // COMPUTE MODIFIED TRANSVERSAL VEKTOR
        // Point zero
        gsMatrix<> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<>> dd_ik_plus, dd_ik_minus;
        dd_ik_minus.resize(auxGeom.size());
        dd_ik_plus.resize(auxGeom.size());
        if (auxGeom.size() > 1)
        {
            gsMatrix<> dd_tilde(2,1);
            dd_tilde.setZero();
            for (size_t i = 0; i < auxGeom.size(); i++)
            {
                gsMatrix<> temp_minus, temp_plus;

                temp_minus = -1/(alpha[2*i].eval(zero.row(0))(0,0)) * (auxGeom[i].getPatch().jacobian(zero).col(1) +
                    beta_S[2*i].eval(zero.row(0))(0,0) * auxGeom[i].getPatch().jacobian(zero).col(0));

                temp_plus = 1/(alpha[2*i + 1].eval(zero.row(0))(0,0)) * (auxGeom[i].getPatch().jacobian(zero).col(0) +
                    beta_S[2*i + 1].eval(zero.row(0))(0,0) * auxGeom[i].getPatch().jacobian(zero).col(1));

                if (isBdy[i][0] == false) // Does not work in case of 3 patches in a single boundary vertex
                {
                    dd_tilde += temp_minus;
                    //dd_ik_minus[i] = temp_minus;
                    dd_ik_plus[i] = temp_plus;
                }
                else
                {
                    dd_ik_minus[i] = temp_minus;
                    dd_tilde += temp_plus;
                    //dd_ik_plus[i] = temp_plus;
                }
            }

            dd_tilde /= auxGeom.size();

            for (size_t i = 0; i < auxGeom.size(); i++)
            {
                if (isBdy[i][0] == false) // Does not work in case of 3 patches in a single boundary vertex
                {
                    dd_ik_minus[i] = dd_tilde;
                }
                else
                {
                    dd_ik_plus[i] = dd_tilde;
                }
            }

        }
        else
        {
            dd_ik_minus[0] = -1 * auxGeom[0].getPatch().jacobian(zero).col(1);
            dd_ik_plus[0] = auxGeom[0].getPatch().jacobian(zero).col(0);
        }

        for (size_t i = 0; i < auxGeom.size(); i++)
        {
            gsMultiPatch<> g1Basis;
            g1BasisVertexVector[i].setG1BasisVertex(g1Basis, dd_ik_minus[i], dd_ik_plus[i]);

            g1BasisVector.push_back(g1Basis);
            auxGeom[i].setG1Basis(g1Basis);
        }


        // Plot alpha
        if (alpha.size() == 2)
        {
            std::string fileName;
            std::string basename = "GluingData";
            gsParaviewCollection collection(basename);
            for (size_t i = 0; i < alpha.size(); i++)
            {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsWriteParaview(alpha[i],fileName,5000);
                collection.addPart(fileName,"0.vts");
            }
            collection.save();
        }


        if (auxGeom.size() == 2 && g1OptionList.getInt("gluingData")==gluingData::global)
        {
            if (auxGeom[0].getGlobalPatchIndex() == 0 && isBdy[0][1])
                g1ConditionRep(alpha[3], alpha[0], beta_S[3], beta_S[0], g1BasisVector[1],  g1BasisVector[0]);
            else if (auxGeom[0].getGlobalPatchIndex() == 0 && isBdy[0][0])
                g1ConditionRep(alpha[1], alpha[2], beta_S[1], beta_S[2], g1BasisVector[0],  g1BasisVector[1]);
            else if (auxGeom[0].getGlobalPatchIndex() == 1 && isBdy[0][1])
                g1ConditionRep(alpha[3], alpha[1], beta_S[3], beta_S[1], g1BasisVector[1],  g1BasisVector[0]);
            else if (auxGeom[0].getGlobalPatchIndex() == 1 && isBdy[0][0])
                g1ConditionRep(alpha[1], alpha[2], beta_S[1], beta_S[2],  g1BasisVector[0],  g1BasisVector[1]);

        }


        if (this->kindOfVertex() == 1) // Interface-Boundary vertex
        {
            gsMatrix<> bigMatrix(0,0);
            gsMatrix<> smallMatrix(0,0);
            for (size_t i = 0; i < auxGeom.size(); i++)
            {
                std::pair<gsMatrix<>, gsMatrix<>> tmp;
                tmp = createSinglePatchSystem(i);
                size_t row_bigMatrix = bigMatrix.rows();
                size_t row_smallMatrix = smallMatrix.rows();

                bigMatrix.conservativeResize(bigMatrix.rows() + tmp.first.rows(), 6);
                smallMatrix.conservativeResize(smallMatrix.rows() + tmp.second.rows(), 6);

                bigMatrix.block(row_bigMatrix, 0, tmp.first.rows(), 6) = tmp.first;
                smallMatrix.block(row_smallMatrix, 0, tmp.second.rows(), 6) = tmp.second;
            }

            //gsInfo << smallMatrix << "\n";

            Eigen::FullPivLU<gsMatrix<>> BigLU(bigMatrix);
            Eigen::FullPivLU<gsMatrix<>> SmallLU(smallMatrix);
            SmallLU.setThreshold(g1OptionList.getReal("threshold"));
            dim_kernel = SmallLU.dimensionOfKernel();

            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());

        }
        else if(this->kindOfVertex() == -1) // Boundary vertex
        {
            Eigen::FullPivLU<gsMatrix<>> BigLU(computeBigSystemMatrix(0));
            Eigen::FullPivLU<gsMatrix<>> SmallLU(computeSmallSystemMatrix(0));
            SmallLU.setThreshold(1e-5);
            dim_kernel = SmallLU.dimensionOfKernel();

            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());
        }


        if (this->kindOfVertex() != 0)
            for (size_t i = 0; i < auxGeom.size(); i++)
            {
                gsMultiPatch<> temp_mp_g1 = g1BasisVector[i];
                for (size_t bf = 0; bf < 6; bf++)
                {
                    gsMatrix<> coef_bf;
                    coef_bf.setZero(temp_mp_g1.patch(bf).coefs().dim().first,1);
                    for (size_t lambda = 0; lambda < 6; lambda++)
                        coef_bf += temp_mp_g1.patch(lambda).coefs() * vertexBoundaryBasis.first(lambda,bf);

                    for (index_t ii = 0; ii < coef_bf.dim().first; ii++)
                        if (coef_bf.at(ii) * coef_bf.at(ii) < 1e-9)
                            coef_bf.at(ii) *= 0;

                    g1BasisVector[i].patch(bf).setCoefs(coef_bf);
                }

                auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);
            }
        else
            for (size_t i = 0; i < auxGeom.size(); i++)
                auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);

    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i){
        return auxGeom[i];
    }

    size_t get_internalDofs() { return dim_kernel; }

    void g1ConditionRep(gsBSpline<> alpha_0, gsBSpline<> alpha_1, gsBSpline<> beta_0, gsBSpline<> beta_1, gsMultiPatch<> g1Basis_0,  gsMultiPatch<> g1Basis_1)
    {
        // BETA
        // first,last,interior,mult_ends,mult_interior,degree
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getPatch().basis().component(1)); // 0 -> v, 1 -> u
        index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

        gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
        gsBSplineBasis<> bsp(kv);

        gsMatrix<> greville = bsp.anchors();
        gsMatrix<> uv1, uv0, ev1, ev0;

        const index_t d = 2;
        gsMatrix<> D0(d,d);

        gsGeometry<>::Ptr beta_temp;

        uv0.setZero(2,greville.cols());
        uv0.bottomRows(1) = greville;

        uv1.setZero(2,greville.cols());
        uv1.topRows(1) = greville;

        const gsGeometry<> & P0 = auxGeom[0].getPatch(); // iFace.first().patch = 1
        const gsGeometry<> & P1 = auxGeom[1].getPatch(); // iFace.second().patch = 0
        // ======================================

        // ======== Determine bar{beta} ========
        for(index_t i = 0; i < uv1.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i),ev0);
            P1.jacobian_into(uv1.col(i),ev1);

            D0.col(1) = ev0.col(0); // (DuFL, *)
            D0.col(0) = ev1.col(1); // (*,DuFR)

            uv0(0,i) = D0.determinant();
        }

        beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
        gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);


        index_t p_size = 8;
        gsMatrix<> points(1, p_size);
        points.setRandom();
        points = points.array().abs();

        gsVector<> vec;
        vec.setLinSpaced(p_size,0,1);
        points = vec.transpose();

        gsMatrix<> points2d_0(2, p_size);
        gsMatrix<> points2d_1(2, p_size);

        points2d_0.setZero();
        points2d_1.setZero();
        points2d_0.row(1) = points; // v
        points2d_1.row(0) = points; // u

        gsVector<> g1Error(6);
        g1Error.setZero();

        gsMatrix<> temp;
        temp = alpha_1.eval(points).cwiseProduct(beta_0.eval(points))
            + alpha_0.eval(points).cwiseProduct(beta_1.eval(points))
            - beta.eval(points);

        gsInfo << "alpha : " << alpha_0.eval(points) << "\n";
        gsInfo << "alpha : " << alpha_1.eval(points) << "\n";
        gsInfo << "alpha : " << beta_0.eval(points) << "\n";;
        gsInfo << "alpha : " << beta_1.eval(points) << "\n";;
        gsInfo << "alpha : " << beta.eval(points) << "\n";

        gsInfo << "Conditiontest Gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";

        for (size_t i = 0; i < g1Basis_0.nPatches(); i++)
        {
            gsMatrix<> temp;
            temp = alpha_1.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).topRows(1))
                + alpha_0.eval(points).cwiseProduct(g1Basis_1.patch(i).deriv(points2d_1).bottomRows(1))
                + beta.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).bottomRows(1));

            if (temp.array().abs().maxCoeff() > g1Error[i])
                g1Error[i] = temp.array().abs().maxCoeff();
        }

        gsInfo << "Conditiontest G1 continuity VERTEX: \n" << g1Error << "\n\n";
    }



protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;
    std::vector<size_t> auxVertexIndices;
    std::vector< std::vector<bool>> isBdy;
    real_t sigma;
    size_t dim_kernel;

    real_t m_zero;


};

}