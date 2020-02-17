//
// Created by afarahat on 2/10/20.
//
#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{


class gsG1AuxiliaryPatch
{
public:

    gsG1AuxiliaryPatch(const gsGeometry<> & singlePatch, const size_t globalPatchIndex):
    auxPatch(singlePatch), patchIndex(globalPatchIndex){
        rotationNum = 0;
        axisOrientation = 0;
        gsInfo << "Single patch created: " << patchIndex << "\n";
    };



    void rotateParamAntiClock(){
        gsMultiBasis<> auxBase(auxPatch);
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
        gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
        gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<> mpar(dimU * dimV, auxPatch.targetDim());

        // Loop over the cols
        for (index_t j = 0; j < dimU; j++)
        {   // Loop over the rows
            for (index_t i = 0; i < dimV; i++)
            {
                mpar.row(i + j * dimV) = auxPatch.patch(0).coefs().row((dimU - 1 - j) + dimU * i);
            }
        }

        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<> newpatch;

        newpatch.addPatch(newgeom1);
        newpatch.computeTopology();

        auxPatch.swap(newpatch);

        // Update the number of rotation of the axis
        rotationNum++;
        this->checkRotation();
    }

    void rotateParamClock(){
        gsMultiBasis<> auxBase(auxPatch);
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
        gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
        gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<> mpar(dimU * dimV, auxPatch.targetDim());

        for (index_t j = (dimU - 1 ) ; j >= 0; j--)
        {   // Loop over the rows
            for (index_t i = 0; i < dimV; i++)
            {
                mpar.row(i + (dimU - j - 1) * dimV) = auxPatch.patch(0).coefs().row((dimV * dimU  -1 - j) - dimU * i);
            }
        }

        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<> newpatch;

        newpatch.addPatch(newgeom1);
        newpatch.computeTopology();

        auxPatch.swap(newpatch);

        // Update the number of rotation of the axis
        rotationNum--;
        this->checkRotation();
    }

    void rotateParamAntiClockTwice(){
        gsMultiBasis<> auxBase(auxPatch);
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
        gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
        gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<> mpar(dimU * dimV, auxPatch.targetDim());

        for (index_t i = 0; i < (dimU * dimV - 1); i++)
        {
            mpar.row(i) = auxPatch.patch(0).coefs().row((dimU * dimV - 1) - i);
        }
        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<2, real_t> newgeom1(temp_basisLU.knots(), temp_basisLV.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<> newpatch;

        newpatch.addPatch(newgeom1);
        newpatch.computeTopology();

        auxPatch.swap(newpatch);

        // Update the number of rotation of the axis (anti-clockwise)
        rotationNum+=2;
        this->checkRotation();
    }

    void swapAxis(){
        gsMultiBasis<> auxBase(auxPatch);
        gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
        gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
        gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<> mpar(dimU * dimV, auxPatch.targetDim());

        for (index_t j = 0; j < dimU; j++)
        {   // Loop over the rows
            for (index_t i = 0; i < dimV; i++)
            {
                mpar.row(i + j * dimV) = auxPatch.patch(0).coefs().row(j + dimU * i);
            }
        }

        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<> newpatch;

        newpatch.addPatch(newgeom1);
        newpatch.computeTopology();

        auxPatch.swap(newpatch);

        axisOrientation = 1;
    }


    gsGeometry<>& getPatch(){
        return auxPatch.patch(0);
    }

    const index_t getGlobalPatchIndex(){
        return patchIndex;
    }

    const index_t getNumberOfRotatioin(){
        return rotationNum;
    }

    void checkRotation(){
        if(rotationNum == 4)
            rotationNum = 0;
    }

protected:

    gsMultiPatch<> auxPatch;

    index_t patchIndex;
    // Stores the orientation of the axis(0 -> x , y while 1 -> -x , y / x , -y)
    index_t axisOrientation;
    // How many rotation of the axis has been executed
    index_t rotationNum;

};

}