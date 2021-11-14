/** @file gsG1AuxiliaryPatch.h
 *
    @brief Reparametrize one Patch

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/
#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>


namespace gismo
{


    class gsG1AuxiliaryPatch
    {
    public:

        gsG1AuxiliaryPatch()
        {}

        gsG1AuxiliaryPatch(const gsGeometry<> & singlePatch, const size_t globalPatchIndex):
                auxPatch(singlePatch), patchIndex(globalPatchIndex){
            rotationNum = 0;
            axisOrientation = 0;
//        gsInfo << "Single patch created: " << patchIndex << "\n";
        };

        void setPlusMinus(index_t plus, index_t minus)
        {
            m_plus = plus;
            m_minus = minus;
        }

        size_t get_n_plus() { return m_plus; }
        size_t get_n_minus() { return m_minus; }

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

            auxPatch.swap(newpatch);
            auxPatch.computeTopology();

            // Update the number of rotation of the axis
            rotationNum++;
        }

        void rotateBasisAntiClock(){
            gsMultiPatch<> newpatch;
            for(size_t np = 0; np < G1repBasis.nPatches(); np++)
            {
                gsMultiBasis<> auxBase(G1repBasis.patch(np));
                gsTensorBSplineBasis<2, real_t>
                        & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
                gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
                gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
                index_t dimU = temp_L.size(0);
                index_t dimV = temp_L.size(1);

                // The number of cols has to match the dimension of the space
                gsMatrix<> mpar(dimU * dimV, G1repBasis.patch(np).targetDim());

                // Loop over the cols
                for (index_t j = 0; j < dimU; j++)
                {   // Loop over the rows
                    for (index_t i = 0; i < dimV; i++)
                    {
                        mpar.row(i + j * dimV) = G1repBasis.patch(np).coefs().row((dimU - 1 - j) + dimU * i);
                    }
                }

                // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
                gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

                newpatch.addPatch(newgeom1);
            }
            G1repBasis.swap(newpatch);
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

            auxPatch.swap(newpatch);
            auxPatch.computeTopology();

            // Update the number of rotation of the axis
            rotationNum--;
            this->checkRotation();
        }


        void rotateBasisClock(){
            gsMultiPatch<> newpatch;
            for(size_t np = 0; np < G1repBasis.nPatches(); np++)
            {
                gsMultiBasis<> auxBase(G1repBasis.patch(np));
                gsTensorBSplineBasis<2, real_t>
                        & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
                gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
                gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
                index_t dimU = temp_L.size(0);
                index_t dimV = temp_L.size(1);

                //gsInfo << " dimU: " << dimU << "\t dimV: " << dimV << "\n";

                // The number of cols has to match the dimension of the space
                gsMatrix<> mpar(dimU * dimV, G1repBasis.patch(np).targetDim());

                for (index_t j = (dimU - 1); j >= 0; j--)
                {   // Loop over the rows
                    for (index_t i = 0; i < dimV; i++)
                    {
                        mpar.row(i + (dimU - j - 1) * dimV) =
                                G1repBasis.patch(np).coefs().row((dimV * dimU - 1 - j) - dimU * i);
                    }
                }

                // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
                gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

                newpatch.addPatch(newgeom1);
            }
            G1repBasis = newpatch;
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

            for (index_t i = 0; i < ( dimU * dimV ); i++)
            {
                mpar.row(i) = auxPatch.patch(0).coefs().row((dimU * dimV - 1) - i);
            }
            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<2, real_t> newgeom1(temp_basisLU.knots(), temp_basisLV.knots(), mpar);

            // Create a new single-patch object
            gsMultiPatch<> newpatch;

            newpatch.addPatch(newgeom1);

            auxPatch.swap(newpatch);
            auxPatch.computeTopology();

            // Update the number of rotation of the axis (anti-clockwise)
            rotationNum+=2;
            this->checkRotation();
        }


        void rotateBasisAntiClockTwice(){
            gsMultiPatch<> newpatch;
            for(size_t np = 0; np < G1repBasis.nPatches(); np++)
            {
                gsMultiBasis<> auxBase(G1repBasis.patch(np));
                gsTensorBSplineBasis<2, real_t>
                        & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
                gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
                gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
                index_t dimU = temp_L.size(0);
                index_t dimV = temp_L.size(1);

                // The number of cols has to match the dimension of the space
                gsMatrix<> mpar(dimU * dimV, G1repBasis.patch(np).targetDim());

                for (index_t i = 0; i < (dimU * dimV); i++)
                {
                    mpar.row(i) = G1repBasis.patch(np).coefs().row((dimU * dimV - 1) - i);
                }
                // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
                gsTensorBSpline<2, real_t> newgeom1(temp_basisLU.knots(), temp_basisLV.knots(), mpar);

                newpatch.addPatch(newgeom1);
            }
            G1repBasis = newpatch;
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

            auxPatch.swap(newpatch);
            //auxPatch.computeTopology();

            axisOrientation = 1;
        }


        void swapBasisAxis(){
            gsMultiPatch<> newpatch;
            for(size_t np = 0; np < G1repBasis.nPatches(); np++)
            {
                gsMultiBasis<> auxBase(G1repBasis.patch(np));
                gsTensorBSplineBasis<2, real_t> & temp_L = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(auxBase.basis(0));
                gsBSplineBasis<> temp_basisLU = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(0));
                gsBSplineBasis<> temp_basisLV = dynamic_cast<gsBSplineBasis<> &>(temp_L.component(1));
                index_t dimU = temp_L.size(0);
                index_t dimV = temp_L.size(1);

                // The number of cols has to match the dimension of the space
                gsMatrix<> mpar(dimU * dimV, G1repBasis.patch(np).targetDim());

                for (index_t j = 0; j < dimU; j++)
                {   // Loop over the rows
                    for (index_t i = 0; i < dimV; i++)
                    {
                        mpar.row(i + j * dimV) = G1repBasis.patch(np).coefs().row(j + dimU * i);
                    }
                }

                // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
                gsTensorBSpline<2, real_t> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

                newpatch.addPatch(newgeom1);
            }
            G1repBasis.swap(newpatch);

        }


        void parametrizeBasisBack( gsMultiPatch<>  g1Basis){
            G1repBasis = g1Basis;

            //gsInfo << "Patch " << patchIndex << " old: " << G1repBasis.patch(0).coefs()<< "\n";

            switch (rotationNum)
            {
                case 2:
                    this->rotateBasisAntiClockTwice();
                    break;
                case -1:
                    this->rotateBasisAntiClock();
                    break;
                case 1:
                    this->rotateBasisClock();
                    break;
                case 0:
                    break;
                default:
                    break;
            }

            if(axisOrientation)
                this->swapBasisAxis();

            //gsInfo << "Patch " << patchIndex << " new: " << G1repBasis.patch(0).coefs() << "\n";
        }

        void setG1Basis(gsMultiPatch<>  g1Ba)
        {
            G1repBasis = g1Ba;
        }

        void computeTopology(){
            this->auxPatch.computeTopology();
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

        const index_t getOrient(){
            return axisOrientation;
        }

        gsMatrix<> getG1BasisCoefs(index_t i){
            return G1repBasis.patch(i).coefs();
        }

        gsMultiPatch<> & getG1Basis(){
            return G1repBasis;
        }

        void checkRotation(){
            if(rotationNum == 4)
                rotationNum = 0;
        }

        void checkOrientation(){
            axisOrientation = ( axisOrientation == 0 ? 1 : 0 );
        }



    protected:

        gsMultiPatch<> auxPatch;
        gsMultiPatch<> G1repBasis;

        // Global patch index in the initial geometry
        index_t patchIndex;
        // Stores the changing of the axis
        // 0 -> axis not changed
        // 1 -> axis swapped (x, y --> y, x)
        bool axisOrientation;
        // How many rotation of the axis has been executed
        // Positive -> anticlockwise    Negative -> clockwise
        index_t rotationNum;

        index_t m_plus, m_minus;

    };

}