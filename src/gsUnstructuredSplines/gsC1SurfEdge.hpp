/** @file gsApproxC1Edge.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1Edge.h>

#include <gsUnstructuredSplines/gsC1AuxiliaryPatch.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>
#include <gsUnstructuredSplines/gsApproxC1EdgeBasisProjection.h>

namespace gismo
{

    template<short_t d,class T>
    gsMultiPatch<> computeAuxTopology(){
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++){
            if(auxGeom[i].getPatch().orientation() == -1)
            {
                auxGeom[i].swapAxis();
//                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        auxTop.computeTopology();
        return auxTop;
    }


    template<short_t d,class T>
    gsMultiPatch<> reparametrizeInterface(){
        gsMultiPatch<> repTop(this->computeAuxTopology());

        if(repTop.interfaces()[0].second().side().index() == 1 && repTop.interfaces()[0].first().side().index() == 3)
            return repTop;

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (repTop.interfaces()[0].second().side().index())
        {
            case 1:
                break;
            case 4: auxGeom[0].rotateParamClock();
                break;
            case 3: auxGeom[0].rotateParamAntiClock();
                break;
            case 2: auxGeom[0].rotateParamAntiClockTwice();
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (repTop.interfaces()[0].first().side().index())
        {
            case 3:
                break;
            case 4: auxGeom[1].rotateParamAntiClockTwice();
                break;
            case 2: auxGeom[1].rotateParamAntiClock();
                break;
            case 1: auxGeom[1].rotateParamClock();
                break;
            default:
                break;
        }

        return this->computeAuxTopology();
    }


    template<short_t d,class T>
    gsMultiPatch<> reparametrizeBoundary(const int bInd){
        gsMultiPatch<> repTop(this->computeAuxTopology());
        if(auxGeom[0].getOrient())
        {
            switch (bInd)
            {
                case 3:
                    break;
                case 2:
                    auxGeom[0].rotateParamClock();
                    break;
                case 4:
                    auxGeom[0].rotateParamAntiClockTwice();
                    break;
                case 1:
                    auxGeom[0].rotateParamAntiClock();
                    break;
            }
        }
        else
        {
            switch (bInd)
            {
                case 1:
                    break;
                case 4:
                    auxGeom[0].rotateParamClock();
                    break;
                case 2:
                    auxGeom[0].rotateParamAntiClockTwice();
                    break;
                case 3:
                    auxGeom[0].rotateParamAntiClock();
                    break;
            }
        }
        return this->computeAuxTopology();
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::computeKernel(gsMultiPatch<> & result_0, gsMultiPatch<> & result_1, index_t side_0)
    {
        index_t n_plus = m_auxPatches[0].getC1BasisRotated().getBasisPlus(side_0).size(); // == patch 1, side 1
        index_t n_minus = m_auxPatches[0].getC1BasisRotated().getBasisMinus(side_0).size();

        index_t dim_U_0 = result_0.basis(0).component(0).size();
        index_t dim_V_0 = result_0.basis(0).component(1).size();

        index_t dim_U_1 = result_1.basis(0).component(0).size();
        index_t dim_V_1 = result_1.basis(0).component(1).size();

        std::vector<index_t> kernel_index;

        for (index_t topBottom = 0; topBottom < 2; ++topBottom)
        {
            index_t shift_0, shift_1;
            std::vector<index_t> bf_indices;
            if (topBottom == 0) // Bottom
            {
                bf_indices.push_back(0);
                bf_indices.push_back(1);
                bf_indices.push_back(n_plus);

                shift_0 = 0;
                shift_1 = 0;
            }
            else
            {
                bf_indices.push_back(n_plus-1);
                bf_indices.push_back(n_plus-2);
                bf_indices.push_back(n_plus+n_minus-1);

                shift_0 = dim_U_0*(dim_V_0-1);
                shift_1 = dim_U_1-1;
            }


            gsMatrix<T> coefs_interface(dim_U_0 + dim_V_1, 3);
            coefs_interface.setZero();

            // Bottom

            for (index_t i = 0; i < dim_U_0; ++i) {
                index_t j = 0;
                for (std::vector<index_t>::iterator it = bf_indices.begin(); it != bf_indices.end(); ++it, ++j) {
                    T coef_temp = result_0.patch(*it).coef(shift_0+i, 0); // v = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_interface(i, j) = coef_temp;
                }
            }
            for (index_t i = 0; i < dim_V_1; ++i) {
                index_t j = 0;
                for (std::vector<index_t>::iterator it = bf_indices.begin(); it != bf_indices.end(); ++it, ++j) {
                    T coef_temp = result_1.patch(*it).coef(i * dim_U_1 + shift_1, 0); // u = 0
                    if (coef_temp * coef_temp > 1e-25)
                        coefs_interface(dim_U_0 + i, j) = coef_temp;
                }
            }

            Eigen::FullPivLU<gsMatrix<>> KernelInterface(coefs_interface);
            KernelInterface.setThreshold(1e-5);
            //gsInfo << "Coefs: " << coefs_interface << "\n";
            //gsInfo << "Kernel: " << KernelInterface.kernel() << "\n";
            //gsInfo << "Kernel dim: " << KernelInterface.dimensionOfKernel() << "\n";

            if (KernelInterface.dimensionOfKernel() == 1) {

                kernel_index.push_back(bf_indices[0]);

                index_t numBF = bf_indices.size();

                gsMatrix<> vertBas;
                vertBas.setIdentity(numBF, numBF);

                gsMatrix<T> kernel = KernelInterface.kernel();

                size_t count = 0;
                while (kernel.cols() < numBF) {
                    kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
                    kernel.col(kernel.cols() - 1) = vertBas.col(count);

                    Eigen::FullPivLU<gsMatrix<>> ker_temp(kernel);
                    ker_temp.setThreshold(1e-5);
                    if (ker_temp.dimensionOfKernel() != 0) {
                        kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
                    }
                    count++;
                }

                gsMatrix<T> kernel_temp = kernel; // Change the kernel position to the second patch! TODO
                kernel.col(0) = kernel_temp.col(1);
                kernel.col(1) = kernel_temp.col(0);
                //gsInfo << "New basis: " << kernel << "\n";


                gsMultiPatch<> temp_result_0 = result_0;
                gsMultiPatch<> temp_result_1 = result_1;

                for (size_t j = 0; j < bf_indices.size(); ++j) {
                    gsMatrix<> coef_bf;
                    coef_bf.setZero(dim_U_0 * dim_V_0, 1);
                    index_t i = 0;
                    for (std::vector<index_t>::iterator it = bf_indices.begin(); it != bf_indices.end(); ++it, ++i)
                        if (kernel(i, j) * kernel(i, j) > 1e-25)
                            coef_bf += temp_result_0.patch(*it).coefs() * kernel(i, j);

                    result_0.patch(bf_indices[j]).setCoefs(coef_bf);

                    coef_bf.setZero(dim_U_1 * dim_V_1, 1);
                    i = 0;
                    for (std::vector<index_t>::iterator it = bf_indices.begin(); it != bf_indices.end(); ++it, ++i)
                        if (kernel(i, j) * kernel(i, j) > 1e-25)
                            coef_bf += temp_result_1.patch(*it).coefs() * kernel(i, j);

                    result_1.patch(bf_indices[j]).setCoefs(coef_bf);
                }
            }
        } //topBottom
    } // ComputeKernel


} // namespace gismo