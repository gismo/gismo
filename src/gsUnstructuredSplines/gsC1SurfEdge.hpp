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

#include <gsUnstructuredSplines/gsG1AuxiliaryPatch.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>
#include <gsUnstructuredSplines/gsApproxC1EdgeBasisProjection.h>

namespace gismo
{

    template<short_t d,class T>
    gsMultiPatch<> gsC1SurfEdge<d,T>::computeAuxTopology(){
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
    gsMultiPatch<> gsC1SurfEdge<d,T>::reparametrizeInterface(){
        gsMultiPatch<> repTop(computeAuxTopology());

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
    gsMultiPatch<> gsC1SurfEdge<d,T>::reparametrizeBoundary(const int bInd){
        gsMultiPatch<> repTop(computeAuxTopology());
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

} // namespace gismo