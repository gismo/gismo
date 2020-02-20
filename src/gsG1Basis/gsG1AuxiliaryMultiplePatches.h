/** @file gsG1AuxiliaryMultiplePatches.h
 *
    @brief Reparametrize the Geometry for one Interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include <gsG1Basis/gsG1AuxiliaryPatch.h>

namespace gismo
{


class gsG1AuxiliaryMultiplePatches
{

public:

    // Constructor for two patches along the common interface
    gsG1AuxiliaryMultiplePatches(const gsMultiPatch<> & mp, const size_t firstPatch, const size_t secondPatch){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), firstPatch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), secondPatch));
        for(unsigned i = 0; i <  auxGeom.size(); i++){
            if(auxGeom[i].getPatch().orientation() == -1)
            {
                auxGeom[i].swapAxis();
                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
        }
    }

    // Constructor for n patches around a common vertex
    gsG1AuxiliaryMultiplePatches(const gsMultiPatch<> & mp, const std::vector<size_t> patchesAroundVertex, std::vector<size_t> vertexIndices){
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(patchesAroundVertex[i]), patchesAroundVertex[i]));
            if(auxGeom[i].getPatch().orientation() == -1)
            {
                auxGeom[i].swapAxis();
                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
            this->reparametrizeG1Vertex(i, vertexIndices[i]);
        }
        gsInfo << "\n";
    }


    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index inside auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4

    gsMultiPatch<> computeAuxTopology(){
        
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++){
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        // After computeTopology() the patch with initial bigger patch-index will have index zero and vice-versa
        auxTop.computeTopology();
        return auxTop;
    }

    gsMultiPatch<> reparametrizeG1Interface(){
        gsMultiPatch<> repTop(this->computeAuxTopology());
        if(repTop.interfaces()[0].second().side().index() == 1 && repTop.interfaces()[0].first().side().index() == 3)
            return repTop;

        switch (repTop.interfaces()[0].second().side().index())
        {
            case 1:
                gsInfo << "Patch: " << repTop.interfaces()[0].second().patch << " not rotated\n";
                break;
            case 4: auxGeom[0].rotateParamClock();
                gsInfo << "Patch: " << repTop.interfaces()[0].second().patch << " rotated clockwise\n";
                break;
            case 3: auxGeom[0].rotateParamAntiClock();
                gsInfo << "Patch: " << repTop.interfaces()[0].second().patch << " rotated anticlockwise\n";
                break;
            case 2: auxGeom[0].rotateParamAntiClockTwice();
                gsInfo << "Patch: " << repTop.interfaces()[0].second().patch << " rotated twice anticlockwise\n";
                break;
            default:
                break;
        }
        switch (repTop.interfaces()[0].first().side().index())
        {
            case 3:
                gsInfo << "Patch: " << repTop.interfaces()[0].first().patch << " not rotated\n";
                break;
            case 4: auxGeom[1].rotateParamAntiClockTwice();
                gsInfo << "Patch: " << repTop.interfaces()[0].first().patch << " rotated twice anticlockwise\n";
                break;
            case 2: auxGeom[1].rotateParamAntiClock();
                gsInfo << "Patch: " << repTop.interfaces()[0].first().patch << " rotated anticlockwise\n";
                break;
            case 1: auxGeom[1].rotateParamClock();
                gsInfo << "Patch: " << repTop.interfaces()[0].first().patch << " rotated clockwise\n";
                break;
            default:
                break;
        }
       return this->computeAuxTopology();
    }

    void parametrizeBack(){

    }




    void reparametrizeG1Vertex(size_t patchInd, size_t vertexIndex){
        if(auxGeom[patchInd].getOrient() == 0)
        {
            switch (vertexIndex)
            {
                case 1:
                    gsInfo << "Patch: " << patchInd << " not rotated\n";
                    break;
                case 4:
                    auxGeom[patchInd].rotateParamAntiClockTwice();
                    gsInfo << "Patch: " << patchInd << " rotated twice anticlockwise\n";
                    break;
                case 2:
                    auxGeom[patchInd].rotateParamAntiClock();
                    gsInfo << "Patch: " << patchInd << " rotated anticlockwise\n";
                    break;
                case 3:
                    auxGeom[patchInd].rotateParamClock();
                    gsInfo << "Patch: " << patchInd << " rotated clockwise\n";
                    break;
            }
        }
        else{
            switch (vertexIndex)
            {
                case 1:
                    gsInfo << "Patch: " << patchInd << " not rotated\n";
                    break;
                case 4:
                    auxGeom[patchInd].rotateParamAntiClockTwice();
                    gsInfo << "Patch: " << patchInd << " rotated twice anticlockwise\n";
                    break;
                case 3:
                    auxGeom[patchInd].rotateParamAntiClock();
                    gsInfo << "Patch: " << patchInd << " rotated anticlockwise\n";
                    break;
                case 2:
                    auxGeom[patchInd].rotateParamClock();
                    gsInfo << "Patch: " << patchInd << " rotated clockwise\n";
                    break;
            }
        }
    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i){
        return auxGeom[i];
    }

protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;


};
}

