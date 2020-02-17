//
// Created by afarahat on 2/6/20.
//

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include "gsG1Basis/gsG1AuxiliaryPatch.h"

namespace gismo
{


class gsG1AuxiliaryMultiplePatches
{

public:

    // Constructor for two patches along the common interface
    gsG1AuxiliaryMultiplePatches(const gsMultiPatch<> & mp, const size_t firstPatch, const size_t secondPatch){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), firstPatch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), secondPatch));
    }

    // Constructor for n patches around a common vertex
    gsG1AuxiliaryMultiplePatches(const gsMultiPatch<> & mp, const std::vector<unsigned> patchesAroundVertex){
        for(unsigned i = 0; i < patchesAroundVertex.size(); i++ )
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(i), patchesAroundVertex[i]));
    }


    // Compute topology 4
    gsMultiPatch<> computeAuxTopology(){
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++){
            if(auxGeom[i].getPatch().orientation() != 1)
            {
                auxGeom[i].swapAxis();
                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        // After computeTopology() the patch with initial bigger patch-index will have index zero and vice-versa
        auxTop.computeTopology();
        return auxTop;
    }

    gsMultiPatch<> reparametrizeG1Interface(){
        gsMultiPatch<> repTop(this->computeAuxTopology());
        if(repTop.interfaces()[0].second().side().index() == 3 && repTop.interfaces()[0].first().side().index() == 1)
            return repTop;

        if(repTop.interfaces()[0].second().side().index() == 1 && repTop.interfaces()[0].first().side().index() == 3)
            return repTop;

        switch (repTop.interfaces()[0].second().side().index())
        {
            case 3:
                break;
            case 1: auxGeom[0].rotateParamClock();
                break;
            case 2: auxGeom[0].rotateParamAntiClock();
                break;
            case 4: auxGeom[0].rotateParamAntiClockTwice();
                break;
        }
        switch (repTop.interfaces()[0].first().side().index())
        {
            case 1:
                break;
            case 2: auxGeom[1].rotateParamAntiClockTwice();
                break;
            case 3: auxGeom[1].rotateParamAntiClock();
                break;
            case 4: auxGeom[1].rotateParamClock();
                break;
        }
       return this->computeAuxTopology();
    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i){
        return auxGeom[i];
    }

protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;


};
}

