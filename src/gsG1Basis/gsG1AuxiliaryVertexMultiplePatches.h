//
// Created by afarahat on 2/25/20.
//

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include <gsG1Basis/gsG1AuxiliaryPatch.h>
# include <gsG1Basis/gsG1BasisEdge.h>


namespace gismo
{

class gsG1AuxiliaryVertexMultiplePatches
{

public:

// Constructor for n patches around a common vertex
    gsG1AuxiliaryVertexMultiplePatches(const gsMultiPatch<> & mp, const std::vector<size_t> patchesAroundVertex, std::vector<size_t> vertexIndices){
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
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    gsMultiPatch<> computeAuxTopology(){
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++){
            if(auxGeom[i].getPatch().orientation() == -1)
            {
                auxGeom[i].swapAxis();
                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        auxTop.computeTopology();
        return auxTop;
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
            // If the orientation is changed
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