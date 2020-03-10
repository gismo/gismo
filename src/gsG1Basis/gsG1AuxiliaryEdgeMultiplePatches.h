/** @file gsG1AuxiliaryEdgeMultiplePatches.h
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
# include <gsG1Basis/gsG1BasisEdge.h>


namespace gismo
{


class gsG1AuxiliaryEdgeMultiplePatches
{

public:

    // Constructor for one patch and it's boundary
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & sp, const size_t patchInd){
        auxGeom.push_back(gsG1AuxiliaryPatch(sp.patch(patchInd), patchInd));
    }

    // Constructor for two patches along the common interface
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & mp, const size_t firstPatch, const size_t secondPatch){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), firstPatch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), secondPatch));
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


    gsMultiPatch<> reparametrizeG1Interface(){
        gsMultiPatch<> repTop(this->computeAuxTopology());

        if(repTop.interfaces()[0].second().side().index() == 1 && repTop.interfaces()[0].first().side().index() == 3)
            return repTop;

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (repTop.interfaces()[0].second().side().index())
        {
            case 1:
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << "\tLocal patch: " << repTop.interfaces()[0].second().patch << " not rotated\n";
                break;
            case 4: auxGeom[0].rotateParamClock();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated clockwise\n";
                break;
            case 3: auxGeom[0].rotateParamAntiClock();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated anticlockwise\n";
                break;
            case 2: auxGeom[0].rotateParamAntiClockTwice();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated twice anticlockwise\n";
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (repTop.interfaces()[0].first().side().index())
        {
            case 3:
                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " not rotated\n";
                break;
            case 4: auxGeom[1].rotateParamAntiClockTwice();
                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated twice anticlockwise\n";
                break;
            case 2: auxGeom[1].rotateParamAntiClock();
                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated anticlockwise\n";
                break;
            case 1: auxGeom[1].rotateParamClock();
                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated clockwise\n";
                break;
            default:
                break;
        }
       return this->computeAuxTopology();
    }


    gsMultiPatch<> reparametrizeG1Boundary(const int bInd){
        gsMultiPatch<> repTop(this->computeAuxTopology());
        switch (bInd)
        {
            case 1:
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " not rotated\n";
                break;
            case 4: auxGeom[0].rotateParamClock();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated clockwise\n";
                break;
            case 2: auxGeom[0].rotateParamAntiClockTwice();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated twice anticlockwise\n";
                break;
            case 3: auxGeom[0].rotateParamAntiClock();
                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated anticlockwise\n";
                break;
        }
        return this->computeAuxTopology();
    }


    void computeG1InterfaceBasis(gsOptionList optionList){

        gsMultiPatch<> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());// Right -> 0 = v along the interface
        mp_init.addPatch(auxGeom[1].getPatch()); // Left -> 1 = u along the interface

        gsMultiPatch<> test_mp(this->reparametrizeG1Interface()); // auxGeom contains now the reparametrized geometry
        gsMultiBasis<> test_mb(test_mp);

//      gsInfo << "p_tilde : " << optionList << "\n";
        gsG1BasisEdge<real_t> g1BasisEdge_0(test_mp.patch(0), test_mb.basis(0), 1, false, optionList);
        gsG1BasisEdge<real_t> g1BasisEdge_1(test_mp.patch(1), test_mb.basis(1), 0, false, optionList);
        gsMultiPatch<> g1Basis_0, g1Basis_1, g1Basis_edge;
        g1BasisEdge_0.constructSolution(g1Basis_0);
        g1BasisEdge_1.constructSolution(g1Basis_1);
        if (optionList.getSwitch("plot"))
            g1BasisEdge_0.plotG1Basis(g1Basis_0,g1Basis_1, test_mp, "G1Basis_old");

//      Patch 0 -> Right
        auxGeom[0].parametrizeBasisBack(g1Basis_0);

//      Patch 1 -> Left
        auxGeom[1].parametrizeBasisBack(g1Basis_1);

        if (optionList.getSwitch("plot"))
            g1BasisEdge_0.plotG1Basis(auxGeom[0].getG1Basis(),auxGeom[1].getG1Basis(), mp_init, "G1Basis");
        //g1BasisEdge_0.g1Condition();
    }


    void computeG1BoundaryBasis(gsOptionList optionList, const int boundaryInd){
        //gsMultiPatch<> mp_init;
        //mp_init.addPatch(auxGeom[0].getPatch());

        gsMultiPatch<> test_mp(this->reparametrizeG1Boundary(boundaryInd));
        gsMultiBasis<> test_mb(test_mp);

        gsG1BasisEdge<real_t> g1BasisEdge(test_mp, test_mb, 1, true, optionList);
        gsMultiPatch<> g1Basis_edge;
        g1BasisEdge.constructSolution(g1Basis_edge);
        //std::string basename_old = "G1BasisBoundary_Old_" + util::to_string(auxGeom[0].getGlobalPatchIndex()) + "_" + util::to_string(boundaryInd);
        //g1BasisEdge.plotG1BasisBoundary(g1Basis_edge, mp_init, basename_old);

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);

        //std::string basename = "G1BasisBoundary_" + util::to_string(auxGeom[0].getGlobalPatchIndex()) + "_" + util::to_string(boundaryInd);
        //if (optionList.getSwitch("plot"))
        //    g1BasisEdge.plotG1BasisBoundary(auxGeom[0].getG1Basis(), mp_init, basename);

    }

    void computeG1EdgeBasis(gsOptionList optionList, const int edgeInd, bool isboundary){
        gsMultiPatch<> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());

        gsMultiPatch<> test_mp;
        test_mp = this->reparametrizeG1Boundary(edgeInd);
        gsMultiBasis<> test_mb(test_mp);

        gsG1BasisEdge<real_t> g1BasisEdge(test_mp, test_mb, 1, isboundary, optionList);
        gsMultiPatch<> g1Basis_edge;
        g1BasisEdge.constructSolution(g1Basis_edge);
        //std::string basename_old = "G1BasisBoundary_Old_" + util::to_string(auxGeom[0].getGlobalPatchIndex()) + "_" + util::to_string(boundaryInd);
        //g1BasisEdge.plotG1BasisBoundary(g1Basis_edge, mp_init, basename_old);

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);
        auxGeom[0].setPlusMinus(g1BasisEdge.get_n_plus(), g1BasisEdge.get_n_minus());

        //std::string basename = "G1BasisBoundary_" + util::to_string(auxGeom[0].getGlobalPatchIndex()) + "_" + util::to_string(boundaryInd);
        //if (optionList.getSwitch("plot"))
        //   g1BasisEdge.plotG1BasisBoundary(auxGeom[0].getG1Basis(), mp_init, basename);

    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i){
        return auxGeom[i];
    }

protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;
};
}

