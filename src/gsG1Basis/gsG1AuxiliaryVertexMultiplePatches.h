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


    void computeG1InternalVertexBasis(gsOptionList optionList){

        gsMultiPatch<> test_mp(this->computeAuxTopology());
        gsMultiBasis<> test_mb(test_mp);

        this->reparametrizeG1Vertex();

        this->computeSigma();

        std::string fileName;
        std::string basename = "singelFunktions";

        gsParaviewCollection collection(basename);


        for(size_t i = 0; i < auxGeom.size(); i++)
        {


            gsG1BasisVertex<real_t> g1BasisVertex_0(auxGeom[i].getPatch(),auxGeom[i].getPatch().basis(), isBdy[i], sigma, optionList);
            gsMultiPatch<> g1Basis;
            g1BasisVertex_0.constructSolution(g1Basis);
            //g1BasisVertex_0.plotG1BasisBoundary(g1Basis, auxGeom[i].getPatch(),"BasisVertex0");





            auxGeom[i].parametrizeBasisBack(g1Basis);

            fileName = basename + "_" + util::to_string(i);
            gsField<> temp_field(auxGeom[i].getPatch(),g1Basis.patch(0));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addPart(fileName,"0.vts");


            //g1BasisVertex_0.plotG1BasisBoundary(auxGeom[i].getG1Basis(), test_mp.patch(0),"BasisVertex_new");
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