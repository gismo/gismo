//
// Created by afarahat on 3/5/20.
//

#pragma once

#include <gismo.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>


namespace gismo
{

class gsG1Mapper
{

public:
    gsG1Mapper( const gsMultiPatch<> & geo )
    {
        std::vector<index_t> aux(geo.nPatches() * 4, 0);
        reducedEdgeMap = aux;
        reducedVertexMap = aux;

        gsMultiBasis<> basis(geo);
        for(size_t i = 0; i < basis.nBases();i++)
        {
            basisPerPatch.push_back(basis.basis(i).size());
            gsTensorBSplineBasis<2, real_t> & temp = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis.basis(i));
            dimensionPerBasis.push_back(std::pair<index_t, index_t> (temp.size(0), temp.size(1)));

            for(index_t j = 0; j < basis.basis(i).size(); j++ )
                globalMap.push_back(i * basis.basis(i).size() + j);

            this->globalEdgeMapper(i);
            this->globalVertexMapper(i);
        }

        this->internalPatchMapper();

        this->reducedEdgeMapper(geo);

        this->reducedVertexmapper(geo);

    }


    gsG1Mapper( const gsMultiPatch<> & geo, const std::vector<size_t> numG1B, std::vector<index_t> nP )
    {


        std::vector<index_t> aux(numG1B[numG1B.size() - 1], 0);
        reducedG1BasisEdgeMap = aux;
        numG1EdgeBasis = numG1B;
        nPlusEdgeDimension = nP;


        this->reducedBasisEdgeMapper(geo);

    }


    index_t localToGlobalPatchMapper(index_t patchInd, index_t i){
        index_t index = i;
        for(index_t p = 0; p < patchInd;p++)
            index += basisPerPatch[p];
        return index;
    }

    std::pair<index_t, index_t> globalToLocalPatchMapper(index_t globalIndex){
        std::pair<index_t, index_t> patchAndInd(0, 0);

            for(size_t i = 0; i < basisPerPatch.size(); i++)
            {
                if( globalIndex - basisPerPatch[i] < 0)
                {
                    patchAndInd.first = i;
                    patchAndInd.second += globalIndex;
                    return patchAndInd;
                }
                else
                {
                    patchAndInd.first = i;
                    patchAndInd.second += basisPerPatch[i];
                }
            }
    }


    void internalPatchMapper()
    {
        for(size_t np = 0; np < basisPerPatch.size(); np++)
        {
            for(index_t j = 2; j < dimensionPerBasis[np].second - 2; j++)
            {
                for(index_t i = 2; i < dimensionPerBasis[np].first - 2; i++)
                {
                    internalPatchMap.push_back(globalMap[j * dimensionPerBasis[np].first + i]);
                }
            }
        }
    }


    void globalEdgeMapper(size_t patchInd)
    {
        globalEdgeMap.push_back(patchInd * 4 + 1);
        globalEdgeMap.push_back(patchInd * 4 + 2);
        globalEdgeMap.push_back(patchInd * 4 + 3);
        globalEdgeMap.push_back(patchInd * 4 + 4);
    }


    void reducedEdgeMapper(const gsMultiPatch<> & geo)
    {
        index_t count = 1;
        for(size_t i = 0; i < reducedEdgeMap.size(); i++)
        {
            if(reducedEdgeMap[i] == 0)
            {
                for(const boundaryInterface &  item : geo.interfaces() )
                {
                    if( ( item.second().patch * 4 + item.second().side() -1) == i  )
                    {
                        std::vector<index_t> temp;
                        reducedEdgeMap[i] = count;
                        reducedEdgeMap[item.first().patch * 4 + item.first().side() - 1] = count;
                        temp.push_back( item.second().patch * 4 + item.second().side() );
                        temp.push_back( item.first().patch * 4 + item.first().side() );
                        interfaceEdgeMap.push_back(temp);

                        count++;
                    }
                }
                if(reducedEdgeMap[i] == 0)
                {
                    reducedEdgeMap[i] = count;
                    reducedBoundaryEdgeMap.push_back(count);
                    count++;
                }
            }
        }
    }

    void printReducedEdgeMapper()
    {
        for(size_t i = 0; i < reducedEdgeMap.size(); i++)
            gsInfo << "Reduced edge mapper: " << reducedEdgeMap[i] << "\n";
    }


    void printInterfaceEdgeMapper()
    {
        for(size_t i = 0; i < interfaceEdgeMap.size(); i++)
            gsInfo << "Interface edge mapper: " << interfaceEdgeMap[i][0] << ", " << interfaceEdgeMap[i][1] << "\n";
    }


    void printReducedBoundaryEdgeMapper()
    {
        for(size_t i = 0; i < reducedBoundaryEdgeMap.size(); i++)
            gsInfo << "Reduced boundary edge mapper: " << reducedBoundaryEdgeMap[i] << "\n";
    }


    index_t localToGlobalEdgeMapper(index_t patchInd, index_t i)
    {
        return reducedEdgeMap[patchInd * 4 + i - 1];
    }



    void globalVertexMapper(size_t patchInd)
    {
        globalVertexMap.push_back(patchInd * 4 + 1);
        globalVertexMap.push_back(patchInd * 4 + 2);
        globalVertexMap.push_back(patchInd * 4 + 3);
        globalVertexMap.push_back(patchInd * 4 + 4);
    }


    void reducedVertexmapper(const gsMultiPatch<> & geo)
    {
        index_t countRed = 1;

        std::vector<std::vector<patchCorner>> allcornerLists = geo.vertices();
        for(size_t i=0; i < allcornerLists.size(); i++)
        {
            gsMultiPatch<> aux;
            std::vector<size_t> patchIndex;
            std::vector<size_t> vertIndex;
            for(size_t j = 0; j < allcornerLists[i].size(); j++)
            {
                patchIndex.push_back(allcornerLists[i][j].patch);
                vertIndex.push_back(allcornerLists[i][j].m_index);
                aux.addPatch(geo.patch(allcornerLists[i][j].patch));
            }
            aux.computeTopology();

            if(patchIndex.size() == 1)
            {
                reducedVertexMap[patchIndex[0] * 4 + vertIndex[0] - 1] = countRed;
                reducedBoundaryVertexMap.push_back( countRed );
                reducedAllBoundaryVertexMap.push_back(countRed);
                countRed++;
            }
            else
            {
                if(patchIndex.size() == aux.interfaces().size())
                {   std::vector<index_t> auxInt;
                    for(size_t np = 0; np < patchIndex.size(); np++)
                    {
                        reducedVertexMap[patchIndex[np] * 4 + vertIndex[np] - 1] = countRed;
                        auxInt.push_back( patchIndex[np] * 4 + vertIndex[np] );
                    }
                    internalVertexMap.push_back(auxInt);
                }
                else
                {
                    std::vector<index_t> auxBdy;
                    for(size_t np = 0; np < patchIndex.size(); np++)
                    {
                        reducedVertexMap[patchIndex[np] * 4 + vertIndex[np] - 1] = countRed;
                        auxBdy.push_back( patchIndex[np] * 4 + vertIndex[np] );
                    }
                    boundaryInterfaceVertexMap.push_back(auxBdy);
                    reducedAllBoundaryVertexMap.push_back(countRed);
                }
                countRed++;
            }
        }
    }


    index_t localToGlobalVertexMapper(index_t patchInd, index_t i){
        return globalVertexMap[patchInd * 4 + i - 1];
    }


    void reducedBasisEdgeMapper(const gsMultiPatch<> & geo)
    {
        index_t count = 1;

        for(size_t i = 0; i < numG1EdgeBasis.size(); i++) // Loop around global edges ( index - 1 )
        {

            for(const boundaryInterface &  item : geo.interfaces() )
            {
                if( ( item.second().patch * 4 + item.second().side() - 1) ==  i ) // If the global edge is an interface
                {
                    if( i == 0) // If the first edge of the 0-patch is an interface
                    {
                        for (size_t pi = 0; pi < numG1EdgeBasis[i]; pi++)
                        {
                            reducedG1BasisEdgeMap[pi] = count;
                            reducedG1BasisEdgeMap[numG1EdgeBasis[item.first().patch * 4 + item.first().side() - 1] + pi] = count;
                            count++;
                        }
                    }
                    else
                    {
                        for (size_t pi = 0; pi < numG1EdgeBasis[i] - numG1EdgeBasis[i-1]; pi++)
                        {
                            reducedG1BasisEdgeMap[numG1EdgeBasis[i-1] + pi] = count;
                            reducedG1BasisEdgeMap[numG1EdgeBasis[item.first().patch * 4 + item.first().side() - 2] + pi] = count;
                            count++;
                        }
                    }
                }
            }


            if( i == 0)
            {
                if(reducedG1BasisEdgeMap[numG1EdgeBasis[0]] == 0 )
                {
                    for (size_t pi = 0; pi < numG1EdgeBasis[i]; pi++)
                    {
                        reducedG1BasisEdgeMap[pi] = count;
                        count++;
                    }
                }
            }
            else
                if(reducedG1BasisEdgeMap[numG1EdgeBasis[i-1]] == 0 )
                {
                    for (size_t pi = 0; pi < numG1EdgeBasis[i] - numG1EdgeBasis[i-1]; pi++)
                    {
                        reducedG1BasisEdgeMap[numG1EdgeBasis[i-1] + pi] = count;
                        count++;
                    }
                }


        }

    }

    void printReducedBasisEdgeMapper()
    {
        for(size_t i = 0; i < reducedG1BasisEdgeMap.size(); i++)
            gsInfo << "Reduced Basis mapper: " << reducedG1BasisEdgeMap[i] << "\n";
    }


protected:

    std::vector<index_t> basisPerPatch;
    std::vector<std::pair<index_t, index_t>> dimensionPerBasis;

    std::vector<index_t> globalMap;
    std::vector<index_t> internalPatchMap;

    std::vector<index_t> globalEdgeMap;
    std::vector<index_t> reducedEdgeMap;
    std::vector<std::vector<index_t>> interfaceEdgeMap;
    std::vector<index_t> reducedBoundaryEdgeMap;

    std::vector<index_t> globalVertexMap;
    std::vector<index_t> reducedVertexMap;
    std::vector<std::vector<index_t>> internalVertexMap;
    std::vector<std::vector<index_t>> boundaryInterfaceVertexMap;
    std::vector<index_t> reducedBoundaryVertexMap;
    std::vector<index_t> reducedAllBoundaryVertexMap;


    std::vector<index_t> reducedG1BasisEdgeMap;
    std::vector<size_t> numG1EdgeBasis;
    std::vector<index_t> nPlusEdgeDimension;


};

}