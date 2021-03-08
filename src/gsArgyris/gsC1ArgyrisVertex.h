/** @file gsC1ArgyrisVertex.h

    @brief Creates the C1 Argyris Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsArgyris/gsC1ArgyrisBasis.h>
#include <gsArgyris/gsC1ArgyrisAuxiliaryPatch.h>



namespace gismo {
template<short_t d, class T>
class gsC1ArgyrisVertex
{

private:
    typedef gsC1ArgyrisBasis<d, T> Basis;
    typedef typename std::vector<Basis> ArgyrisBasisContainer;
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

    /// Shared pointer for gsC1Argyris
    typedef memory::shared_ptr<gsC1ArgyrisVertex> Ptr;

    /// Unique pointer for gsC1Argyris
    typedef memory::unique_ptr<gsC1ArgyrisVertex> uPtr;


public:
    /// Empty constructor
    ~gsC1ArgyrisVertex() { }


    gsC1ArgyrisVertex(gsMultiPatch<T> const & mp,
                ArgyrisBasisContainer & bases,
                const std::vector<size_t> & patchesAroundVertex,
                const std::vector<size_t> & vertexIndices,
                const index_t & numVer,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_patchesAroundVertex(patchesAroundVertex),
                m_vertexIndices(vertexIndices), m_optionList(optionList)
    {
        m_auxPatches.clear();
        basisVertexResult.clear();

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            index_t vertex_1 = m_vertexIndices[i];
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsC1ArgyrisAuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        reparametrizeVertexPatches();

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            if (m_optionList.getSwitch("twoPatch")) {
                gsMultiPatch<> test_mp(m_auxPatches[0].getPatch());

                gsMultiPatch<> g1Basis;
                gsTensorBSplineBasis<d, T> basis_edge = m_auxPatches[0].getArygrisBasisRotated().getVertexBasis(
                        m_auxPatches[0].side()); // 0 -> u, 1 -> v

                index_t dim_u = basis_edge.component(0).size();
                index_t dim_v = basis_edge.component(1).size();
                for (index_t j = 0; j < 2; j++) // v
                {
                    for (index_t i = 0; i < 2; i++) // u
                    {
                        gsMatrix<> coefs;
                        coefs.setZero(dim_u * dim_v, 1);

                        coefs(j * dim_u + i, 0) = 1;

                        g1Basis.addPatch(basis_edge.makeGeometry(coefs));
                    }
                }
                m_auxPatches[0].parametrizeBasisBack(g1Basis);
                basisVertexResult.push_back(g1Basis);
            }
        }

/*
        std::string fileName;
        std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
        gsParaviewCollection collection(basename);

        for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
        {
            for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
            {
                fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(m_patchesAroundVertex[np]), basisVertexResult[np].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");

            }
        }
        collection.save();
*/
    }

    void saveBasisVertex(gsSparseMatrix<T> & system)
    {
        for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
        {
            index_t patch_1 = m_patchesAroundVertex[np];
            index_t corner = m_vertexIndices[np];

            index_t shift_row = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
            {
                shift_row += m_bases[np].size_rows();
                shift_col += m_bases[np].size_cols();
            }

            index_t ii = 0;
            for (index_t i = m_bases[patch_1].rowBegin("vertex",corner); i < m_bases[patch_1].rowEnd("vertex",corner); ++i, ++ii)
            {
                index_t jj = 0;
                for (index_t j = m_bases[patch_1].colBegin("vertex",corner); j < m_bases[patch_1].colEnd("vertex",corner); ++j, ++jj)
                    if (basisVertexResult[0].patch(ii).coef(jj,0)*basisVertexResult[0].patch(ii).coef(jj,0)>1e-25)
                        system.insert(shift_row+i,shift_col+j) = basisVertexResult[0].patch(ii).coef(jj,0);
            }



        }
    }

    void reparametrizeVertexPatches()
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            switch (m_auxPatches[i].side()) // == vertex
            {
                case 1:
//                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " not rotated\n";
                    break;
                case 4:
                    m_auxPatches[i].rotateParamAntiClockTwice();
//                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " rotated twice anticlockwise\n";
                    break;
                case 2:
                    m_auxPatches[i].rotateParamAntiClock();
//                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " rotated anticlockwise\n";
                    break;
                case 3:
                    m_auxPatches[i].rotateParamClock();
//                    gsInfo << "Patch: " << auxGeom[i].getGlobalPatchIndex() << " rotated clockwise\n";
                    break;
            }
        }
    }

    void checkOrientation(size_t i)
    {
        if (m_auxPatches[i].getPatch().orientation() == -1)
        {
            m_auxPatches[i].swapAxis();
//            gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";

            //this->swapBdy(i); //Swap boundary edge bool-value

            // Swap vertices index after swapping axis
            if(m_auxPatches[i].side() == 2)
                m_auxPatches[i].setSide(3);
            else if(m_auxPatches[i].side() == 3)
                m_auxPatches[i].setSide(2);
        }
    }


protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    ArgyrisBasisContainer & m_bases;

    const std::vector<size_t> & m_patchesAroundVertex;
    const std::vector<size_t> & m_vertexIndices;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    ArgyrisAuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

}; // gsC1ArgyrisVertex

} // namespace gismo
