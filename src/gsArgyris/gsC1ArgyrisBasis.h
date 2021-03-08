/** @file gsC1ArgyrisBasis.h

    @brief Creates the C1 Argyris basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsCore/gsBasis.h>

namespace gismo
{
template<short_t d,class T>
class gsC1ArgyrisBasis : public gsBasis<T>
{

public:
    /// Shared pointer for gsC1Argyris
    typedef memory::shared_ptr< gsC1ArgyrisBasis > Ptr;

    /// Unique pointer for gsDPatch
    typedef memory::unique_ptr< gsC1ArgyrisBasis > uPtr;

public:

    gsC1ArgyrisBasis(gsMultiPatch<T> const & mp,
                const index_t & patchID,
                const gsOptionList & optionList)
                : m_mp(mp), m_patchID(patchID)
    {
        basisG1Container.resize(9);

        // For each side:
        basisMinusContainer.resize(4);
        basisPlusContainer.resize(4);
        basisGeoContainer.resize(4);
        basisGluingDataContainer.resize(4);

        // For size
        rowContainer.resize(9);
        colContainer.resize(9);
    }

    static uPtr make(   const gsC1ArgyrisBasis& other)
    { return uPtr( new gsC1ArgyrisBasis( other ) ); }

    gsC1ArgyrisBasis() { } //destructor


public:
    // [0] : inner, [1] : west, [2] : east, [3] : south, [4] : north,
    // [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
    void setInnerBasis(gsTensorBSplineBasis<d, T> & innerBasis) { basisG1Container[0] = innerBasis; }
    gsTensorBSplineBasis<d, T> & getInnerBasis() { return basisG1Container[0]; }

    void setEdgeBasis(gsTensorBSplineBasis<d, T> & edgeBasis, index_t side) { basisG1Container[side] = edgeBasis; }
    gsTensorBSplineBasis<d, T> & getEdgeBasis(index_t side) { return basisG1Container[side]; }

    void setVertexBasis(gsTensorBSplineBasis<d, T> & vertexBasis, index_t corner) { basisG1Container[4+corner] = vertexBasis; }
    gsTensorBSplineBasis<d, T> & getVertexBasis(index_t corner) { return basisG1Container[4+corner]; }

    void setBasisPlus(gsBSplineBasis<> & basisPlus, index_t side) { basisPlusContainer[side-1] = basisPlus; }
    gsBSplineBasis<> & getBasisPlus(index_t side) { return basisPlusContainer[side-1]; }

    void setBasisMinus(gsBSplineBasis<> & basisMinus, index_t side) { basisMinusContainer[side-1] = basisMinus; }
    gsBSplineBasis<> & getBasisMinus(index_t side) { return basisMinusContainer[side-1]; }

    void setBasisGeo(gsBSplineBasis<> & basisGeo, index_t side) { basisGeoContainer[side-1] = basisGeo; }
    gsBSplineBasis<> & getBasisGeo(index_t side) { return basisGeoContainer[side-1]; }

    void setBasisGluingData(gsBSplineBasis<> & basisGD, index_t side) { basisGluingDataContainer[side-1] = basisGD; }
    gsBSplineBasis<> & getBasisGluingData(index_t side) { return basisGluingDataContainer[side-1]; }

    void uniformRefine()
    {
        for (size_t i=0; i< basisG1Container.size(); ++i)
            basisG1Container[i].uniformRefine();

        for (size_t i=0; i< basisPlusContainer.size(); ++i)
            basisPlusContainer[i].uniformRefine();
        for (size_t i=0; i< basisMinusContainer.size(); ++i)
            basisMinusContainer[i].uniformRefine();
        for (size_t i=0; i< basisGeoContainer.size(); ++i)
            basisGeoContainer[i].uniformRefine();
        for (size_t i=0; i< basisGluingDataContainer.size(); ++i)
            basisGluingDataContainer[i].uniformRefine();
    }

    void swapAxis()
    {
        std::vector<gsTensorBSplineBasis<d, T>> temp1 = basisG1Container;

        std::vector<gsBSplineBasis<T>> temp2 = basisPlusContainer;
        std::vector<gsBSplineBasis<T>> temp3 = basisMinusContainer;
        std::vector<gsBSplineBasis<T>> temp4 = basisGeoContainer;

        std::vector<gsBSplineBasis<T>> temp5 = basisGluingDataContainer;

        for (size_t i=0; i< basisG1Container.size(); ++i)
        {
            gsTensorBSplineBasis<d, T> newTensorBasis(basisG1Container[i].knots(1),basisG1Container[i].knots(0));
            basisG1Container[i].swap(newTensorBasis);
        }

        /*
         * W -> S; E -> N; S -> W; N -> E
         */
/*        for (size_t i=0; i< basisMinusContainer.size(); ++i)
        {
            basisG1Container[i+1].swap(temp1[(i+2)%4+1]);
            basisG1Container[i+5].swap(temp1[(i+2)%4+5]);

            basisPlusContainer[i].swap(temp2[(i+2)%4]);
            basisMinusContainer[i].swap(temp3[(i+2)%4]);
            basisGeoContainer[i].swap(temp4[(i+2)%4]);
            basisGluingDataContainer[i].swap(temp5[(i+2)%4]);
        }

        for (index_t i = 0; i < 4; i++)
            gsInfo << "i " << i+1 << ": " << (i+2)%4+5 << "\n";
*/
    }

    void init()
    {
        // Col == number of coefs
        // Row == number of basis functions

        // Cols:
        for (size_t i = 0; i<basisG1Container.size(); ++i)
            if (basisG1Container[i].size() == 1)
                colContainer[i] = 0;
            else
                colContainer[i] = basisG1Container[i].size();

        // Inner basis functions
        rowContainer[0] = basisG1Container[0].size(); // TODO change
        index_t dim_u = basisG1Container[0].component(0).size();
        index_t dim_v = basisG1Container[0].component(1).size();
        rowContainer[0] = (dim_u - 4)*(dim_v - 4);

        // Interface basis functions
        for (index_t i = 0; i<4; ++i)
        {
            if (m_mp.isBoundary(m_patchID,i+1)) // +1 of side index
                rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size() - 8;
            else
                rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size();
        }


        // Vertex basis functions
        for (index_t i = 0; i<4; ++i)
            if (basisG1Container[4+i+1].size() == 1)
                rowContainer[1+4+i] = 0;
            else
                rowContainer[1+4+i] = 4;

        gsInfo << "Patch: " << m_patchID << "\n";
        for (size_t i = 0; i<colContainer.size(); ++i)
            gsInfo << colContainer[i] << ", ";
        gsInfo << "\n";
        for (size_t i = 0; i<rowContainer.size(); ++i)
            gsInfo << rowContainer[i] << ", ";
        gsInfo << "\n";
    }

    index_t cols(std::string type, index_t side = 0) const
    {
        index_t col_index = -1;
        if (type == "inner")
            col_index = colContainer[0];
        else if (type == "edge")
            col_index = colContainer[side];
        else if (type == "vertex")
            col_index = colContainer[4+side];

        return col_index;
    }

    index_t rows(std::string type, index_t side = 0) const
    {
        index_t row_index = -1;
        if (type == "inner")
            row_index = rowContainer[0];
        else if (type == "edge")
            row_index = rowContainer[side];
        else if (type == "vertex")
            row_index = rowContainer[4+side];

        return row_index;
    }

    index_t colBegin(std::string type, index_t side = 0) const
    {
        index_t col_index = 0;
        if (type == "inner")
            col_index = 0; // Nothing happens
        else if (type == "edge")
        {
            col_index += colContainer[0];
            for (index_t i = 1; i < side; ++i)
                col_index += colContainer[i];
        }
        else if (type == "vertex")
        {
            col_index += colContainer[0];
            for (index_t i = 1; i < 5; ++i) // add all sides
                col_index += colContainer[i];

            for (index_t i = 1; i < side; ++i)
                col_index += colContainer[4+i];
        }

        return col_index;
    }

    index_t colEnd(std::string type, index_t side = 0) const
    {
        index_t col_index = 0;
        if (type == "inner")
            col_index = colContainer[0]; // Nothing happens
        else if (type == "edge")
        {
            col_index += colContainer[0];
            for (index_t i = 1; i < side+1; ++i)
                col_index += colContainer[i];
        }
        else if (type == "vertex")
        {
            col_index += colContainer[0];
            for (index_t i = 1; i < 5; ++i) // add all sides
                col_index += colContainer[i];

            for (index_t i = 1; i < side+1; ++i)
                col_index += colContainer[4+i];
        }

        return col_index;
    }

    index_t rowBegin(std::string type, index_t side = 0) const
    {
        index_t row_index = 0;
        if (type == "inner")
            row_index = 0; // Nothing happens
        else if (type == "edge")
        {
            row_index += rowContainer[0];
            for (index_t i = 1; i < side; ++i)
                row_index += rowContainer[i];
        }
        else if (type == "vertex")
        {
            row_index += rowContainer[0];
            for (index_t i = 1; i < 5; ++i) // add all sides
                row_index += rowContainer[i];

            for (index_t i = 1; i < side; ++i)
                row_index += rowContainer[4+i];
        }

        return row_index;
    }

    index_t rowEnd(std::string type, index_t side = 0) const
    {
        index_t row_index = 0;
        if (type == "inner")
            row_index = rowContainer[0]; // Nothing happens
        else if (type == "edge")
        {
            row_index += rowContainer[0];
            for (index_t i = 1; i < side+1; ++i)
                row_index += rowContainer[i];
        }
        else if (type == "vertex")
        {
            row_index += rowContainer[0];
            for (index_t i = 1; i < 5; ++i) // add all sides
                row_index += rowContainer[i];

            for (index_t i = 1; i < side+1; ++i)
                row_index += rowContainer[4+i];
        }

        return row_index;
    }

    gsMatrix<index_t> boundaryOffset(boxSide const & side , index_t offset = 0) const
    {
        const short_t side_id = side.index();
        index_t num = basisPlusContainer[side_id-1].size() - 4;
        index_t num_vert = 0;

        std::vector<index_t> corner_id;
        switch (side_id) {
           case 1:
               corner_id.push_back(1);
               corner_id.push_back(3);
               break;
           case 2:
               corner_id.push_back(2);
               corner_id.push_back(4);
               break;
           case 3:
               corner_id.push_back(1);
               corner_id.push_back(2);
               break;
           case 4:
               corner_id.push_back(3);
               corner_id.push_back(4);
               break;
           default:
               break;
        }
        // Special case two Patch TODO
        if (rows("vertex", corner_id[0]) != 0)
           num_vert += 3;
        if (rows("vertex", corner_id[1]) != 0)
           num_vert += 3;

        index_t ii = 0;
        gsMatrix<index_t> indizes(num+num_vert,1);
        //for(index_t of = 0;of<=offset;++of)
        {
            index_t start = rowBegin("edge", side_id); // The first num basis functions

            for (index_t i = start; i < start+num; i++, ii++) // Single basis function
                indizes(ii,0) = i;
        }

        // TODO Better
        //for(index_t of = 0;of<=offset;++of)
        {

            for (size_t j = 0; j < corner_id.size(); j++)
            {
                index_t start = rowBegin("vertex", corner_id[j]); // The first 3 basis functions
                if (rows("vertex", corner_id[j]) != 0)
                {
                    for (index_t i = start; i < start+3; i++, ii++) // Single basis function
                        indizes(ii,0) = i;
                }
            }

        }

        if (offset > 0)
            gsInfo << "Not implemented! \n";

        return indizes;
    }

public:

    GISMO_CLONE_FUNCTION(gsC1ArgyrisBasis)

    short_t domainDim() const
    {
      return d;
    }

    void connectivity(const gsMatrix<T> & nodes,
                          gsMesh<T>   & mesh) const
    {
      GISMO_UNUSED(nodes); GISMO_UNUSED(mesh);
      GISMO_NO_IMPLEMENTATION;
    }

    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
      GISMO_UNUSED(coefs);
      GISMO_NO_IMPLEMENTATION;
    }

    std::ostream &print(std::ostream &os) const
    {
      GISMO_UNUSED(os);
      GISMO_NO_IMPLEMENTATION;
    }

    short_t degree(short_t deg) const
    {
        return 0;
    }

    index_t size() const {
        index_t sz = 0;
        for (size_t i = 0; i < rowContainer.size(); ++i)
            sz += rowContainer[i];
        return sz;
    }

    index_t size_rows() const { return size(); }

    index_t size_cols() const {
        index_t sz = 0;
        for (size_t i = 0; i < colContainer.size(); ++i)
            sz += colContainer[i];
        return sz;
    }

protected:

    gsMultiPatch<T> m_mp;
    index_t m_patchID;

    std::vector<gsTensorBSplineBasis<d, T>> basisG1Container;

    std::vector<gsBSplineBasis<T>> basisPlusContainer;
    std::vector<gsBSplineBasis<T>> basisMinusContainer;
    std::vector<gsBSplineBasis<T>> basisGeoContainer;

    std::vector<gsBSplineBasis<T>> basisGluingDataContainer;

    std::vector<index_t> colContainer;
    std::vector<index_t> rowContainer;


}; // Class gsC1ArgyrisBasis

} // Namespace gismo
