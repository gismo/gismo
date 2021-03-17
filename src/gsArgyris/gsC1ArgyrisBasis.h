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
        info = optionList.getSwitch("info");
        twoPatch = optionList.getSwitch("twoPatch");

        basisG1Container.resize(9);

        // For each side:
        basisMinusContainer.resize(4);
        basisPlusContainer.resize(4);
        basisGeoContainer.resize(4);
        basisGluingDataContainer.resize(4);

        // For size
        rowContainer.resize(9);
        colContainer.resize(9);

        // For topology
        kindOfEdge.resize(4);
        kindOfVertex.resize(4);

        // For boundary
        numDofsVertex.resize(4);
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

    std::vector<gsTensorBSplineBasis<d, T>> & getBasisG1Container() { return basisG1Container; }
    std::vector<index_t> & getRowContainer() { return rowContainer; }
    std::vector<index_t> & getColContainer() { return colContainer; }


    // Kind of edge
    // true == interface
    // false == boundary edge
    bool isInterface(index_t side) const { return kindOfEdge[side-1]; }

    // Kind of vertex
    // -1 Boundary vertex
    // 0 Internal vertex
    // 1 Interface boundary vertey
    void setKindOfVertex(index_t i, index_t corner) { kindOfVertex[corner-1] = i; }
    index_t getKindOfVertex(index_t corner) { return kindOfVertex[corner-1]; }

    void setNumDofsVertex(index_t i, index_t corner) { numDofsVertex[corner-1] = i; }

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
        for (size_t i=0; i< basisG1Container.size(); ++i)
        {
            gsTensorBSplineBasis<d, T> newTensorBasis(basisG1Container[i].knots(1),basisG1Container[i].knots(0));
            basisG1Container[i].swap(newTensorBasis);
        }
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
        rowContainer[0] = basisG1Container[0].size();
        index_t dim_u = basisG1Container[0].component(0).size();
        index_t dim_v = basisG1Container[0].component(1).size();
        rowContainer[0] = (dim_u - 4)*(dim_v - 4);

        // Interface basis functions
        for (index_t i = 0; i<4; ++i)
        {
            if (m_mp.isBoundary(m_patchID,i+1)) // +1 of side index
            {
                if (twoPatch)
                    rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size() - 8;
                else
                    rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
                kindOfEdge[i] = false;
            }
            else
            {
                if (twoPatch)
                    rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size();
                else
                    rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
                kindOfEdge[i] = true;
            }

        }

        // Vertex basis functions
        for (index_t i = 0; i<4; ++i)
        {
            if (twoPatch)
            {
                if (basisG1Container[4+i+1].size() == 1)
                    rowContainer[1+4+i] = 0;
                else
                    rowContainer[1+4+i] = 4;
            }
            else
                rowContainer[1+4+i] = 6;
        }



        if (info)
        {
            gsInfo << "Patch: " << m_patchID << "\n";
            for (size_t i = 0; i < colContainer.size(); ++i)
                gsInfo << colContainer[i] << ", ";
            gsInfo << "\n";
            for (size_t i = 0; i < rowContainer.size(); ++i)
                gsInfo << rowContainer[i] << ", ";
            gsInfo << "\n";
            gsInfo << "Kind of Vertex\n";
            for (size_t i = 0; i < kindOfVertex.size(); ++i)
                gsInfo << kindOfVertex[i] << ", ";
            gsInfo << "\n";

        }


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

    index_t rowBegin(index_t side = 0) const
    {
        index_t row_index = 0;
        for (index_t i = 0; i < side; ++i)
            row_index += rowContainer[i];
        return row_index;
    }

    index_t rowEnd(index_t side = 0) const
    {
        index_t row_index = 0;
        for (index_t i = 0; i < side+1; ++i)
            row_index += rowContainer[i];
        return row_index;
    }

    index_t colBegin(index_t side = 0) const
    {
        index_t col_index = 0;
        for (index_t i = 0; i < side; ++i)
            col_index += colContainer[i];
        return col_index;
    }

    index_t colEnd(index_t side = 0) const
    {
        index_t col_index = 0;
        for (index_t i = 0; i < side+1; ++i)
            col_index += colContainer[i];
        return col_index;
    }

    gsMatrix<index_t> boundaryOffset(boxSide const & side , index_t offset = 0) const
    {
        if (side.index() < 5) // Edge
        {
            short_t side_id = side.index();
            index_t num = 0;
            if (offset == 0)
            {
                index_t bdy_shift = twoPatch ? 4 : 6;
                if (!kindOfEdge[side_id - 1])
                    num = basisPlusContainer[side_id - 1].size() - bdy_shift; // Boundary
                else
                    num = basisPlusContainer[side_id - 1].size() - (twoPatch ? 0 : 6); // Interface
            }
            else if (offset == 1)
            {
                index_t bdy_shift = twoPatch ? 4 : 6;
                if (!kindOfEdge[side_id - 1])
                    num = basisMinusContainer[side_id - 1].size() - bdy_shift; // Boundary might not used and wrong
                else
                    num = basisMinusContainer[side_id - 1].size() - (twoPatch ? 0 : 4); // Interface
            } else
                gsInfo << "Offset > 1 is not implemented! \n";


            index_t ii = 0;
            gsMatrix<index_t> indizes(num , 1);
            //for(index_t of = 0;of<=offset;++of)
            {
                index_t start = rowBegin(side_id); // The first num basis functions

                if (offset == 1) {
                    index_t bdy_shift = twoPatch ? 4 : 6;
                    if (!kindOfEdge[side_id - 1])
                        start += basisPlusContainer[side_id - 1].size() - bdy_shift; // Boundary
                    else
                        start += basisPlusContainer[side_id - 1].size() - (twoPatch ? 0 : 6); // Interface
                }

                for (index_t i = start; i < start + num; i++, ii++) // Single basis function
                    indizes(ii, 0) = i;
            }
            return indizes;
        }
        else if (side.index() > 4)
        {
            index_t corner_id = side.index(); // + 4 already included!
            if (offset == 0 && rows("vertex", corner_id - 4) != 0) {

                if (twoPatch) {
                    index_t ii = 0;
                    gsMatrix<index_t> indizes(3, 1);
                    index_t start = rowBegin(corner_id); // The first 3 basis functions

                    for (index_t i = start; i < start + 3; i++, ii++) // Single basis function
                        indizes(ii, 0) = i;

                    return indizes;
                } else
                {
                    index_t ii = 0;
                    gsMatrix<index_t> indizes(6 - numDofsVertex[corner_id - 4 - 1], 1);
                    index_t corner_id = side.index(); // + 4 already included!
                    index_t start = rowBegin(corner_id); // The first 3 basis functions
                    for (index_t i = start + numDofsVertex[corner_id - 4 - 1];
                         i < start + 6; i++, ii++) // Single basis function
                        indizes(ii, 0) = i;
                    return indizes;
                }
            }
            else if (offset == 1 && !twoPatch)
            {
                index_t ii = 0;
                gsMatrix<index_t> indizes(6, 1);
                index_t start = rowBegin(corner_id); // The first 3 basis functions
                for (index_t i = start;
                     i < start + 6; i++, ii++) // Single basis function
                    indizes(ii, 0) = i;
                return indizes;
            }
            else {
                    gsMatrix<index_t> null(1, 1);
                    null(0, 0) = -1;
                    return null;
            }
        }

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

    // Returm max degree of all the spaces, otherwise i =
    short_t degree(short_t dir) const
    {
        short_t deg = 0;
        for (size_t i=0; i< basisG1Container.size(); ++i)
            if (basisG1Container[i].degree(dir) > deg)
                deg = basisG1Container[i].degree(dir);

        return deg;
    }

    index_t size_rows() const {
        index_t sz = 0;
        for (size_t i = 0; i < rowContainer.size(); ++i)
            sz += rowContainer[i];
        return sz;
    }

    index_t size_cols() const { return size(); }

    index_t size() const {
        index_t sz = 0;
        for (size_t i = 0; i < colContainer.size(); ++i)
            sz += colContainer[i];
        return sz;
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & side) const
    {
        // Using the inner basis for iterating
        return basisG1Container[0].makeDomainIterator(side);
    }

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t> & result) const
    {
        GISMO_ASSERT(u.rows() == d, "Dimension of the points in active_into is wrong");
        //if (u.cols() > 1)
        //    gsInfo << "Active_into only for one point computed\n";

        result.resize(0,1);
        //for (index_t u_i = 0; u_i < u.cols(); ++u_i) // For each points
        {
            index_t u_i = 0; // Check if the points are in the same element TODO
            index_t shift = 0;
            gsMatrix<index_t> result_single(0,1);
            for (size_t i=0; i< basisG1Container.size(); ++i)
            {
                if (rowContainer[i] != 0)
                {
                    gsMatrix<index_t> result_temp(0,1);
                    basisG1Container[i].active_into(u.col(u_i), result_temp);
                    result_temp.array() += shift;
                    result_single.conservativeResize(result_single.rows()+result_temp.rows(), 1 );
                    result_single.bottomRows(result_temp.rows()) = result_temp;

                    shift += basisG1Container[i].size();
                }
            }
            result.conservativeResize(result.rows()+result_single.rows(), 1 );
            result.bottomRows(result_single.rows()) = result_single;
        }
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        result.resize(0, u.cols());
        for (size_t i=0; i< basisG1Container.size(); ++i)
        {
            if (rowContainer[i] != 0)
            {
                gsMatrix<T> result_temp;
                basisG1Container[i].eval_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;
            }
        }
    }

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        result.resize(0, u.cols());
        for (size_t i=0; i< basisG1Container.size(); ++i)
        {
            if (rowContainer[i] != 0)
            {
                gsMatrix<T> result_temp;
                basisG1Container[i].deriv_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;
            }
        }
    }

    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        result.resize(0, u.cols());
        for (size_t i=0; i< basisG1Container.size(); ++i)
        {
            if (rowContainer[i] != 0)
            {
                gsMatrix<T> result_temp;
                basisG1Container[i].deriv2_into(u, result_temp);
                result.conservativeResize(result.rows()+result_temp.rows(), result.cols());
                result.bottomRows(result_temp.rows()) = result_temp;
            }
        }
    }

protected:


    // Input
    gsMultiPatch<T> m_mp;
    index_t m_patchID;

    bool info, twoPatch;

    std::vector<gsTensorBSplineBasis<d, T>> basisG1Container;

    std::vector<gsBSplineBasis<T>> basisPlusContainer;
    std::vector<gsBSplineBasis<T>> basisMinusContainer;
    std::vector<gsBSplineBasis<T>> basisGeoContainer;

    std::vector<gsBSplineBasis<T>> basisGluingDataContainer;

    std::vector<index_t> colContainer;
    std::vector<index_t> rowContainer;

    std::vector<bool> kindOfEdge;
    std::vector<index_t> kindOfVertex;

    std::vector<index_t> numDofsVertex;

}; // Class gsC1ArgyrisBasis

} // Namespace gismo
