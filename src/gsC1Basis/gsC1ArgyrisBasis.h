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
                gsOptionList & optionList)
                : m_mp(mp), m_patchID(patchID)
    {
        info = optionList.getSwitch("info");
        twoPatch = optionList.getSwitch("twoPatch");
        simplified = optionList.getSwitch("simplified");

        // 9 Subspaces for the single patch
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

        // For C1 Basis at vertex
        valenceOfVertex.resize(4);
        for (size_t i = 0; i<valenceOfVertex.size(); i++)
            valenceOfVertex[i] = 6; // to get for boundary vertex == 6

        // For boundary
        numDofsVertex.resize(4);
        for (size_t i = 0; i<numDofsVertex.size(); i++)
            numDofsVertex[i] = 1; // to get for boundary vertex == 1
    }

    static uPtr make(   const gsC1ArgyrisBasis& other)
    { return uPtr( new gsC1ArgyrisBasis( other ) ); }

    gsC1ArgyrisBasis() { } //destructor


public:
    // basisG1Container:
    // - Interior space: [0] : inner,
    // - Edge spaces:    [1] : west, [2] : east, [3] : south, [4] : north,
    // - Vertex spaces:  [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
    void setInnerBasis(gsTensorBSplineBasis<d, T> & innerBasis) { basisG1Container[0] = innerBasis; }
    gsTensorBSplineBasis<d, T> & getInnerBasis() { return basisG1Container[0]; }

    // side index: 1 == west, 2 == east, 3 == south, 4 == north
    void setEdgeBasis(gsTensorBSplineBasis<d, T> & edgeBasis, index_t side) { basisG1Container[side] = edgeBasis; }
    gsTensorBSplineBasis<d, T> & getEdgeBasis(index_t side) { return basisG1Container[side]; }

    // corner index: 1 == sw, 2 == se, 3 == nw, 4 == ne
    void setVertexBasis(gsTensorBSplineBasis<d, T> & vertexBasis, index_t corner) { basisG1Container[4+corner] = vertexBasis; }
    gsTensorBSplineBasis<d, T> & getVertexBasis(index_t corner) { return basisG1Container[4+corner]; }
    // basisG1Container END

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

    void print_spaces()
    {
        // Some tests:
        for (index_t i = 0; i < 9; i++)
            if(basisG1Container[i].getMinCellLength() != basisG1Container[i].getMaxCellLength())
                gsInfo << "Different mesh-sizes is not implemented! \n";


        gsInfo << "-------------------------- Spaces for patch " << m_patchID << " --------------------------\n";
        gsInfo << "Interior space: S_1(" << basisG1Container[0].degree(0) << ", [";
        std::vector<index_t> kv_mult = basisG1Container[0].knots(0).multiplicities();
        for (size_t j = 1; j < kv_mult.size()-1; j++)
            gsInfo << kv_mult[j] << " ";
        gsInfo << "], " << basisG1Container[0].getMinCellLength() <<") ";
        gsInfo << "x S_2(" << basisG1Container[0].degree(1) << ", [";
        std::vector<index_t> kv_mult2 = basisG1Container[0].knots(1).multiplicities();
        for (size_t j = 1; j < kv_mult2.size()-1; j++)
            gsInfo << kv_mult2[j] << " ";
        gsInfo << "], " << basisG1Container[0].getMinCellLength() <<")\n";
        gsInfo << "\n------ Edge space:\n";
        for (index_t i = 1; i < 5; i++)
        {
            gsInfo << (kindOfEdge[i-1] ? "Interface-edge" : "Boundary-edge") << " space: S_1(" << basisG1Container[i].degree(0) << ", [";
            std::vector<index_t> kv_mult = basisG1Container[i].knots(0).multiplicities();
            for (size_t j = 1; j < kv_mult.size()-1; j++)
                gsInfo << kv_mult[j] << " ";
            gsInfo << "], " << basisG1Container[i].getMinCellLength() <<") ";
            gsInfo << "x S_2(" << basisG1Container[i].degree(1) << ", [";
            std::vector<index_t> kv_mult2 = basisG1Container[i].knots(1).multiplicities();
            for (size_t j = 1; j < kv_mult2.size()-1; j++)
                gsInfo << kv_mult2[j] << " ";
            gsInfo << "], " << basisG1Container[i].getMinCellLength() <<")\n";
        }
        gsInfo << "\n------ Vertex space:\n";
        for (index_t i = 5; i < 9; i++)
        {
            gsInfo << (kindOfVertex[i-5] == -1 ? "Boundary-vertex" : (kindOfVertex[i-5] == 0 ? "Internal-vertex" : "Interface-vertex"))
                << " space: S_1(" << basisG1Container[i].degree(0) << ", [";
            std::vector<index_t> kv_mult = basisG1Container[i].knots(0).multiplicities();
            for (size_t j = 1; j < kv_mult.size()-1; j++)
                gsInfo << kv_mult[j] << " ";
            gsInfo << "], " << basisG1Container[i].getMinCellLength() <<") ";
            gsInfo << "x S_2(" << basisG1Container[i].degree(1) << ", [";
            std::vector<index_t> kv_mult2 = basisG1Container[i].knots(1).multiplicities();
            for (size_t j = 1; j < kv_mult2.size()-1; j++)
                gsInfo << kv_mult2[j] << " ";
            gsInfo << "], " << basisG1Container[i].getMinCellLength() <<")\n";
        }
        gsInfo << "\n------ Plus/Minus space:\n";
        for (index_t i = 0; i < 4; i++)
        {
            gsInfo << "Plus space: S_1(" << basisPlusContainer[i].degree() << ", [";
            std::vector<index_t> kv_mult = basisPlusContainer[i].knots().multiplicities();
            for (size_t j = 1; j < kv_mult.size()-1; j++)
                gsInfo << kv_mult[j] << " ";
            gsInfo << "], " << basisPlusContainer[i].getMinCellLength() <<") ";
            gsInfo << "Minus space: S_1(" << basisMinusContainer[i].degree() << ", [";
            std::vector<index_t> kv_mult2 = basisMinusContainer[i].knots().multiplicities();
            for (size_t j = 1; j < kv_mult2.size()-1; j++)
                gsInfo << kv_mult2[j] << " ";
            gsInfo << "], " << basisMinusContainer[i].getMinCellLength() <<")\n";
        }
        gsInfo << "\n------ Gluing data/Geo space:\n";
        for (index_t i = 0; i < 4; i++)
        {
            gsInfo << "Gluing data space: S_1(" << basisGluingDataContainer[i].degree() << ", [";
            std::vector<index_t> kv_mult = basisGluingDataContainer[i].knots().multiplicities();
            for (size_t j = 1; j < kv_mult.size()-1; j++)
                gsInfo << kv_mult[j] << " ";
            gsInfo << "], " << basisGluingDataContainer[i].getMinCellLength() <<") ";
            gsInfo << "Geo space: S_1(" << basisGeoContainer[i].degree() << ", [";
            std::vector<index_t> kv_mult2 = basisGeoContainer[i].knots().multiplicities();
            for (size_t j = 1; j < kv_mult2.size()-1; j++)
                gsInfo << kv_mult2[j] << " ";
            gsInfo << "], " << basisGeoContainer[i].getMinCellLength() <<")\n";
        }
    }

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

    void setValenceOfVertex(index_t i, index_t corner) { valenceOfVertex[corner-1] = i+3; } // valence plus 3

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
        {
            index_t dim_u = basisG1Container[0].component(0).size();
            index_t dim_v = basisG1Container[0].component(1).size();
            rowContainer[0] = (dim_u - 4) * (dim_v - 4);
        }

        // Interface basis functions
        for (index_t i = 0; i<4; ++i)
        {
            index_t dim_u = basisG1Container[i+1].component(0).size();
            index_t dim_v = basisG1Container[i+1].component(1).size();

            if (m_mp.isBoundary(m_patchID,i+1)) // +1 of side index
            {
                if (twoPatch)
                    rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size() - 8;
                else if (simplified)
                    rowContainer[1+i] = math::max( ((i+1 > 2) ? dim_u*2 : dim_v*2 ) - 10, 0);
                else
                    rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
                kindOfEdge[i] = false; // bdy
            }
            else
            {
                if (twoPatch)
                    rowContainer[1+i] = basisPlusContainer[i].size()+basisMinusContainer[i].size();
                else
                    rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
                kindOfEdge[i] = true; // interface
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
                rowContainer[1+4+i] = valenceOfVertex[i];
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

    index_t getPatchID() const { return m_patchID; }

    index_t cols(index_t side = 0) const { return colContainer[side]; }

    index_t rows(index_t side = 0) const { return rowContainer[side]; }

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

    // TODO REALLY BAD IMPLEMENTATION: NEED A HUGE NEW IMPLEMENTATION
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
                    if (simplified)
                        num = basisG1Container[side_id].component(side_id < 3 ? 1 : 0).size() - bdy_shift; // Boundary
                    else
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

            num = num < 0 ? 0 : num; // if there are no bf at the interface

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
            if (offset == 0 && rows(corner_id ) != 0) {

                if (twoPatch) {
                    index_t ii = 0;
                    gsMatrix<index_t> indizes(3, 1);
                    index_t start = rowBegin(corner_id); // The first 3 basis functions

                    for (index_t i = start; i < start + 3; i++, ii++) // Single basis function
                        indizes(ii, 0) = i;

                    return indizes;
                } else if (kindOfVertex[corner_id - 4 - 1] != 0)
                {
                    index_t ii = 0;
                    gsMatrix<index_t> indizes(valenceOfVertex[corner_id - 4 - 1] - numDofsVertex[corner_id - 4 - 1], 1);
                    index_t corner_id = side.index(); // + 4 already included!
                    index_t start = rowBegin(corner_id); // The first 3 basis functions
                    for (index_t i = start + numDofsVertex[corner_id - 4 - 1];
                         i < start + valenceOfVertex[corner_id - 4 - 1]; i++, ii++) // Single basis function
                        indizes(ii, 0) = i;
                    return indizes;
                }
                else
                {
                    gsMatrix<index_t> null(1, 1);
                    null(0, 0) = -1;
                    return null;
                }
            }
            else if (offset == 1 && !twoPatch)
            {
                index_t ii = 0;
                gsMatrix<index_t> indizes(valenceOfVertex[corner_id - 4 - 1], 1);
                index_t start = rowBegin(corner_id); // The first 3 basis functions
                for (index_t i = start;
                     i < start + valenceOfVertex[corner_id - 4 - 1]; i++, ii++) // Single basis function
                    indizes(ii, 0) = i;
                return indizes;
            }
        }
        gsMatrix<index_t> null(1, 1);
        null(0, 0) = -1;

        return null;
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

    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        // Using the inner basis for iterating
        return basisG1Container[0].makeDomainIterator();
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

    bool info, twoPatch, simplified;

    std::vector<gsTensorBSplineBasis<d, T>> basisG1Container;

    std::vector<gsBSplineBasis<T>> basisPlusContainer;
    std::vector<gsBSplineBasis<T>> basisMinusContainer;
    std::vector<gsBSplineBasis<T>> basisGeoContainer;

    std::vector<gsBSplineBasis<T>> basisGluingDataContainer;

    std::vector<index_t> colContainer;
    std::vector<index_t> rowContainer;

    std::vector<bool> kindOfEdge;
    std::vector<index_t> kindOfVertex;
    std::vector<index_t> valenceOfVertex;

    std::vector<index_t> numDofsVertex;

}; // Class gsC1ArgyrisBasis

} // Namespace gismo
