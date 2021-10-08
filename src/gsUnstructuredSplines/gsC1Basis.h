/** @file gsC1Basis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

/*
    TO DO (@Pascal)
        - Put everything from gsC1AgyrisBasis.h in here and in gsC1Basis.hpp
 */

#pragma once

#include<gsCore/gsBasis.h>
#include<gsIO/gsOptionList.h>

namespace gismo
{

template<short_t d,class T>
class gsC1Basis  : public gsBasis<T>
{

public:
    /// Shared pointer for gsC1Basis
    typedef memory::shared_ptr< gsC1Basis > Ptr;

    /// Unique pointer for gsC1Basis
    typedef memory::unique_ptr< gsC1Basis > uPtr;

    // gsC1Basis() : m_something(nullptr)
    // { }

    gsC1Basis(gsMultiPatch<T> & mp, index_t patchID );

    // gsC1Basis( const gsC1Basis& other );

    ~gsC1Basis() {};

    gsOptionList options() {return m_options;}
    void defaultOptions();
    void setOptions(gsOptionList opt) {m_options.update(opt, gsOptionList::addIfUnknown); };

    void print_spaces();

    void uniformRefine();
    void swapAxis();

    void init();

    gsMatrix<index_t> boundaryOffset(boxSide const & side , index_t offset = 0) const;

    // boundary(boxSide const & side)

// implementations of gsBasis
public:

    static uPtr make(   const gsC1Basis& other)
    { return uPtr( new gsC1Basis( other ) ); }

    //gsC1Basis() { } //destructor

public:

GISMO_CLONE_FUNCTION(gsC1Basis)

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

    index_t size_cols() const {
        index_t sz = 0;
        for (size_t i = 0; i < colContainer.size(); ++i)
            sz += colContainer[i];
        return sz;
    }

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

    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t> & result) const
    {
        GISMO_ASSERT(u.rows() == d, "Dimension of the points in active_into is wrong");
        //if (u.cols() > 1)
        //    gsInfo << "Active_into only for one point computed\n";

        result.resize(0,1);
        for (index_t u_i = 0; u_i < u.cols(); ++u_i) // For each points
        {
            //index_t u_i = 0; // Check if the points are in the same element TODO
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
            result.conservativeResize(result.rows()+result_single.rows(), u.cols() );
            result.block(result.rows()-result_single.rows(),u_i, result_single.rows(), 1) = result_single;
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



// implementations of gsC1Basis
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
    index_t getBasisPlusSize(index_t side) const { return basisPlusContainer[side-1].size(); }

    void setBasisMinus(gsBSplineBasis<> & basisMinus, index_t side) { basisMinusContainer[side-1] = basisMinus; }
    gsBSplineBasis<> & getBasisMinus(index_t side) { return basisMinusContainer[side-1]; }
    index_t getBasisMinusSize(index_t side) const { return basisMinusContainer[side-1].size(); }

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

    void setValenceOfVertex(index_t i, index_t corner) { valenceOfVertex[corner-1] = i+3; } // valence plus 3

    void setNumDofsVertex(index_t i, index_t corner) { numDofsVertex[corner-1] = i; }

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

    // Data members
protected:
    /// Class options
    gsOptionList m_options;

    /// The multipatch
    gsMultiPatch<T> & m_patches;

    /// The ID of the single basis
    index_t m_patchID;

    // Collection of the subspaces
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

};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1Basis.hpp)
#endif
