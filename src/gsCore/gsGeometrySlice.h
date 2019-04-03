/** @file gsGeometrySlice.h

    @brief Provides declaration of the gsGeometrySlice class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): F. Buchegger
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

/** 
    \brief gsGeometrySlice is a class representing an iso parametric slice of a geometry object.
    It stores a pointer to the geometry object, which is only valid as long as this object is alive.
    Methods for printing to paraview are available in gsWriteToParaview.h
*/
template<class T>
class gsGeometrySlice : public gsFunction<T>
{
    typedef gsFunction<T> Base;
public:
    /// Shared pointer for gsGeometrySlice
    typedef memory::shared_ptr< gsGeometrySlice > Ptr;

    /// Unique pointer for gsGeometrySlice
    typedef memory::unique_ptr< gsGeometrySlice > uPtr;

    gsGeometrySlice(const gsFunction<T>* geo, index_t fixed_dir,T par)
    : m_geo(geo),m_fixed_dir(fixed_dir),m_par(par)
    {
        GISMO_ASSERT(fixed_dir>=0 && geo->domainDim()>static_cast<int>(fixed_dir),"Geometry has not big enough dimension to fix the given fixed_dim.");
        GISMO_ASSERT(geo->domainDim()!=1,"Cannot take a slice of a curve.");
    }

    /// Copyconstructor for new gsGeometrySlice(*this)
    gsGeometrySlice(const gsGeometrySlice* geoSlice)
        : m_geo(geoSlice->m_geo),m_fixed_dir(geoSlice->m_fixed_dir),m_par(geoSlice->m_par)
    {
        GISMO_ASSERT(geoSlice->m_fixed_dir>=0 && geoSlice->m_geo->domainDim()>static_cast<int>(geoSlice->m_fixed_dir),"Geometry has not big enough dimension to fix the given fixed_dim.");
        GISMO_ASSERT(geoSlice->m_geo->domainDim()!=1,"Cannot take a slice of a curve.");
    }

    /// \brief Gives back the domain dimension of this slice
    /// Note that this is one less than the domain dimension of the
    /// underlying geometry.
    short_t domainDim() const
    {
        return m_geo->domainDim()-1;
    }

    /// \brief Gives back the target dimension of this slice
    /// Note that this is the same as the target dimension of the
    /// underlying geometry.
    short_t targetDim() const
    {
        return m_geo->targetDim();
    }

    GISMO_CLONE_FUNCTION(gsGeometrySlice)

    /// \brief Gives back the values of this slice at points \a u in \a result
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> fullU;
        getFullParMatrix(u,fullU);
        m_geo->eval_into(fullU,result);
    }

    /// \brief Gives back the derivative of this slice at points \a u in \a result
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> fullU;
        getFullParMatrix(u,fullU);
        m_geo->deriv_into(fullU,result);
    }

    /// \brief Gives back the parameterRange of this slice in a Matrix
    gsMatrix<T> support() const
    {
        const gsMatrix<T> fullRange = m_geo->support();
        const index_t rows = fullRange.rows()-1;
        const index_t cols = fullRange.cols();
        gsMatrix<T> range(rows,cols);
        range.topRows(m_fixed_dir) = fullRange.topRows(m_fixed_dir);
        range.bottomRows(rows - m_fixed_dir) = fullRange.bottomRows(rows - m_fixed_dir);
        return range;
    }

    inline gsMatrix<T> parameterRange() const { return support();}

private:

    /// \brief This function takes a point matrix u and adds the row of the
    /// fixed direction filled with the value for it.
    void getFullParMatrix(const gsMatrix<T>& u, gsMatrix<T>& fullU) const
    {
        const index_t rows = u.rows()+1;
        const index_t cols = u.cols();
        fullU.resize(rows,cols);
        fullU.topRows(m_fixed_dir) = u.topRows(m_fixed_dir);
        fullU.row(m_fixed_dir).setConstant(m_par);
        fullU.bottomRows(rows - m_fixed_dir-1) = u.bottomRows(rows - m_fixed_dir-1);
    }

private:
    const gsFunction<T> * m_geo; // pointer to the function object
    const index_t m_fixed_dir; // fixed parameter direction
    const T m_par; // value for the fixed direction

}; // class gsGeometrySlice


}
