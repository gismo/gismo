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
public:
    gsGeometrySlice(const gsGeometry<T>* geo, index_t fixed_dir,T par)
    : m_geo(geo),m_fixed_dir(fixed_dir),m_par(par)
    {
        GISMO_ASSERT(geo->domainDim()>fixed_dir,"Geometry has not big enough dimension to fix the given fixed_dim.");
    }

    /// \brief Gives back the domain dimension of this slice
    /// Note that this is one less than the domain dimension of the
    /// underlying geometry.
    int domainDim() const
    {
        return m_geo->domainDim()-1;
    }

    /// \brief Gives back the target dimension of this slice
    /// Note that this is the same as the target dimension of the
    /// underlying geometry.
    int targetDim() const
    {
        return m_geo->targetDim();
    }

    /// Clone function. Makes a deep copy of the geometry object.
    gsGeometrySlice * clone() const
    {
        return new gsGeometrySlice(m_geo,m_fixed_dir,m_par);
    }

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
    gsMatrix<T> parameterRange() const
    {
        const gsMatrix<T> fullRange = m_geo->parameterRange();
        const index_t rows = fullRange.rows()-1;
        const index_t cols = fullRange.cols();
        gsMatrix<T> range(rows,cols);
        range.topRows(m_fixed_dir) = fullRange.topRows(m_fixed_dir);
        range.bottomRows(rows - m_fixed_dir) = fullRange.bottomRows(rows - m_fixed_dir);
        return range;
    }

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
    const gsGeometry<T> * m_geo; // pointer to the goemetry object
    const index_t m_fixed_dir; // fixed parameter direction
    const T m_par; // value for the fixed direction

}; // class gsGeometrySlice


}
