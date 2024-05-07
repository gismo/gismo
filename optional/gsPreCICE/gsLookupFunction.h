/** @file gsLookupFunction.h

    @brief Provides declaration of gsLookupFunction class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-...)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

namespace gismo
{

/**
 * @brief      Class defining a function that looks up registered data on points
 *
 * @tparam     T     Number format
 */
template <class T>
class gsLookupFunction : public gsFunction<T>
{

    struct Compare
    {
        bool operator()(const gsVector<T> & a, const gsVector<T> & b) const
        {
            return std::lexicographical_compare(  a.begin(), a.end(), b.begin(), b.end());
        }
    };

public:
    typedef gsGeometry<T> Base;

    /// Shared pointer for gsLookupFunction
    typedef memory::shared_ptr< gsLookupFunction > Ptr;

    /// Unique pointer for gsLookupFunction
    typedef memory::unique_ptr< gsLookupFunction > uPtr;

    // /// Default constructor
    // gsLookupFunction() { }

    /**
     * @brief      Constructs a new instance of the gsLookupFunction
     *
     * @param      interface   The precice::SolverInterface (see \a gsPreCICE)
     * @param[in]  meshName      The ID of the mesh on which the data is located
     * @param[in]  dataName      The ID of the data
     * @param[in]  patches     The geometry
     * @param[in]  parametric  Specifies whether the data is defined on the parametric domain or not
     */
    gsLookupFunction(   const gsMatrix<T> & points,
                        const gsMatrix<T> & data   )
    :
    m_points(points),
    m_data(data),
    {
        GISMO_ASSERT(m_points.cols()==m_data.cols(),"Points and data must have the same number of columns");
        for (index_t k = 0; m_points.cols(); k++)
            m_map.insert({m_points.col(k),k}); // m_map.at(vector) returns the column index of vector
    }

    /// Constructs a function pointer
    static uPtr make(   const gsMatrix<T> & points,
                        const gsMatrix<T> & data   )
    { return uPtr(new gsLookupFunction(points, data)); }

    GISMO_CLONE_FUNCTION(gsLookupFunction)

    /// Access a piece
    const gsLookupFunction<T> & piece(const index_t) const
    {
        return *this;
    }

    /// See \a gsFunction
    virtual short_t domainDim() const
    { return m_points.rows(); }

    /// Gives the targetDomain, currently only scalar functions (todo)
    virtual short_t targetDim() const
    { return m_data.rows(); }

    /// See \a gsFunction
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        index_t col;
        result.resize(this->targetDim(),u.cols());
        result.setZero();
        for (index_t k = 0; u.cols(); k++)
        {
            GISMO_ASSERT(m_map.find(u.col(k))!=m_map.end(),"Coordinate " + u.col(k).transpose() + " not registered in the table");
            col = m_map.at(u.col(k));
            result.col(k) = m_data.col(col);
        }
    }

    /// See \a gsFunction
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result) const
    {
        result.resize(1);
        // This would be nice to have with higher-order (IGA) coupling of precice
        gsMatrix<T> tmp;
        this->eval_into(u,tmp);
        result[0]= tmp;
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsLookupFunction\n";
        return os;
    }

protected:

    const gsMatrix<T> & m_points;
    const gsMatrix<T> & m_data;

    std::map<gsVector<T>,index_t,Compare> m_map;

};

}
