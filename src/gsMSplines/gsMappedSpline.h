/** @file gsMappedSpline.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once
#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsMSplines/gsMappedSingleBasis.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{
template<short_t d,class T>
class gsMappedSpline : public gsGeoTraits<d,T>::GeometryBase
{

private:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsMappedSingleBasis<d,T> Basis;

    /// Shared pointer for gsMappedSpline
    typedef memory::shared_ptr< gsMappedSpline > Ptr;

    /// Unique pointer for gsMappedSpline
    typedef memory::unique_ptr< gsMappedSpline > uPtr;

public:
    /// Default empty constructor
    gsMappedSpline() : Base(), m_compBasis(nullptr) { }

    gsMappedSpline( gsMultiPatch<T> const & mp,std::string pathToMap );

    /// Construct Geom by basis and coefficient matrix
    gsMappedSpline( const gsMappedBasis<d,T> & basis, const gsMatrix<T> & coefs );

    /// Construct Geom by (porjecting) a multipatch together with a coefficient matrix
    gsMappedSpline( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & bmap );

    gsMappedSpline( const gsMappedSpline& other );

    gsMappedSpline<d,T> & operator=( const gsMappedSpline& other );

    ~gsMappedSpline()
    { delete m_compBasis; } //destructor

    GISMO_BASIS_ACCESSORS

public:

//////////////////////////////////////////////////
// Virtual base members with a new implementation
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    GISMO_CLONE_FUNCTION(gsMappedSpline)

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        m_compBasis->print(os);
        return os;
    }

    const gsMappedBasis<d,T> & getCompBasis() const
    { return *m_compBasis; }

    gsMappedBasis<d,T> & getCompBasis()
    { return *m_compBasis; }

    /// The gsBasisFun points to the i-th basis function of m_basis
    /// after calling this setter.
    void setBasis( unsigned const & i ) const
    { this->basis().setBasis(i); }


    /// Point to the first basis function of the basis
    void first()
    { this->basis().first(); }

    /// Points to the first basis function of the basis
    /// or return false if this is the last basis function
    bool next() const
    {  return this->basis().next(); }

    gsMultiPatch<T> exportToPatches() const;

    gsGeometry<T> * exportPatch(int i,gsMatrix<T> const & localCoef) const
    { return m_compBasis->exportPatch(i,localCoef); }

// Data members
protected:
    gsMappedBasis<d,T> * m_compBasis;

}; // class gsMappedSpline

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMappedSpline.hpp)
#endif