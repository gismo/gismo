/** @file gsCompositeGeom.h

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
class gsCompositeGeom : public gsGeoTraits<d,T>::GeometryBase
{

private:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsMappedSingleBasis<d,T> Basis;

    /// Shared pointer for gsCompositeGeom
    typedef memory::shared_ptr< gsCompositeGeom > Ptr;

    /// Unique pointer for gsCompositeGeom
    typedef memory::unique_ptr< gsCompositeGeom > uPtr;

public:
    /// Default empty constructor
    gsCompositeGeom() : Base(), m_compBasis(nullptr) { }

    gsCompositeGeom( gsMultiPatch<T> const & mp,std::string pathToMap );

    /// Construct Geom by basis and coefficient matrix
    gsCompositeGeom( const gsMappedBasis<d,T> & basis, const gsMatrix<T> & coefs );

    /// Construct Geom by (projecting) a multipatch together with a coefficient matrix
    gsCompositeGeom( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & bmap );

    gsCompositeGeom( const gsCompositeGeom& other );

    gsCompositeGeom<d,T> & operator=( const gsCompositeGeom<d,T> & other );

    ~gsCompositeGeom()
    { delete m_compBasis; } //destructor

    GISMO_BASIS_ACCESSORS

public:

//////////////////////////////////////////////////
// Virtual base members with a new implementation
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    // GISMO_CLONE_FUNCTION(gsCompositeGeom)

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

}; // class gsCompositeGeom

} // namespace gismo
