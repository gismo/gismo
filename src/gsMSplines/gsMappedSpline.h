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
#include <gsMSplines/gsMappedSingleSpline.h>
#include <gsMSplines/gsMappedSingleBasis.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{

// forward declarations of the mapper classes
template<short_t d,class T> class gsMappedSingleSpline;

template<short_t d,class T>
class gsMappedSpline : public gsGeoTraits<d,T>::GeometryBase
{
    friend class gsMappedSingleSpline<d,T>;

private:
    typedef typename gsGeoTraits<d,T>::GeometryBase GeoType;
    typedef gsBasis<T> BasisType;

    /// Shared pointer for gsMappedSpline
    typedef memory::shared_ptr< gsMappedSpline > Ptr;

    /// Unique pointer for gsMappedSpline
    typedef memory::unique_ptr< gsMappedSpline > uPtr;

public:
    /// Default empty constructor
    gsMappedSpline() : m_mbases(nullptr) { }

    /// Construct Geom by multipatch and transformation matrix.
    /// The original coefficients are projected in the gsMappedBasis
    gsMappedSpline( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & m );

    /// Construct Geom by basis and coefficient matrix
    gsMappedSpline( const gsMappedBasis<d,T> & mbases, const gsMatrix<T> & coefs );

    gsMappedSpline( const gsMappedSpline& other );

    gsMappedSpline<d,T> & operator=( const gsMappedSpline& other );

    ~gsMappedSpline()
    { delete m_mbases; } //destructor

    void init(const gsMappedBasis<d,T> & mbasis)
    {
        GISMO_ASSERT(mbasis.domainDim()==d, "Error in dimensions");

        m_ss.clear();
        m_ss.reserve(mbasis.nPieces());
        for ( index_t k=0; k!=mbasis.nPieces(); k++ )
        {
            m_ss.push_back( gsMappedSingleSpline<d,T>(this,k) );
        }
    }

    void init(const gsMappedBasis<d,T> & mbasis, const gsMatrix<T> & coefs)
    {
        GISMO_ASSERT(mbasis.domainDim()==d, "Error in dimensions");

        m_global.clear();
        m_global = coefs;

        m_mbases=mbasis.clone().release();

        m_ss.clear();
        m_ss.reserve(mbasis.nPieces());
        for ( index_t k=0; k!=mbasis.nPieces(); k++ )
        {
            m_ss.push_back( gsMappedSingleSpline<d,T>(this,k) );
        }
    }

public:

    index_t nPieces() const {return m_mbases->nPieces();}

    /// getter for (const) m_bases[i]
    BasisType const & getBase(int i) const
    { return m_mbases->getBase(i); }

    /// getter for m_bases[i]
    BasisType & getBase(int i)
    { return m_mbases->getBase(i); }

     /// getter for m_bases
    const std::vector<BasisType*> getBases() const
    { return m_mbases->getBases(); }

    /// getter for m_mapper
    gsWeightMapper<T> const & getMapper() const
    { return m_mbases->getMapper(); }

    /// getter for m_mapper
    gsWeightMapper<T> * getMapPointer() const
    { return m_mbases->getMapPointer(); }

    /// getter for m_topol
    gsBoxTopology const & getTopol() const
    { return m_mbases->getTopol(); }

//////////////////////////////////////////////////
// Virtual base members with a new implementation
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    typedef gsMappedSingleBasis<d, T> Basis;

    GISMO_BASIS_ACCESSORS

    GISMO_CLONE_FUNCTION(gsMappedSpline)

    short_t domainDim() const
    { return m_mbases->domainDim(); }

    short_t targetDim() const
    { return m_global.cols(); }

    /// returns the amount of patches of the multi patch
    size_t nPatches() const
    { return m_mbases->nPatches(); }

    const gsMappedSingleSpline<d,T> & piece(const index_t k) const { return m_ss[k]; }

    index_t size() const {return nPieces();}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { return m_mbases->print(os); }

    const gsMappedBasis<d,T> & getMappedBasis() const
    { return *m_mbases; }

    const gsMatrix<T> & getMappedCoefs() const
    { return m_global; }

    gsMappedBasis<d,T> & getMappedBasis()
    { return *m_mbases; }

    gsMultiPatch<T> exportToPatches() const;

    // support (domain of definition)
    gsMatrix<T> support(index_t k) const
    { return m_mbases->getBase(k).support(); }

    gsGeometry<T> * exportPatch(int i,gsMatrix<T> const & localCoef) const;

public:
    //////////////////////////////////////////////////
    // functions for evaluating and derivatives
    //////////////////////////////////////////////////

    /// \brief Evaluates nonzero basis functions of \a patch at point \a u into \a result.
    ///
    /// Let...\n
    /// \a d denote the dimension of the parameter domain.\n
    /// \a k denote the number of active (i.e., non-zero) basis functions (see active_into()).
    /// \a n denote the number of evaluation points.\n
    ///
    /// The \a n <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>n</em>.
    /// Each column of \a u represents one evaluation point.\n
    /// \n
    /// The gsMatrix <b>\a result</b> contains the computed function values in the following form:\n
    /// Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    /// The column contains the values of all active functions "above" each other.\n
    ///
    /// For example, for scalar basis functions \a Bi : (x,y,z)-> R, a colum represents\n
    /// (B1, B2, ... , Bn)^T,\n
    /// where the order the basis functions \a Bi is as returned by active() and active_into().
    ///
    /// \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>n</em>.
    /// See above for details.
    /// \param[in,out] result gsMatrix of size <em>k</em> x <em>n</em>.
    /// See above for details.
    ///
    void eval_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void deriv_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void deriv2_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// @brief Evaluate the nonzero basis functions of \a patch and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const unsigned patch, const gsMatrix<T> & u,
                          const int n, std::vector<gsMatrix<T> >& result ) const;



// Data members
protected:
    /// Underlying gsMappedBasis
    gsMappedBasis<d,T> * m_mbases;
    /// Coefficients on the mapped basis
    gsMatrix<T> m_global;
    /// Underlying gsMappedSpline per patch
    std::vector<gsMappedSingleSpline<d,T> > m_ss;

}; // class gsMappedSpline

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMappedSpline.hpp)
#endif
