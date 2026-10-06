/** @file gsNurbsBasis.h

    @brief Represents a NURBS basis with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsRationalBasis.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBoehm.h>


namespace gismo
{

/** \brief
    A univariate NURBS basis.

    This is the rational version of gsBSplineBasis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
    \ingroup Nurbs
*/

template<class T>
class gsTensorNurbsBasis<1,T> : public gsRationalBasis<gsBSplineBasis<T> >
{
public:

    typedef gsKnotVector<T> KnotVectorType;

    typedef memory::unique_ptr<gsGeometry<T> > gsGeoPtr;
    typedef memory::unique_ptr<gsBasis<T> >    gsBasisPtr;
    
    /// The family type
    typedef gsBSplineBasis<T> Family_t;

    /// Associated geometry type
    typedef gsNurbs<T> GeometryType;

    typedef gsRationalBasis<gsBSplineBasis<T> > Base;

    /// Associated Boundary basis type
    typedef gsConstantBasis<T> BoundaryBasisType;

    typedef memory::shared_ptr< gsTensorNurbsBasis > Ptr;
    typedef memory::unique_ptr< gsTensorNurbsBasis > uPtr;

    /// Dimension of the parameter domain
    static const int Dim = 1;

public:

    /// @brief Construct a NURBS basis with unit weights.
    /// \param u0 starting parameter
    /// \param u1 end parameter parameter
    /// \param interior number of interior knots
    /// \param degree degree of the spline space
    /// \param mult_interior multiplicity at the interior knots
    gsTensorNurbsBasis(T u0, T u1, unsigned interior, int degree, unsigned mult_interior=1):
    Base( gsBSplineBasis<T>(u0,u1,interior,degree,mult_interior) )
    { 
        // if( ! check()  )
        //   gsWarn << "Warning: Inconsistent "<< *this<< "\n";
    }

    /// Default empty constructor
    gsTensorNurbsBasis() : Base() { }

    /// Construct NURBS basis by a Bspline basis pointer (consumed) plus weights
    gsTensorNurbsBasis( gsBSplineBasis<T> *bs, gsMatrix<T> w) :
    Base( bs, give(w) )  { }

    /// Construct NURBS basis by a Bspline basis plus weights
    gsTensorNurbsBasis( gsBSplineBasis<T> bs, gsMatrix<T> w) :
    Base(new gsBSplineBasis<T>(), give(w))
    {
        *this->m_src = give(bs);
    }

    /// Construct NURBS basis of a knot vector
    explicit gsTensorNurbsBasis( gsKnotVector<T> KV ) :
    Base( new gsBSplineBasis<T>(give(KV)) )
    { }

    /// Construct a rational counterpart of B-spline basis given by knots and weights
    gsTensorNurbsBasis(gsKnotVector<T> KV, gsMatrix<T> w) :
    Base( new gsBSplineBasis<T>(give(KV)), give(w) ) { }

    /// Copy Constructor 
    gsTensorNurbsBasis( const gsTensorNurbsBasis & o) : Base(o) { }

    virtual ~gsTensorNurbsBasis();

public:

    /// Clone function. Used to make a copy of a derived basis
    GISMO_CLONE_FUNCTION(gsTensorNurbsBasis)

    gsGeoPtr makeGeometry( gsMatrix<T>coefs ) const override;

    static gsBasisPtr create(std::vector<KnotVectorType> cKV, gsMatrix<T> weights);
    using Base::create;

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const override
    {
        os << "NURBS Basis: deg=" << this->degree()
           << ", size=" << this->size() << ", knot vector:\n";
        os << this->knots() << ", weights: [ ";   
        os << this->weights().transpose() << " ]";   
        return os;
    }

    using Base::source;

// ***********************************************
// Additional members which forward to gsBSplineBasis
// ***********************************************

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return this->source().domainStart(); };

    /// Returns the starting value of the domain of the basis
    T domainEnd() const { return this->source().domainEnd(); };

    /// Returns the index of the first active (ie. non-zero) basis function at point u
    inline unsigned firstActive(const T & u) const { return this->source().firstActive(u); }

    /** \brief Number of active basis functions at an arbitrary parameter value.
     *
     * This assumes that this number doesn't change for different parameters.
     */
    inline index_t numActive() const { return this->source().numActive(); }

    /// Returns the index of the first active (ie. non-zero) basis
    /// function at all columns (points) of u
    inline index_t firstActive(const gsMatrix<T,1> & u) const
    { return this->source().firstActive(u.value()); };

    // /// Returns the knot vector of the basis
    const KnotVectorType & knots() const { return this->source().knots(); }
    KnotVectorType & knots()       { return this->source().knots(); }

    /// Insert a new knot \a val with multiplicity \a i.
    void insertKnot( T val, int i = 1)
    { 
        // TO DO: There is also Oslo Algorithm and others
        gsBoehm( this->knots(), this->weights(), val, i ); 
    }

    /// Insert the new knots given by the range \a [begin..end).
    template <class It>
    void insertKnots( It begin, It end )
    {
        gsBoehmRefine( this->knots(), this->weights(), this->degree(), begin, end);
    }
  
    /// Refine the basis uniformly by inserting \a numKnots new knots per knot span.
    void uniformRefine(int numKnots = 1, int mul=1, int dir = -1) override
    { 
        GISMO_UNUSED(dir);
        // TO DO ; replace this with global refinemnt by
        // Lane-Riesenfeld-like  Algorithm
        std::vector<T> newKnots;
        this->knots().getUniformRefinementKnots(numKnots, newKnots,mul);
        this->insertKnots( newKnots.begin(), newKnots.end() );
    }
  
    /// Apply k-refinement to the basis i times
    void uniform_k_refine(int const & i = 1) 
    {
        GISMO_NO_IMPLEMENTATION
        //m_p += i;
        //m_knots->degreeElevate(i);
        //m_knots->uniformRefine();
    };

}; // class gsTensorNurbsBasis<1,T>

} // namespace gismo

// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNurbsBasis.hpp)
#else
#ifdef gsNurbsBasis_EXPORT
#include GISMO_HPP_HEADER(gsNurbsBasis.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsTensorNurbsBasis<1,real_t>;
}
#endif
// *****************************************************************
