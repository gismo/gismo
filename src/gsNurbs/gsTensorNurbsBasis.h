/** @file gsTensorNurbsBasis.h

    @brief Provides declaration of TensorNurbsBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsRationalBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

/** \brief 
    A tensor product Non-Uniform Rational B-spline (NURBS) basis.

    This is the rational version of gsTensorBSplineBasis.

    \tparam d dimension of the parameter domain
    \tparam T coefficient type
    \tparam KnotVectorType  the knot vector type the underlying NURBS bases use

    \ingroup basis
    \ingroup Nurbs
*/
template<unsigned d, class T>
class gsTensorNurbsBasis : public gsRationalBasis<typename gsBSplineTraits<d,T>::Basis>
{

public: 
    typedef gsKnotVector<T> KnotVectorType;

    /// Base type
    typedef gsRationalBasis<typename gsBSplineTraits<d,T>::Basis> Base;

    /// Family type
    typedef gsBSplineBasis<T>  Family_t;

    /// Source basis type
    typedef typename gsBSplineTraits<d,T>::Basis Src_t;

    /// Coordinate basis type
    typedef typename Src_t::Basis_t Basis_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef typename gsBSplineTraits<d,T>::RatGeometry GeometryType;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d-1,T>::RatBasis BoundaryBasisType;

    /// Shared pointer for gsTensorNurbsBasis
    typedef memory::shared_ptr< gsTensorNurbsBasis > Ptr;

    /// Unique pointer for gsTensorNurbsBasis
    typedef memory::unique_ptr< gsTensorNurbsBasis > uPtr;
    
    //typedef typename Base::iterator iterator;
    //typedef typename Base::const_iterator const_iterator;

public:

    /// Constructors for gsTensorNurbsBasis
    gsTensorNurbsBasis( const KnotVectorType& KV1, const KnotVectorType& KV2 )
    : Base( new gsBSplineBasis<T>(KV1, KV1.degree()), new gsBSplineBasis<T>(KV2, KV2.degree()) )
    { }

    gsTensorNurbsBasis( const KnotVectorType& KV1, const KnotVectorType& KV2, const KnotVectorType& KV3 )
    : Base( new gsBSplineBasis<T>(KV1, KV1.degree()),
            new gsBSplineBasis<T>(KV2, KV2.degree()),
            new gsBSplineBasis<T>(KV3, KV3.degree()) )
    { }

    // TO DO: more constructors
    //gsTensorNurbsBasis( gsBSplineBasis * x,  gsBSplineBasis* y, Basis_t* z ) : Base(x,y,z) { };
    //gsTensorNurbsBasis( std::vector<Basis_t* > const & bb ) : Base(bb) { };


    // Constructors forwarded from the base class
    gsTensorNurbsBasis() : Base() { };

    gsTensorNurbsBasis( Src_t* basis ) : Base(basis) { }

    gsTensorNurbsBasis( const Src_t & basis ) : Base(basis) { }

    gsTensorNurbsBasis( std::vector<gsBasis<T>* > const & bb, 
                        gsMatrix<T> w )
    : Base(Src_t(bb), give(w)) { }

    gsTensorNurbsBasis( Src_t* basis, gsMatrix<T> w ) : Base(basis, give(w)) { }

    gsTensorNurbsBasis(const gsTensorNurbsBasis & o) : Base(o) { }

    GISMO_CLONE_FUNCTION(gsTensorNurbsBasis)
  
    GISMO_MAKE_GEOMETRY_NEW

public:

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "TensorNurbsBasis: dim=" << this->dim()<< ", size="<< this->size() << ".";
        for ( unsigned i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< this->m_src->component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    BoundaryBasisType * boundaryBasis(boxSide const & s ) const   
    { 
        typename Src_t::BoundaryBasisType * bb = m_src->boundaryBasis(s);
        gsMatrix<unsigned> ind = m_src->boundary(s);
        
        gsMatrix<T> ww( ind.size(),1);
        for ( index_t i=0; i<ind.size(); ++i)
            ww(i,0) = m_weights( (ind)(i,0), 0);
        
        return new BoundaryBasisType(bb, give(ww));// note: constructor consumes the pointer
    }

protected:
    using Base::m_src;
    using Base::m_weights;

};


} // namespace gismo
