/** @file gsGenericTensorBasis.h

    @brief Provides declaration of QuasiTensorNurbsBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsTensor/gsTensorBasis.h>


namespace gismo
{

/*template <short_t d, class T, class KnotVectorType>
  struct TensorVersionOf<d, gsNurbsBasis<T,KnotVectorType> >
  {
  // Tensor basis for gsBSplineBasis
  typedef gsQuasiTensorNurbsBasis<d,T,KnotVectorType> R;
  };*/

/** 
    @brief Class for a quasi-tensor B-spline basis

    \param T coefficient type
    \param d dimension of the parameter domain
    
    \ingroup Tensor
*/
template<short_t d, class T>
class gsGenericTensorBasis : public gsTensorBasis<d,T>
{

public: 
    /// Base type
    typedef gsTensorBasis<d,T> Base;

    /// Coefficient type
    typedef T Scalar_t;

    typedef gsBasis<T> Basis_t;

    /// Shared pointer for gsGenericTensorBasis
    typedef memory::shared_ptr< gsGenericTensorBasis > Ptr;

    /// Unique pointer for gsGenericTensorBasis
    typedef memory::unique_ptr< gsGenericTensorBasis > uPtr;

    // Associated geometry type
    typedef gsGenericGeometry<d,T> GeometryType;

    // Associated Boundary basis type
    //typedef typename Base::BoundaryBasisType BoundaryBasisType;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;

public:
    // Constructors forwarded from the base class
    gsGenericTensorBasis() : Base() { }

    gsGenericTensorBasis( Basis_t* x,  Basis_t*  y) : Base(x,y) { }

    gsGenericTensorBasis( const Basis_t & x, const Basis_t & y) :
    Base(x.clone().release(), y.clone().release()) { }

    gsGenericTensorBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) : Base(x,y,z) { }

    gsGenericTensorBasis( const Basis_t & x, const Basis_t & y, const Basis_t & z ) :
        Base(x.clone().release(), y.clone().release(), z.clone().release()) { }

    gsGenericTensorBasis( std::vector<Basis_t*> & bb ) : Base(bb.data()) 
    { 
        GISMO_ENSURE( d == bb.size(), "Wrong d in the constructor of gsTensorBSplineBasis." );
    }



public:

    GISMO_MAKE_GEOMETRY_NEW

    std::ostream &print(std::ostream &os) const
    {
        os << "GenericTensorNurbsBasis<" << this->dim()<< ">, size "<< this->size() <<".";
        for ( unsigned i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< this->component(i) <<" ";
        return os;
    }

    GISMO_CLONE_FUNCTION(gsGenericTensorBasis)

};


} // namespace gismo
