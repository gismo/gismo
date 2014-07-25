
#pragma once

#include <gsCore/gsTensorBasis.h>

namespace gismo
{

// forward declaration
template<unsigned d, class T>
class gsTensorBernsteinBasis;

  /** 
      Class for a tensor product Bernstein basis

      \param T coefficient type
      \param d dimension of the parameter domain

      \ingroup basis
  */

  
template<unsigned d, class T>
class gsTensorBernsteinBasis : public gsTensorBasis< d, gsBernsteinBasis<T> >  
{

public: 
  /// Base type
  typedef gsTensorBasis< d, gsBernsteinBasis<T> > Base;

  /// Coordinate basis type
  typedef gsBernsteinBasis<T> Basis_t;

  /// Coefficient type
  typedef T Scalar_t;

  /// Associated geometry type
  typedef gsTensorBezier<d,T> GeometryType;

  /// Associated Boundary basis type
  //typedef typename TensorBasisBoundaryType<d, T>::R BoundaryBasisType;
  using typename Base::BoundaryBasisType;

  using typename Base::iterator;
  using typename Base::const_iterator;

public:

  // TO DO: temporary, REMOVE
  gsTensorBernsteinBasis( const gsTensorBasis<d,gsBernsteinBasis<T> > & bb)
    : Base( bb )
  { }

  // Constructors forwarded from the base class
  gsTensorBernsteinBasis() : Base() { };

  gsTensorBernsteinBasis( Basis_t* x,  Basis_t*  y) : Base(x,y) { };

  gsTensorBernsteinBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) : Base(x,y,z) { };

  gsTensorBernsteinBasis( std::vector<Basis_t* > const & bb ) : Base(bb) { };



public:

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    {
      os << "TensorBernsteinBasis of dimension " << this->dim()<< ", size "<< this->size() <<".";
      return os;
    }

  /// Clone function. Used to make a copy of the object
  gsTensorBernsteinBasis * clone() const
    { return new gsTensorBernsteinBasis(*this); }

    
    GISMO_MAKE_GEOMETRY_NEW

};


} // namespace gismo
