
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsTensorBernsteinBasis.h>

namespace gismo {


template <class T> class gsBernsteinBasis;
template <unsigned d, class T> class gsTensorBezier;


/** \brief
    A tensor product of \em d piecewise Bezier functions, with arbitrary target dimension.

    This is the geometry type associated with gsTensorBernsteinBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    \ingroup geometry
*/

template<unsigned d, class T>
class gsTensorBezier : public gsGenericGeometry< gsTensorBernsteinBasis<d,T > >
{ 

public: 
    typedef T Scalar_t;
    typedef gsTensorBernsteinBasis<d,T> Basis;
    typedef gsGenericGeometry< Basis >  Base;

public: 

  // TO DO : temporary, REMOVE
  gsTensorBezier( const gsTensorBasis<d,gsBernsteinBasis<T> > * basis, const gsMatrix<T> * coefs ) 
  { 
    this->m_basis =  new Basis(*basis) ;
    this->m_coefs = *coefs ;    
  }


  gsTensorBezier( const Basis * basis, const gsMatrix<T> * coefs )
      : Base (basis,coefs) { }

  gsTensorBezier( const Basis & basis, const gsMatrix<T> & coefs )
      : Base (basis,coefs) { }

  /// Clone function. Used to make a copy of the (derived) geometry
  gsTensorBezier * clone() const { return new gsTensorBezier(*this); };

  /// Prints the object as a string.
  virtual std::ostream &print(std::ostream &os) const { os<<"gsTensorBezier\n"; return os; };

};

}; // namespace gismo
