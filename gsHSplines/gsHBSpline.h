
#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsHSplines/gsHBSplineBasis.h>

namespace gismo
{

/** \brief
    A hierarchical B-Spline function, in \em d dimensions.

    This is the geometry type associated with gsHBSplineBasis.

    \tparam T is the coefficient type

    \ingroup geometry
*/    
    
template<unsigned d, class T>
class gsHBSpline : public gsGenericGeometry<gsHBSplineBasis<d,T> >
{   
public: 
  typedef gsHBSplineBasis<d,T> Basis;
  typedef gsGenericGeometry< gsHBSplineBasis<d,T> > Base;

  /// Shared pointer for gsHBSpline
  typedef memory::shared_ptr< gsHBSpline<2,T> > Ptr;
    
public:
    
  /// Default empty constructor
  gsHBSpline() { }
  
  /// Construct HB-Spline by basis functions and coefficient matrix
  gsHBSpline( const Basis * basis, const gsMatrix<T> * coefs ) :
    Base( basis, coefs ) { }

  /// Construct HB-Spline by basis functions and coefficient matrix
  gsHBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs ) { }
  
  /// Construct B-Spline from a Tensor B-Spline
  gsHBSpline( const gsTensorBSpline<2,T> & tbsp )
  { 
    this->m_basis = new Basis(tbsp);
    this->m_coefs = tbsp->coefs();
  }
  
  /// Copy constructor
  gsHBSpline( const gsHBSpline & other )
  { 
    this->m_basis = other.basis().clone();
    this->m_coefs = other.coefs();
  }
  
  /// Clone the gsHBspline
  virtual gsHBSpline * clone() const
  { return new gsHBSpline(*this); };
  
  ~gsHBSpline() { }; //destructor   
  
public:
  
    //////////////////////////////////////////////////
    // Virtual member functions required by the base class
    //////////////////////////////////////////////////
    
    /// Returns the degree wrt direction i
    unsigned degree(const unsigned & i) const 
    { return this->basisComponent(i).degree(); };
    
  
}; // class gsHBSpline
  
    
//////////////////////////////////////////////////
//////////////////////////////////////////////////
    

    
    
}; // namespace gismo
    
