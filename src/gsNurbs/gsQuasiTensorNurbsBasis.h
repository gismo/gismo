
#pragma once

#include <gsTensor/gsTensorBasis.h>
#include <gsNurbs/gsNurbsBasis.h>


namespace gismo
{

// forward declaration
template<unsigned d, class T, class KnotVectorType>
class gsQuasiTensorNurbsBasis;


/*template <unsigned d, class T, class KnotVectorType>
struct TensorVersionOf<d, gsNurbsBasis<T,KnotVectorType> >
{
    // Tensor basis for gsBSplineBasis
    typedef gsQuasiTensorNurbsBasis<d,T,KnotVectorType> R;
};*/

  /** 
      Class for a tensor B-spline basis

      \param T coefficient type
      \param d dimension of the parameter domain
  */

  
template<unsigned d, class T, class KnotVectorType = gsKnotVector<T> >
class gsQuasiTensorNurbsBasis : public gsTensorBasis< d, gsNurbsBasis<T> >  
{

public: 
  /// Base type
  typedef gsTensorBasis< d, gsNurbsBasis<T> > Base;

  /// Coordinate basis type
  typedef gsNurbsBasis<T> Basis_t;

  /// Coefficient type
  typedef T Scalar_t;

  /// Associated geometry type
  typedef gsTensorBSpline<d, T> GeometryType;

  /// Associated Boundary basis type
  //typedef typename TensorBasisBoundaryType<d, T>::R BoundaryBasisType;
  using typename Base::BoundaryBasisType;

  using typename Base::iterator;
  using typename Base::const_iterator;

public:
  /// Constructors for gsQuasiTensorNurbsBasis
  gsQuasiTensorNurbsBasis( const KnotVectorType& KV1, const KnotVectorType& KV2 )
      : Base( new gsNurbsBasis<T>(KV1, KV1.degree()), new gsNurbsBasis<T>(KV2, KV2.degree()) )
    { }

  gsQuasiTensorNurbsBasis( const KnotVectorType& KV1, const KnotVectorType& KV2, const KnotVectorType& KV3 )
      : Base( new gsNurbsBasis<T>(KV1, KV1.degree()),
              new gsNurbsBasis<T>(KV2, KV2.degree()),
              new gsNurbsBasis<T>(KV3, KV3.degree()) )
    { }

    // TO DO: more constructors



  // Constructors forwarded from the base class
  gsQuasiTensorNurbsBasis() : Base() { };

  gsQuasiTensorNurbsBasis( Basis_t* x,  Basis_t*  y) : Base(x,y) { };

  gsQuasiTensorNurbsBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) : Base(x,y,z) { };

  gsQuasiTensorNurbsBasis( std::vector<Basis_t* > const & bb ) : Base(bb) { };



public:

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    {
      os << "QuasiTensorNurbsBasis<" << this->dim()<< ">, size "<< this->size() <<".";
      for ( unsigned i = 0; i!=d; ++i )
	  os << "\n  Direction "<< i <<": "<< this->m_src->component(i).knots() <<" ";
      return os;
    }

  /// Clone function. Used to make a copy of the object
  gsQuasiTensorNurbsBasis * clone() const
    { return new gsQuasiTensorNurbsBasis(*this); }

  virtual gsGeometry<T> * makeGeometry( const gsMatrix<T> * coefs ) const
    { return NULL; }
    //    { return new GeometryType(this, coefs); }

};


} // namespace gismo
