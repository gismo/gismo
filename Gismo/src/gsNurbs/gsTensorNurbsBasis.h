
#pragma once

#include <gsCore/gsForwardDeclarations.h>


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
  */
  
template<unsigned d, class T, class KnotVectorType >
class gsTensorNurbsBasis : public gsRationalBasis< gsTensorBSplineBasis<d,T,KnotVectorType> >  
{

public: 
  /// Base type
  typedef gsRationalBasis< gsTensorBSplineBasis<d,T,KnotVectorType> > Base;

  /// Source basis type
  typedef gsTensorBSplineBasis<d,T,KnotVectorType> Src_t;

  /// Coordinate basis type
  typedef typename Src_t::Basis_t Basis_t;

  /// Coefficient type
  typedef T Scalar_t;

  /// Associated geometry type
  typedef gsTensorNurbs<d, T> GeometryType;

  /// Associated Boundary basis type
  using typename Base::BoundaryBasisType;

  typedef memory::shared_ptr< gsTensorNurbsBasis > Ptr;
    
  //using typename Base::iterator;
  //using typename Base::const_iterator;

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

  gsTensorNurbsBasis( Src_t* basis, gsMovable< gsMatrix<T> > w ) : Base(basis,w) { }

  gsTensorNurbsBasis(const gsTensorNurbsBasis & o) : Base(o) { }


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

  /// Clone function. Used to make a copy of the object
  gsTensorNurbsBasis * clone() const
    { return new gsTensorNurbsBasis(*this); }
  
    GISMO_MAKE_GEOMETRY_NEW

};


} // namespace gismo
