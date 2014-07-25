
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsNurbs/gsBezier.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

// forward declaration
template <class T> class gsBezier;

/// Traits for BernsteinBasis in more dimensions
template<unsigned d, class T>
struct gsTraits<gsBernsteinBasis<T>,d>
{
  typedef gsTensorBernsteinBasis<d,T>        TensorBasisType;
  typedef gsTensorBezier<d, T>               TensorGeometryType;
  typedef typename 
  gsTraits<gsBernsteinBasis<T>,d-1>::TensorBasisType
    TensorBoundaryType;

  // typedef gsRatTensorBernsteinBasis<d,T>     RationalBasisType;
  // typedef gsRatBezier<d, T>                  RationalGeometryType;
  // typedef typename 
  // gsTraits<gsBernsteinBasis<T>,d-1>::RationalBasisType
  //   RationalBoundaryType;
};

/// Traits for BernsteinBasis in 1 dimension: specialization for d=1
template<class T>
struct gsTraits<gsBernsteinBasis<T>,1>
{
  typedef gsBernsteinBasis<T>               TensorBasisType;
  typedef gsBezier<T>                       TensorGeometryType;
  typedef gsBernsteinBasis<T>               TensorBoundaryType;

  // typedef gsRatBernsteinBasis<T>            RationalBasisType;
  // typedef gsRatBezier<T>                    RationalBoundaryType;
  // typedef gsRatBezier<T>                    RationalGeometryType;
};

  /** @brief
      A univariate Bernstein basis.

      \tparam T coefficient type

      \ingroup basis
  */
  
template<class T>
class gsBernsteinBasis : public gsBasis<T>
{
public:
  typedef gsBasis<T> Base;

  /// Coefficient type
  typedef T Scalar_t;

  /// Associated geometry type
  typedef gsBezier<T> GeometryType;

  /// Associated Boundary basis type
  typedef gsBernsteinBasis<T> BoundaryBasisType;

  /// Dimension of the parameter domain
  static const int Dim = 1;

  /// Shared pointer for gsBernsteinBasis
  typedef memory::shared_ptr< gsBernsteinBasis > Ptr;

  static Ptr Make ( const gsKnotVector<T> & KV, const int & p)
    { return Ptr( new gsBernsteinBasis(KV,p) ); };

public:

  /// Default empty constructor
  gsBernsteinBasis()  : Base() { };

  /// Construct Bernstein basis along the knots KV and degree p
  gsBernsteinBasis ( const gsKnotVector<T> & KV, const int & p) :
      m_p(p), m_breaks(0, KV.unique())
    {
	if(p != KV.degree()  )
	    std::cout << "gsBernsteinBasis Warning: Knots deg="<< KV.degree()<< " different than "<< p <<"\n";
    };
    
    gsBernsteinBasis<T>(T const& u0, T const& u1, int const& p, unsigned const & interior= 0)
    { 
	m_breaks= gsKnotVector<T>(u0,u1,interior);
	m_p = p;
    };

    gsBernsteinBasis<T>( gsBSplineBasis<T> const & bsb)
    { 
	m_breaks= gsKnotVector<T>(0, bsb.knots().unique());
	m_p = bsb.degree();
    };

  ~gsBernsteinBasis() { }; //destructor

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Clone function. Used to make a copy of a derived basis
  gsBernsteinBasis * clone() const;

    GISMO_MAKE_GEOMETRY_NEW

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    {
      typename gsKnotVector<T>::const_iterator itr;
      os << "Bernstein Basis: deg=" << this->degree()
	 << ", size="<< this->size() 
         << ", domain= [";
      for ( itr= m_breaks.begin(); itr != m_breaks.end(); ++itr )
	  os << *itr << " ";
      os << "].\n";
      return os; };
    
  int dim() const { return Dim; }

  /// Returns the number of basis functions in the basis
  int size() const { return m_p * (m_breaks.size()-1) +1; }

  /// Returns the number of elements.
  int numElements() const { return m_breaks.size() - 1; }

  // Look at gsBasis class for a description
  gsBernsteinBasis<T>& component(unsigned i) const;

  // Look at gsBasis class for a description
  void anchors_into(gsMatrix<T> & result) const;

  // Look at gsBasis class for a description
  void connectivity(const gsMatrix<T> & nodes, 
                    gsMesh<T> & mesh) const;

  // Look at gsBasis class for a description
  void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const;

  /// Returns the indices of the basis functions that touch the domain boundary
  gsMatrix<unsigned> * boundary( ) const ;

  /// Returns the indices of the basis functions that touch the domain boundary
  gsMatrix<unsigned> * boundary(boundary::side const & s ) const;

  /// Returns the boundary basis for side s
  gsBernsteinBasis<T> * boundaryBasis(boundary::side const & s ) const 
  { 
    gsKnotVector<T> kv(0,1,0,1); 
    return new gsBernsteinBasis<T>(kv,0); 
  }

  /// Returns a bounding box for the basis' domain
  gsMatrix<T> support() const ;

  /// Returns a bounding box for the basis' domain
  gsMatrix<T> support(const unsigned & i) const ;

  /// Evaluates the non-zero basis functions at value u.  
  /// Adapted from Algorithm A2.2 from 'The NURBS BOOK' pg70.
  virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

  /// Evaluates i-th basis functions at value u.  
  virtual void evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

  /// Evaluates a Bernstein given by coefs at points u
  virtual void eval_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

  /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
  void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial)derivatives of the i-th basis function at (the columns of) u.
  void derivSingle_into(unsigned & i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial) derivatives of a Bernstein given by coefs at (the columns of) u.
  void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial) second derivatives of a Bernstein given by coefs at (the columns of) u.
  gsMatrix<T> * deriv2(const gsMatrix<T> & u ) const ;

  /// Evaluates the laplacians of non-zero basis functions at (the columns of) u.
  gsMatrix<T> * laplacian(const gsMatrix<T> & u ) const ;
   
  /// Check the BernsteinBasis for consistency
  bool check() const
    {   // TO DO
	return true;
    };
    

//////////////////////////////////////////////////
// Additional members for univariate Bernstein basis
//////////////////////////////////////////////////

  /// Evaluates the non-zero basis functions and their
  /// first k derivatives at value u.

  /// Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72.
  virtual void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

  // Look at gsBasis class for a description
  int degree(int i = 0) const 
  { 
      GISMO_ASSERT(i==0,"Asked for degree(i) in 1D basis.");
      return m_p; 
  }

  /// Sets the degree of the basis as a reference
  inline void setDegree(int const & i) { m_p=i; };

  /// Returns the order of the basis
  inline unsigned order() const { return m_p+1; };

  /// Returns the starting value of the domain of the basis
  T domainStart() const { return m_breaks.first(); };

  /// Returns the ending value of the domain of the basis
  T domainEnd() const { return m_breaks.last(); };

  /// Returns the index of the first active (ie. non-zero) basis function at point u
  /// Takes into account non-clamped knots.
  inline unsigned firstActive(const T & u) const { 
     return m_p * m_breaks.findspan(u); 
  };

  /** \brief Number of active basis functions at an arbitrary parameter value.
   *
   * This assumes that this number doesn't change for different parameters.
   */
  inline unsigned numActive() const { return m_p + 1; }

  /// Returns the index of the first active (ie. non-zero) basis
  /// function at all columns (points) of u
  inline gsMatrix<unsigned,1> * firstActive(const gsMatrix<T,1> & u) const { 
      gsMatrix<unsigned,1> * res = m_breaks.findspan(u);
      res->array() *= m_p;
    return res;
  };

  /// Returns the knot vector of the basis
  gsKnotVector<T> * domain() const { return const_cast<gsKnotVector<T>*>(&m_breaks); };

  /// Insert a knot
  void insertKnot(T const knot)
    { m_breaks.insert(knot); };
    
  void insertKnot_withCoefs(T const knot, gsMatrix<T> & coefs);


  /** Refine uniformly the basis by adding \a numKnots
   *  knots between every two distinct knots.
   */
  void uniformRefine(int numKnots = 1)
    { m_breaks.uniformRefine(numKnots); }

  void uniformRefine_withCoefs(gsMatrix<T> & coefs, int numKnots = 1);

  void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1);

  /// Apply k-refinement to the basis i times
  void uniform_k_refine(int const & i = 1) 
    { 
        m_p += i; 
        m_breaks.degreeElevate(i);
        m_breaks.uniformRefine();
    };
  
  void degreeElevate(int const & i = 1) { m_p+=i; };

  inline int trueSize(){ return size(); }

// Data members
private:

  // Degree
  int m_p;
  // Knot vector
  gsKnotVector<T> m_breaks;


}; // class gsBernsteinBasis


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBernsteinBasis.hpp)
#endif

 
