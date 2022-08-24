/** @file gsCompositeSingleBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsMSplines/gsMappedBasis.h>

namespace gismo
{

template<short_t d,class T> class gsMappedBasis;

  /**
     @brief Class gsMappedSingleBasis represents an indivisual .....of a

      Note that it does not own the underlying basis: if you delete
      the object the basis is not deleted. The lifetime of the
      underlying basis should be at least the lifetime of gsMappedSingleBasis.

  */
template<short_t d,class T>
class gsMappedSingleBasis : public gsBasis<T>
{
private:
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;
    typedef gsBasis<T> Base;

public:
    /// Shared pointer for gsMappedSingleBasis
    typedef memory::shared_ptr< gsMappedSingleBasis > Ptr;

    /// Unique pointer for gsMappedSingleBasis
    typedef memory::unique_ptr< gsMappedSingleBasis > uPtr;

private:
  /// Default empty constructor
  gsMappedSingleBasis() { }

public:

  /// Dimension of the parameter domain
  static const int Dim = d;

  /// Construct a basis function by a pointer to a basis and an index i
  gsMappedSingleBasis(gsMappedBasis<d,T> * basis, unsigned const & i = 0)
      : m_basis(basis),m_index(i)
  {
      GISMO_ASSERT( i<unsigned(m_basis->nPatches()),"Invalid basis function index" );
  }

  gsMappedSingleBasis( const gsMappedSingleBasis& other ) : Base( other )
  {
      m_basis = other.m_basis;
      m_index = other.m_index;
  }

  ~gsMappedSingleBasis() { } //destructor

public:

  short_t domainDim() const
  {
      return d;
  }

  void connectivity(const gsMatrix<T> & nodes,
                          gsMesh<T>   & mesh) const
  {
      GISMO_UNUSED(nodes); GISMO_UNUSED(mesh);
      GISMO_NO_IMPLEMENTATION;
      //m_basis->connectivity(nodes,mesh);
  }

  // Look at gsBasis class for a description
  size_t numElements() const { return m_basis->getBase(m_index).numElements(); }

  /*
  void refine(gsMatrix<T> const & boxes)
  {  m_basis->refine(m_index,boxes);  }

  void refineElements(std::vector<unsigned> const & boxes)
  {  m_basis->refineElements(m_index,boxes);  }
  */

  /// Returns the indices of active (non zero) basis functions at points (columns of) u, as a list of indices, in result
  void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
  {
          m_basis->active_into(m_index,u,result);
  }

  /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
  void numActive_into(const gsMatrix<T> & u, gsVector<unsigned>& result) const
  {
      // Assuming all patches have the same degree
      m_basis->numActive_into(m_index,u,result);
  }

  /// Returns a bounding box for the basis' domain
  gsMatrix<T> support() const
  {
        return m_basis->getBase(m_index).support();
  }


  /// Returns a bounding box for the basis' domain
  gsMatrix<T> support(const index_t & i) const
  {
      return m_basis->getBase(m_index).support(i);
  }

  /// Evaluates the non-zero basis functions at value u.
  void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
  {
       // m_basis->evalGlobal_into(m_index,u,result);
       m_basis->eval_into(m_index,u,result);
  }

  /// Evaluates i-th basis functions at value u.
  void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
  {
      m_basis->evalSingle_into(m_index,i,u,result);
  }

  /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
  void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
  {
      m_basis->deriv_into(m_index,u,result);
  }

  /// Evaluates the (partial)derivatives of the i-th basis function at (the columns of) u.
  void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
  {
      GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(result);
      GISMO_NO_IMPLEMENTATION;
  }

  /// Evaluates the (partial) derivatives of the nonzero basis functions at points \a u into \a result.
  void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
  {
      m_basis->deriv2_into(m_index,u,result);
  }

  /// @brief Evaluate the (partial) derivatives of the \a i-th basis function
  /// at points \a u into \a result.
  void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
  {
      GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(result);
      GISMO_NO_IMPLEMENTATION;
  }

  /// @brief Evaluate the nonzero basis functions and their derivatives up
  /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, std::vector<gsMatrix<T> >& result) const
  {
      m_basis->evalAllDers_into(m_index,u,n,result);
  }

  /// @brief Evaluate the basis function \a i and its derivatives up
  /// to order \a n at points \a u into \a result.
  void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
  {
      GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(n); GISMO_UNUSED(result);
      GISMO_NO_IMPLEMENTATION;
  }

  /// @brief Evaluate the (partial) derivative(s) of order \a n the \a i-th basis function
  /// at points \a u into \a result.
  void evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
  {
      GISMO_UNUSED(i); GISMO_UNUSED(u); GISMO_UNUSED(n); GISMO_UNUSED(result);
      GISMO_NO_IMPLEMENTATION;
  }

  GISMO_CLONE_FUNCTION(gsMappedSingleBasis)

  memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
  {
      GISMO_UNUSED(coefs);
      GISMO_NO_IMPLEMENTATION;
  }

  std::ostream &print(std::ostream &os) const
  {
      GISMO_UNUSED(os);
      GISMO_NO_IMPLEMENTATION;
  }

  /// Prints the object as a string with extended details.
  std::string detail() const
  {
      // By default just uses print(..)
      std::ostringstream os;
      print(os);
      return os.str();
  }

  //////////////////////////////////////////////////
  // Virtual member that may be implemented or not by the derived class
  //////////////////////////////////////////////////

  index_t size() const
  {
      return m_basis->size(m_index);
  }

  /// Returns the polynomial degree.
  short_t maxDegree() const
  {
      return degree();
  }

  /// Returns the polynomial degree.
  short_t minDegree() const
  {
      return degree();
  }

  /// Returns the polynomial degree.
  short_t degree() const
  {
      return m_basis->maxDegree();                                   // must fix this (just took max_degree)
  }

  /// Returns the polynomial degree.
  short_t degree(short_t i) const
  {
      return m_basis->degree(m_index,i);
  }

  /// Check the LagrangeBasis for consistency
  bool check() const
  {
      GISMO_NO_IMPLEMENTATION;
  }

  /// The gsBasisFun points to the i-th basis function of m_basis
  /// after calling this setter.
  void setBasis( unsigned const & i ) const
  {
    GISMO_ASSERT( i<unsigned(m_basis->nPatches()),"Invalid basis index" );
    m_index = i;
  }


  /// Point to the first basis function of the basis
  void first() const
  {
    m_index = 0;
  }

  /// Points to the first basis function of the basis
  /// or return false if this is the last basis function
  bool next() const
  {
      return ( ++m_index < (index_t)m_basis->nPatches() );
  }

  /// Return the 1-d basis of the underlying tensor product basis for the \a i-th parameter component.
  const gsBasis<T>& component(short_t i) const
  {
      return m_basis->getBase(m_index).component(i);
  }

  gsBasis<T>& component(short_t i)
  {
      return m_basis->getBase(m_index).component(i);
  }

  typename gsBasis<T>::domainIter makeDomainIterator() const
  {
      return m_basis->getBase(m_index).makeDomainIterator();
  }

  typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
  {
      return m_basis->getBase(m_index).makeDomainIterator(s);
  }

// Data members
private:
  gsMappedBasis<d,T> * m_basis;
  mutable index_t m_index;

}; // class gsMappedSingleBasis


} // namespace gismo
