/** @file gsTensorNurbs.h

    @brief Represents a tensor-product NURBS patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsCore/gsGeometry.h>

#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo
{

template<unsigned d, class T, class KnotVectorType> class gsTensorNurbs;

/** \brief 
    A tensor product Non-Uniform Rational B-spline function
    (NURBS) of parametric dimension \em d, with arbitrary target
    dimension.

    This is the geometry type associated with gsTensorNurbsBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector the NURBS bases use

    \ingroup geometry
    \ingroup Nurbs
*/

template<unsigned d, class T, class KnotVectorType>
class gsTensorNurbs : public gsGeoTraits<d,T>::GeometryBase
{

public: 
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef T Scalar_t;
    
    typedef gsTensorBSplineBasis<d,T,KnotVectorType> TBasis;      // underlying tensor basis
    
    /// Family type
    typedef gsBSplineBasis<T,KnotVectorType>  Family_t;
    
    // rational version of tensor basis (basis for this geometry)
    typedef gsTensorNurbsBasis<d,T,KnotVectorType>   Basis;
    
    /// Shared pointer for gsTensorNurbs
    typedef memory::shared_ptr< gsTensorNurbs<d,T,KnotVectorType> > Ptr;

public:

  /// Default empty constructor
  gsTensorNurbs() : Base() { }

  gsTensorNurbs( const Basis & basis, const gsMatrix<T> & coefs ) :
      Base( basis, coefs ) { }

  gsTensorNurbs( const Basis & basis, gsMovable< gsMatrix<T> > coefs ) :
      Base( basis, coefs ) { }

  /// Construct 2D tensor NURBS by knot vectors, degrees and coefficient matrix
  gsTensorNurbs( gsKnotVector<T> const& KV1, gsKnotVector<T> const & KV2,
                 gsMovable< gsMatrix<T> > tcoefs)
    {
      assert(d==2);

      gsBSplineBasis<T>    * Bu    = new gsBSplineBasis<T>(KV1);
      gsBSplineBasis<T>    * Bv    = new gsBSplineBasis<T>(KV2);

      TBasis   *tbasis = new TBasis(Bu,Bv) ;//d==2
      
      this->m_basis = new Basis(tbasis) ;
      this->m_coefs = tcoefs;
    };

  /// Construct 2D tensor NURBS by knot vectors, degrees and coefficient matrix
  gsTensorNurbs( gsKnotVector<T> const& KV1, gsKnotVector<T> const & KV2,
                 gsMovable< gsMatrix<T> > tcoefs, gsMovable< gsMatrix<T> > wgts)
    {
      assert(d==2);

      gsBSplineBasis<T>    * Bu    = new gsBSplineBasis<T>(KV1);
      gsBSplineBasis<T>    * Bv    = new gsBSplineBasis<T>(KV2);

      TBasis   *tbasis = new TBasis(Bu,Bv) ;//d==2
      
      this->m_basis = new Basis(tbasis , wgts) ;
      this->m_coefs = tcoefs;
    };

  /// Construct 3D tensor NURBS by knot vectors, degrees and coefficient matrix
  gsTensorNurbs( gsKnotVector<T> const & KV1, 
                 gsKnotVector<T> const & KV2, 
                 gsKnotVector<T> const & KV3,
                 gsMovable< gsMatrix<T> > tcoefs, 
                 gsMovable< gsMatrix<T> > wgts )
    {
      assert(d==3);
      
      gsBSplineBasis<T> * Bu= new gsBSplineBasis<T>(KV1);
      gsBSplineBasis<T> * Bv= new gsBSplineBasis<T>(KV2);
      gsBSplineBasis<T> * Bw= new gsBSplineBasis<T>(KV3);
      TBasis *tbasis = new TBasis(Bu,Bv,Bw) ;//d==3
      
      Basis *rbasis;

      this->m_basis = new Basis(tbasis, wgts) ;
      this->m_coefs = tcoefs;
    };

  gsTensorNurbs( gsKnotVector<T> const & KV1, 
                 gsKnotVector<T> const & KV2, 
                 gsKnotVector<T> const & KV3,
                 gsMovable< gsMatrix<T> > tcoefs)
    {
      assert(d==3);
      
      gsBSplineBasis<T> * Bu= new gsBSplineBasis<T>(KV1);
      gsBSplineBasis<T> * Bv= new gsBSplineBasis<T>(KV2);
      gsBSplineBasis<T> * Bw= new gsBSplineBasis<T>(KV3);
      TBasis *tbasis = new TBasis(Bu,Bv,Bw) ;//d==3
      
      this->m_basis = new Basis(tbasis) ;
      this->m_coefs = tcoefs;
    };

  /// Construct 3D tensor B-Spline by knot vectors, degrees and coefficient matrix
  //gsTensorNurbs( gsTensorBasis<T,d> * const basis, gsMatrix<T> * const coefs ) :
  //  gsGeometry<T,d>( basis, coefs ) { };

  /// Construct nD tensor B-Spline 
  //gsTensorNurbs( gsTensorBasis<T,d> * const basis, gsMatrix<T> * const coefs ) :
  //  gsGeometry<T,d>( basis, coefs ) { };

    GISMO_BASIS_ACCESSORS

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Clone function. Used to make a copy of the (derived) geometry
  virtual gsTensorNurbs * clone() const
    { return new gsTensorNurbs(*this); }

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
  { os << "Tensor-NURBS geometry "<< "R^"<< this->parDim() << 
      " --> R^"<< this->geoDim()<< ", #control pnts= "<< this->coefsSize() <<": "
       << this->coef(0) <<" ... "<< this->coef(this->coefsSize()-1); 
    os << "\nweights: "
       << this->basis().weights().at(0) <<" ... "
       << this->basis().weights().at(this->coefsSize()-1)
       <<"\n" ;
    return os; }

//////////////////////////////////////////////////
// Additional members for tensor NURBS
//////////////////////////////////////////////////

  /// Access to i-th weight
  T & weight(int i) const { return this->basis().weight(i); }

  /// Returns the weights of the rational basis
  gsMatrix<T> & weights() const { return this->basis().weights(); }

  /// Returns the degree of the basis wrt direction i 
  unsigned degree(unsigned i) const 
    { return this->basis().source().component(i).degree(); }

  using Base::basis;

  /// Returns a basis component of the tensor BSpline geometry
  const gsBSplineBasis<T>& basis(unsigned i) const 
    {
      if ( i >= d )
        std::cout<<"gsTensorNurbs error: Asked for basis component bigger than "<< d <<".\n";
      return this->basis().source().component(i);
    }

/// Toggle orientation wrt coordinate k
void reverse(unsigned k)
{ 
    GISMO_ASSERT(d==2, "only 2D for now");

    gsVector<int> str(d); 
    gsVector<int> sz (d); 

    gsTensorBSplineBasis<d,T> & tbsbasis = this->basis().source();

    sz[0]  = tbsbasis.component(k ).size();
    sz[1]  = tbsbasis.component(!k).size();
    str[0] = tbsbasis.source().stride( k );
    str[1] = tbsbasis.source().stride(!k );
    
    gsMatrix<T> & w  = tbsbasis.weights();

    for  ( int i=0; i< sz[0]; i++ )
        for  ( int j=0; j< sz[1]/2; j++ )
        {
            this->m_coefs.row(i*str[0] + j*str[1] ).swap(
                this->m_coefs.row(i*str[0] + (sz[1]-j-1)*str[1] )
                );

            w.row(i*str[0] + j*str[1] ).swap(
                w.row(i*str[0] + (sz[1]-j-1)*str[1] )
                );
        }
    tbsbasis.component(k).reverse();
}


protected:
 // TODO Check function
 // check function: check the coefficient number, degree, knot vector ...


// Data members
private:

}; // class gsTensorNurbs


//////////////////////////////////////////////////
//////////////////////////////////////////////////


} // namespace gismo
