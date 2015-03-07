/** @file gsTensorBSpline.h

    @brief Represents a tensor-product B-spline patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsNurbs/gsBoehm.h>
#include <gsUtils/gsMultiIndexIterators.h>

namespace gismo
{


/*
template <unsigned d, typename T>
struct gsTenBSplTraits
{

};

template <typename T>
struct gsTenBSplTraits<1,T>
{
    typedef gsBSpline<T> BoundaryGeo;
};

template <typename T>
struct gsTenBSplTraits<2,T>
{
    typedef gsBSpline<T> BoundaryGeo;
};

template <typename T>
struct gsTenBSplTraits<3,T>
{
    typedef gsTensorBSpline<2, T> BoundaryGeo;
};
*/

/** \brief
    A tensor product of \em d B-spline functions, with arbitrary target dimension.

    This is the geometry type associated with gsTensorBSplineBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector the B-spline bases use

    \ingroup geometry
    \ingroup Nurbs
*/

template<unsigned d, class T, class KnotVectorType>
class gsTensorBSpline : 
    public gsGenericGeometry< gsTensorBSplineBasis<d,T,KnotVectorType> >
{

public: 
  typedef T Scalar_t;
  typedef gsTensorBSplineBasis<d,T,KnotVectorType> Basis;
  typedef gsGenericGeometry< Basis > Base;

  /// Shared pointer for gsTensorBSpline
  typedef memory::shared_ptr< gsTensorBSpline<d,T> > Ptr;

public:

    /// Default empty constructor
    gsTensorBSpline() : Base() { }

    /// Construct B-Spline by basis functions and coefficient matrix
    gsTensorBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
        Base( basis, coefs ) 
    { }
    
    /// Construct B-Spline by basis functions and coefficient matrix
    gsTensorBSpline( const Basis & basis, gsMovable< gsMatrix<T> > coefs ) :
        Base( basis, coefs ) 
    { }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and coefficient matrix
    gsTensorBSpline( KnotVectorType const & KV1, KnotVectorType const & KV2,
                     gsMovable< gsMatrix<T> > tcoefs)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 2 knot-vectors.");
        
        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        Basis *tbasis = new Basis(Bu,Bv) ;//d==2
        
        GISMO_ASSERT(tbasis->size()== tcoefs.ref().rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );

        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( KnotVectorType const & KV1, KnotVectorType const & KV2,
                     const gsMatrix<T> & tcoefs)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 2 knot-vectors.");
        
        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        Basis *tbasis = new Basis(Bu,Bv) ;//d==2
        
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and 4 corner vertices,
    /// the four boundary curves are linear interpolations of the corners,
    /// the size of the matrix *corner* is 4 by 3
    gsTensorBSpline(gsMatrix<T> const    & corner, 
                    KnotVectorType const & KV1, 
                    KnotVectorType const & KV2);

    /// Construct 3D tensor B-Spline by knot vectors, degrees and coefficient matrix
    gsTensorBSpline( KnotVectorType const & KV1, 
                     KnotVectorType const & KV2, 
                     KnotVectorType const & KV3,
                     gsMovable< gsMatrix<T> > tcoefs )
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 3 knot-vectors.");
        
        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        gsBSplineBasis<T,KnotVectorType> * Bw= new gsBSplineBasis<T,KnotVectorType>(KV3);
        Basis *tbasis = new Basis(Bu,Bv,Bw) ;//d==3
    
        GISMO_ASSERT(tbasis->size()== tcoefs.ref().rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }

    /// Construct 3D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( KnotVectorType const & KV1, 
                     KnotVectorType const & KV2, 
                     KnotVectorType const & KV3,
                     const gsMatrix<T> & tcoefs )
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 3 knot-vectors.");

        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        gsBSplineBasis<T,KnotVectorType> * Bw= new gsBSplineBasis<T,KnotVectorType>(KV3);
        Basis *tbasis = new Basis(Bu,Bv,Bw) ;//d==3
        
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }

    /// Construct 4D tensor B-Spline by knot vectors, degrees and coefficient matrix
    gsTensorBSpline( KnotVectorType const & KV1, 
                     KnotVectorType const & KV2, 
                     KnotVectorType const & KV3,
                     KnotVectorType const & KV4,
                     gsMovable< gsMatrix<T> > tcoefs )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 4 knot-vectors.");
        
        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        gsBSplineBasis<T,KnotVectorType> * Bw= new gsBSplineBasis<T,KnotVectorType>(KV3);
        gsBSplineBasis<T,KnotVectorType> * Bz= new gsBSplineBasis<T,KnotVectorType>(KV4);        
        Basis *tbasis = new Basis(Bu,Bv,Bw,Bz) ;//d==4
    
        GISMO_ASSERT(tbasis->size()== tcoefs.ref().rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }

    /// Construct 4D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( KnotVectorType const & KV1, 
                     KnotVectorType const & KV2, 
                     KnotVectorType const & KV3,
                     KnotVectorType const & KV4,
                     const gsMatrix<T> & tcoefs )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 4 knot-vectors.");

        gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
        gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
        gsBSplineBasis<T,KnotVectorType> * Bw= new gsBSplineBasis<T,KnotVectorType>(KV3);
        gsBSplineBasis<T,KnotVectorType> * Bz= new gsBSplineBasis<T,KnotVectorType>(KV4);
        Basis *tbasis = new Basis(Bu,Bv,Bw,Bz) ;//d==3
        
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
        "Coefficient matrix for the tensor B-spline does not have the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }

    
    
    gsTensorBSpline( gsTensorBSpline const & o ) : Base(o) { }
    
    /// Clone function. Used to make a copy of the geometry
    gsTensorBSpline * clone() const
    { 
        return new gsTensorBSpline( *this );
    };

    // Destructor
    ~gsTensorBSpline() { }

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
  { 
    os << "Tensor BSpline geometry "<< "R^"<< this->parDim() << 
        " --> R^"<< this->geoDim()<< ", #control pnts= "<< this->coefsSize() <<
      ": "<< this->coef(0) <<" ... "<< this->coef(this->coefsSize()-1); 
    os<<"\nBasis:\n" << this->basis() ;
    return os; 
  }


//////////////////////////////////////////////////
// Additional members for tensor B-Splines
//////////////////////////////////////////////////

  /// Returns the degree of the basis wrt direction i
  unsigned degree(const unsigned & i) const 
    { return this->basis().component(i).degree(); }

/// Return the face of side s



// Overlaps/reimplements functionality of gsGeometry::boundary(),
// and adds orientation to the result.
/* LOOK AT THE gsGeometry::boundary
 *
 * this function is commented in case someone needs this version
 * with additional orientation
gsGeometry<T> * boundary( boxSide const& s ) const
{
    // this sould be a vector
    gsMatrix<unsigned> * ind = this->m_basis->boundary(s);

    gsMatrix<T> * fs = new gsMatrix<T>(ind->rows(),this->geoDim());

    // CORRECT ORIENTATION
    if ( (s == boundary::north) || (s == boundary::west) )
    {
        ind->reverseInPlace();
    }


    for ( index_t i = 0; i< ind->rows(); ++i)
    {
        fs->row(i) = this->m_coefs.row((*ind)(i, 0));
    }


    if (this->parDim() == 1)
    {
        return new typename gsTenBSplTraits<d, T>::BoundaryGeo();
    }
    else if (this->parDim() == 2 || this->parDim() == 3)
    {
        return new typename gsTenBSplTraits<d, T>::BoundaryGeo(
                    *this->m_basis->boundaryBasis(s), *fs);
    }

    GISMO_ERROR("Function face is only implemented for dimension 1, 2 and 3!");
    return new gsBSpline<T>();
}
*/

/// Toggle orientation wrt coordinate k
void reverse(unsigned k)
{ 
    GISMO_ASSERT(d==2, "only 2D for now");

    gsVector<int> str(d); 
    gsVector<int> sz (d); 
    gsTensorBSplineBasis<d,T> & tbsbasis = this->basis();

    sz[0]  = tbsbasis.component(k).size();
    sz[1]  = tbsbasis.component(!k).size();
    str[0] = tbsbasis.stride( k );
    str[1] = tbsbasis.stride(!k );
    
    for  ( int i=0; i< sz[0]; i++ )
        for  ( int j=0; j< sz[1]/2; j++ )
        {
            this->m_coefs.row(i*str[0] + j*str[1] ).swap(
                this->m_coefs.row(i*str[0] + (sz[1]-j-1)*str[1] )
                );
        }
    tbsbasis.component(k).reverse();
}

/// Sets the resulting BSpline to be periodic in direction \param dir.
inline void setPeriodic( int dir )
{
    this->basis().setPeriodic( dir );
    this->m_coefs = this->basis().perCoefs( this->m_coefs, dir );
}

private:
/// Helper function for the slice function
/// selects the row of coefficients from fullCoefs that are suitable
/// for the isoparametric slice in \a dir_fixed with \a par.
void constructCoefsForSlice(index_t dir_fixed,T par,const gsTensorBSpline<d,T>& geo,gsMatrix<T>& result) const;

public:

/// Constucts an isoparametric slice of this tensorBSpline by fixing
/// \a par in direction \a dir_fixed. The resulting tensorBSpline has
/// one less dimension and is given back in \a result.
void slice(index_t dir_fixed,T par,gsTensorBSpline<d-1,T>& result) const;

protected:
// TODO Check function
// check function: check the coefficient number, degree, knot vector ...


// Data members
private:

}; // class gsBSpline

//////////////////////////////////////////////////
//////////////////////////////////////////////////

template<unsigned d, class T, class KnotVectorType>
gsTensorBSpline<d,T,KnotVectorType>::gsTensorBSpline(gsMatrix<T> const & corner, KnotVectorType const& KV1, KnotVectorType const & KV2)
{
  assert(d==2);

  gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
  gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
  Basis *tbasis = new Basis(Bu,Bv) ;//d==2

  int n1 = KV1.size() - KV1.degree() - 1;
  int n2 = KV2.size() - KV2.degree() - 1;

  GISMO_ASSERT( (corner.rows()==4) && (corner.cols()==3),
           "gsTensorBSpline: Please make sure that the size of *corner* is 4-by-3");

  gsMatrix<T> pcp (n1*n2, 3);
  // set up CPs on boundary first. The inner CPs on each boundary curve are
  // uniformly linear dependent on the two corner CPs
  int j=0; // boundary v=0
  for (int i=0; i<=n1-1; i++)
  {
      for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
      {
          pcp(i+j*n1,xi)=corner(0,xi) + i/((T)(n1-1))*( corner(1,xi) - corner(0,xi) );
      }
  }
  j=n2-1; // boundary v=1
  for (int i=0; i<=n1-1; i++)
  {
      for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
      {
          pcp(i+j*n1,xi)=corner(3,xi) + i/((T)(n1-1))*( corner(2,xi) - corner(3,xi) );
      }
  }
  int i=0; // boundary u=0;
  for (j=0; j<=n2-1; j++)
  {
      for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
      {
          pcp(i+j*n1,xi)=corner(0,xi) + j/((T)(n2-1))*( corner(3,xi) - corner(0,xi) );
      }
  }
  i=n1-1; // boundary u=1;
  for (j=0; j<=n2-1; j++)
  {
      for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
      {
          pcp(i+j*n1,xi)=corner(1,xi) + j/((T)(n2-1))*( corner(2,xi) - corner(1,xi) );
      }
  }
  // uniformly linear dependent in horizontal direction
  for (j=1; j<=n2-2; j++)
  {
    for (i=1; i<=n1-2; i++)
    {
      for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
      {
          pcp(i+j*n1,xi)=pcp(0+j*n1,xi) + i/((T)(n1-1))*( pcp(n1-1+j*n1,xi)-pcp(0+j*n1,xi) );
      }
    }
  }

  this->m_basis = tbasis;
  this->m_coefs.swap( pcp );
}

// todo: move to hpp
template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::slice(index_t dir_fixed,T par,gsTensorBSpline<d-1,T>& result) const
{
    typedef gsTensorBSplineBasis<d-1,T,KnotVectorType> newBasis;
    GISMO_ASSERT(d-1>0,"cannot take iso slice of a curve");
    // construct the d-1 basis
    std::vector<gsBSplineBasis<T,KnotVectorType>* > bases;
    for(unsigned i=0;i<d;++i)
        if(i!=dir_fixed)
            bases.push_back(new gsBSplineBasis<T,KnotVectorType>(this->basis().knots(i)));
    newBasis *tbasis = new newBasis(bases) ;

    const int mult = this->basis().knots(dir_fixed).multiplicity(par);
    const unsigned degree = this->basis().degree(dir_fixed);

    gsMatrix<T> coefs;
    if(mult>=degree)
    {
        // no knot insertion needed, just extract the right coefficients
        constructCoefsForSlice(dir_fixed,par,*this,coefs);
    }
    else
    {
        // clone the basis and inserting upto degree knots at par
        gsTensorBSpline<d,T>* clone = this->clone();

        gsVector<unsigned> strides;
        gsVector<int> intStrides;
        this->basis().stride_cwise(intStrides);
        strides=intStrides.cast<unsigned>();
        gsTensorBoehm<T,KnotVectorType,gsMatrix<T> >(
                    clone->basis().knots(dir_fixed),clone->coefs(),par,dir_fixed,
                    strides,degree-mult,true);

        // extract right ceofficients
        constructCoefsForSlice(dir_fixed,par,*clone,coefs);
        delete clone;
    }

    // construct the object
    result = gsTensorBSpline<d-1,T>(*tbasis, give(coefs) );
    delete tbasis;
}

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::constructCoefsForSlice(index_t dir_fixed,T par,const gsTensorBSpline<d,T>& geo,gsMatrix<T>& result) const
{
    const gsTensorBSplineBasis<d,T,KnotVectorType>& base = geo.basis();
    const gsMatrix<T>& fullCoefs=geo.coefs();
    // pick the right coefficients and store them in coefs
    const unsigned degree = base.degree(dir_fixed);
    const KnotVectorType& knots = base.knots(dir_fixed);
    const int index = knots.findspan(par)-degree;
    gsVector<index_t,d> sizes,lowerCorner,upperCorner;

    gsVector<unsigned,d> ssizes;
    base.size_cwise( ssizes );
    sizes=ssizes.template cast<index_t>();
    lowerCorner.setZero();
    lowerCorner(dir_fixed)=index;
    upperCorner=sizes;
    upperCorner(dir_fixed)=index+1;

    // to do: gsMatrix<index_t> ind = gsTensorBasis::coefSlice(dim_fixed, index) ?

    gsTensorGridIterator<index_t> gridIter(sizes);
    gsTensorGridIterator<index_t> *iter = gridIter.makeSubGridIterator(lowerCorner,upperCorner);
    index_t size=1;
    for(index_t i = 0;i<d;++i)
        if(dir_fixed!=i)
            size*=sizes(i);
    result.resize(size,fullCoefs.cols());
    index_t i=0;
    for(iter->first();iter->good();iter->next())
    {
        result.row(i)=fullCoefs.row(iter->flatIndex());
        ++i;
    }
}

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSpline.hpp)
#endif
