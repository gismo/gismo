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
/// \todo use flipTensor to generalize to adn dimension
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
    this->m_coefs = this->basis().perCoefs( this->m_coefs, dir );
    this->basis().setPeriodic( dir );
}

private:

/// Helper function for the slice function
/// selects the row of coefficients from fullCoefs that are suitable
/// for the isoparametric slice in \a dir_fixed with \a par.
void constructCoefsForSlice(index_t dir_fixed,T par,
                            const gsTensorBSpline<d,T>& geo,
                            gsMatrix<T>& result) const;

public:

/// Constucts an isoparametric slice of this tensorBSpline by fixing
/// \a par in direction \a dir_fixed. The resulting tensorBSpline has
/// one less dimension and is given back in \a result.
void slice(index_t dir_fixed,T par,gsTensorBSpline<d-1,T>& result) const;

protected:
// TODO Check function
// check function: check the coefficient number, degree, knot vector ...

}; // class gsBSpline

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSpline.hpp)
#endif
