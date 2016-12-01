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

namespace gismo
{

namespace internal
{
  using Eigen::internal::enable_if;
}

/** \brief
    A tensor product of \em d B-spline functions, with arbitrary target dimension.

    This is the geometry type associated with gsTensorBSplineBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    \ingroup geometry
    \ingroup Nurbs
*/
template<unsigned d, class T>
class gsTensorBSpline GISMO_FINAL : public gsGeoTraits<d,T>::GeometryBase
{
public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef T Scalar_t;

    typedef typename gsBSplineTraits<d,T>::Basis Basis;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    /// Family type
    typedef gsBSplineBasis<T>  Family_t;

    /// Shared pointer for gsTensorBSpline
    typedef memory::shared_ptr< gsTensorBSpline<d,T> > Ptr;

    /// Unique pointer for gsTensorBSpline
    typedef memory::unique_ptr< gsTensorBSpline> uPtr;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d-1,T>::Geometry BoundaryGeometryType;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d-1,T>::Basis BoundaryBasisType;

public:

    /// Default empty constructor
    gsTensorBSpline() { }

    // Construct B-Spline by basis functions and coefficient matrix
    //gsTensorBSpline( const gsConstantBasis<T> & basis, const gsMatrix<T> & coefs )
    //{ }

    /// Construct B-Spline by basis functions and coefficient matrix
    gsTensorBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs ) 
    { }

    /// Construct B-Spline by basis functions and coefficient matrix
    gsTensorBSpline( const Basis & basis, gsMovable< gsMatrix<T> > coefs ) :
    Base( basis, coefs ) 
    { }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and coefficient matrix
    template<typename U>
    gsTensorBSpline( gsKnotVector<T> const & KV1, gsKnotVector<U> const & KV2,
                     gsMovable< gsMatrix<T> > tcoefs,
                     typename internal::enable_if<d==2,U>::type * = NULL)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d
                     <<"D tensor B-spline using 2 knot-vectors.");
        
        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        Basis * tbasis = Basis::New(cbases); //d==2
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
        GISMO_ASSERT(tbasis->size()==m_coefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
    }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    template<typename U>
    gsTensorBSpline( gsKnotVector<T> const & KV1, gsKnotVector<U> const & KV2,
                     const gsMatrix<T> & tcoefs,
                     typename internal::enable_if<d==2,U>::type * = NULL)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 2 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        Basis * tbasis = Basis::New(cbases); //d==2
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and 4 corner vertices,
    /// the four boundary curves are linear interpolations of the corners,
    /// the size of the matrix *corner* is 4 by 3
    gsTensorBSpline(gsMatrix<T> const    & corner, 
                    gsKnotVector<T> const & KV1, 
                    gsKnotVector<T> const & KV2);

    /// Construct 3D tensor B-Spline by knot vectors, degrees and coefficient matrix
    gsTensorBSpline( gsKnotVector<T> const & KV1, 
                     gsKnotVector<T> const & KV2, 
                     gsKnotVector<T> const & KV3,
                     gsMovable< gsMatrix<T> > tcoefs)
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 3 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        cbases.push_back(new gsBSplineBasis<T>(KV3) );
        Basis * tbasis = Basis::New(cbases); //d==3
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
        
        GISMO_ASSERT(tbasis->size()==m_coefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
    }

    /// Construct 3D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( gsKnotVector<T> const & KV1, 
                     gsKnotVector<T> const & KV2, 
                     gsKnotVector<T> const & KV3,
                     const gsMatrix<T> & tcoefs )
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 3 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        cbases.push_back(new gsBSplineBasis<T>(KV3) );
        Basis * tbasis = Basis::New(cbases); //d==3
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }

    /// Construct 4D tensor B-Spline by knot vectors, degrees and coefficient matrix
    gsTensorBSpline( gsKnotVector<T> const & KV1, 
                     gsKnotVector<T> const & KV2, 
                     gsKnotVector<T> const & KV3,
                     gsKnotVector<T> const & KV4,
                     gsMovable< gsMatrix<T> > tcoefs )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 4 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.reserve(4);
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        cbases.push_back(new gsBSplineBasis<T>(KV3) );
        cbases.push_back(new gsBSplineBasis<T>(KV4) );
        Basis * tbasis = Basis::New(cbases); //d==4
            
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
        GISMO_ASSERT(tbasis->size()==m_coefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
    }

    /// Construct 4D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( gsKnotVector<T> const & KV1, 
                     gsKnotVector<T> const & KV2, 
                     gsKnotVector<T> const & KV3,
                     gsKnotVector<T> const & KV4,
                     const gsMatrix<T> & tcoefs )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 4 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.reserve(4);
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        cbases.push_back(new gsBSplineBasis<T>(KV3) );
        cbases.push_back(new gsBSplineBasis<T>(KV4) );
        Basis * tbasis = Basis::New(cbases); //d==4
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(), 
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = tcoefs;
    }
    
    /// Clone function. Used to make a copy of the geometry
    gsTensorBSpline * clone() const
    { 
        return new gsTensorBSpline( *this );
    }

    GISMO_BASIS_ACCESSORS

public:

    // Look at gsGeometry class for a description
    void degreeElevate(int const i = 1, int const dir = -1);

    /// Inserts knot \a knot at direction \a dir, \a i times
    void insertKnot( T knot, int dir, int i = 1);

    /// Returns a reference to the knot vector in direction \a i
    gsKnotVector<T> & knots(const int i) { return this->basis().knots(i); }

    /// Returns a reference to the knot vector \a i
    const gsKnotVector<T> & knots(const int i) const { return this->basis().knots(i); }

    /*** Virtual member functions required by the base class ***/

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    /*** Additional members for tensor B-Splines ***/

    /// Returns the degree of the basis wrt direction i
    unsigned degree(const unsigned & i) const 
    { return this->basis().component(i).degree(); }

    /// Toggle orientation wrt coordinate k
    /// \todo use flipTensor to generalize to any dimension
    void reverse(unsigned k);

    /// Toggle orientation wrt coordinate k
    /// \todo use flipTensor to generalize to any dimension
    void swapDirections(const unsigned i, const unsigned j);
    
    /// \brief Return true if point \a u is a corner of
    /// the patch with tolerance \a tol
    bool isPatchCorner(gsMatrix<T> const &v, T tol = 1e-3) const;

    /// \brief returns the tensor-index \a curr of the corner control
    /// point \a v, or an invalid index if the corner is not found
    /// within the tolerance \a tol
    void findCorner(const gsMatrix<T>   & v, 
                    gsVector<index_t,d> & curr,
                    T tol = 1e-3);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the origin of the parametrization of the patch. Assumes
    /// that \a v is either input is indeed a corner of this patch
    void setOriginCorner(gsMatrix<T> const &v);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the ending corner of the parametrization of the
    /// patch. Assumes that \a v is either input is indeed a corner of
    /// this patch
    void setFurthestCorner(gsMatrix<T> const &v);
    
    /// Sets the resulting BSpline to be periodic in direction \em dir.
    /// \param dir
    inline void setPeriodic( int dir )
    {
        this->m_coefs = this->basis().perCoefs( this->m_coefs, dir );
        this->basis().setPeriodic( dir );
    }
    
    /// Returns a local representation of the geometry in the cell
    /// containing the point \a u
    gsGeometry<T> * localRep(const gsMatrix<T> & u) const;

private:

    /// Helper function for the slice function
    /// selects the row of coefficients from coefficients of geo that are suitable
    /// for the isoparametric slice in \a dir_fixed with \a par.
    /// Note that geo has to have already C^0 continuity at \a par in direction \a dir.
    void constructCoefsForSlice(unsigned dir_fixed,T par,
                                const gsTensorBSpline<d,T> & geo,
                                gsMatrix<T>& result) const;

public:

    /// Constucts an isoparametric slice of this tensorBSpline by fixing
    /// \a par in direction \a dir_fixed. The resulting tensorBSpline has
    /// one less dimension and is given back in \a result.
    void slice(index_t dir_fixed,T par,BoundaryGeometryType & result) const;

protected:
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...

    using Base::m_basis;
    using Base::m_coefs;

}; // class gsBSpline

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSpline.hpp)
#endif
