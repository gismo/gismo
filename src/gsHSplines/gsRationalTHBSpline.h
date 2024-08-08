/** @file gsRationalTHBSpline.h

    @brief Represents a rational truncated hierarchical B-Spline patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/
 
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsHSplines/gsRationalTHBSplineBasis.h>
#include <gsCore/gsForwardDeclarations.h>


namespace gismo
{

/** \brief 
    A rational truncated hierarchical B-Spline function
    of parametric dimension \em d, with arbitrary target
    dimension.

    This is the geometry type associated with gsRationalTHBSplineBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    \ingroup geometry
    \ingroup HSplines
*/

template<short_t d, class T>
class gsRationalTHBSpline : public gsGeoTraits<d,T>::GeometryBase
{

public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef T Scalar_t;
    
    typedef gsTensorBSplineBasis<d,T> TBasis;      // underlying tensor basis
    
    /// Family type
    typedef gsBSplineBasis<T>  Family_t;
    
    // rational version of tensor basis (basis for this geometry)
    typedef gsRationalTHBSplineBasis<d,T>   Basis;

    /// Associated boundary geometry type
    // typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::RatGeometry BoundaryGeometryType;
    // typedef gsRationalTHBSpline<static_cast<short_t>(d-1),T> BoundaryGeometryType;
    typedef typename
    util::conditional<d==1, gsConstantFunction<T>, gsRationalTHBSpline<static_cast<short_t>(d-1),T>
                      >::type BoundaryGeometryType;

    /// Associated boundary basis type
    typedef gsRationalTHBSplineBasis<d-1,T> BoundaryBasisType;

    /// Shared pointer for gsRationalTHBSpline
    typedef memory::shared_ptr< gsRationalTHBSpline > Ptr;

    /// Unique pointer for gsRationalTHBSpline
    typedef memory::unique_ptr< gsRationalTHBSpline > uPtr;

public:

    /// Default empty constructor
    gsRationalTHBSpline() : Base() { }

    gsRationalTHBSpline( const Basis & basis,  gsMatrix<T> coefs ) :
    Base( basis, give(coefs) ) { }

    // gsRationalTHBSpline( gsTensorNurbs<d> nurbs): Base(new gsRationalTHBSplineBasis<d>(nurbs.basis().source(), nurbs.weights() ), nurbs.coefs() ){ }
    /// Construct B-Spline from a Tensor B-Spline
    explicit gsRationalTHBSpline( const gsTensorNurbs<d,T> & nurbs )
    {
        // gsTHBSplineBasis<d> thbBasis( nurbs.basis().source().clone() );
        // this->m_basis = new Basis( *thbBasis.clone(), nurbs.weights() );
        this->m_basis = new Basis( new gsTHBSplineBasis<d>( nurbs.basis().source() ), nurbs.weights() );
        this->m_coefs = nurbs.coefs();
    }



    GISMO_BASIS_ACCESSORS
    
    public:

// ***********************************************
// Virtual member functions required by the base class
// ***********************************************

    GISMO_CLONE_FUNCTION(gsRationalTHBSpline)

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

// ***********************************************
// Additional members
// ***********************************************

    // /// Inserts knot \a knot at direction \a dir, \a i times
    // void insertKnot( T knot, int dir, int i = 1)
    // {
    //     GISMO_ASSERT( i>0, "multiplicity must be at least 1");
    //     GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
    //                   "Invalid basis component "<< dir <<" requested for degree elevation" );

    //     // Combine control points and weights to n+1 dimensional "control points"
    //     gsMatrix<T> projectiveCoefs = basis().projectiveCoefs(m_coefs, weights());

    //     // Dimension of physical space + 1 for the weights
    //     const index_t n = projectiveCoefs.cols();

    //     Basis & tbs = this->basis().source();
    //     gsVector<index_t,d> sz;
    //     tbs.size_cwise(sz); // Size of basis per (parametric) direction
        
    //     swapTensorDirection(0, dir, sz, projectiveCoefs  );

    //     const index_t nc = sz.template tail<static_cast<short_t>(d-1)>().prod();
    //     projectiveCoefs.resize( sz[0], n * nc );
        
    //     // Insert knot and update projective control points
    //     gsBoehm(tbs.knots(dir), projectiveCoefs  , knot, i, true );
    //     sz[0] = projectiveCoefs.rows();

    //     // New number of coefficients
    //     const index_t ncoef = sz.prod();
    //     projectiveCoefs.resize(ncoef, n );

    //     swapTensorDirection(0, dir, sz, projectiveCoefs  );

    //     // Unpack projective control points and update the coefs and weights of the original NURBS
    //     basis().setFromProjectiveCoefs(projectiveCoefs, m_coefs, weights() );
    // }

    /// Access to i-th weight
    T & weight(int i) const { return this->basis().weight(i); }

    /// Returns the NURBS weights
    const gsMatrix<T> & weights() const { return this->basis().weights(); }

    /// Returns the NURBS weights as non-const reference
    gsMatrix<T> & weights() { return this->basis().weights(); }
    
    // /// Returns the degree of the basis wrt direction i 
    // short_t degree(unsigned i) const
    // { return this->basis().source().component(i).degree(); }

protected:
    // todo: check function: check the coefficient number, degree, knot vector ...

    using gsGeometry<T>::m_coefs;
    using gsGeometry<T>::m_basis;

}; // class gsRationalTHBSpline


// ***********************************************
// ***********************************************


} // namespace gismo
