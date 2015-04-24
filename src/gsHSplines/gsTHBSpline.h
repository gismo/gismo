/** @file gsTHBSpline.h

    @brief Provides declaration of THBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/** \brief
    A truncated hierarchical B-Spline function, in \em d dimensions.

    This is the geometry type associated with gsTHBSplineBasis.

    R^d -> R

    \tparam T is the coefficient type

    \ingroup geometry
    \ingroup HSplines
*/

template<unsigned d, class T>
class gsTHBSpline : public gsGenericGeometry<gsTHBSplineBasis<d,T> >
{
public:
    typedef gsTHBSplineBasis<d,T> Basis;
    typedef gsGenericGeometry< gsTHBSplineBasis<d,T> > Base;

    /// Shared pointer for gsHBSpline
    typedef memory::shared_ptr< gsTHBSpline<d,T> > Ptr;

    typedef typename
    choose<d==1, gsConstantFunction<T>, gsTHBSpline<d-1,T>
           >::type BoundaryGeometryType;

    typedef typename gsTHBSplineBasis<d,T>::BoundaryBasisType BoundaryBasisType;

public:

    /// Default empty constructor
    gsTHBSpline() { }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
        Base( basis, coefs ) 
    {
        
    }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis& basis, gsMatrix<T> & coefs ) :
        Base( basis, coefs )
    {

    }

    /// Construct B-Spline from a Tensor B-Spline
    gsTHBSpline( const gsTensorBSpline<d,T> & tbsp )
    {
        this->m_basis = new Basis(tbsp.basis(), 3);// 3 levels
        this->m_coefs = tbsp.coefs();
    }

    /// Copy constructor
    gsTHBSpline( const gsTHBSpline & other )
    {
        this->m_basis = other.basis().clone();
        this->m_coefs = other.coefs();
    }

    /// Clone the gsHBspline
    virtual gsTHBSpline * clone() const
    { return new gsTHBSpline(*this); }

    ~gsTHBSpline() { } //destructor

    //void deriv2_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    /// Returns the degree wrt direction i
    unsigned degree(const unsigned & i) const
     //{ return this->basisComponent(i).degree(); };
     { return this->basis().component(i).degree(); }

    /** Refines the basis and adjusts the coefficients to keep the geometry the same.
     * The indices of the boxes are the same as in gsHTensorBasis<>::refineElements_withCoefs.
     */
    void refineElements( std::vector<unsigned> const & boxes )
    {
        gsMatrix<> & coefs = this->m_coefs;
        this->basis().refineElements_withCoefs( coefs, boxes );
    }


//////////////////////////////////////////////////
/// Other member functions
//////////////////////////////////////////////////
public:
    ///get all the B-spline patches out of a THB-spline geometry
    //void getBsplinePatches(gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level, std::vector< gsTensorBSpline<2> > & bpatches) const;
    void getBsplinePatches(gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level) const;

    /// Refines the whole domain to the finest level present in the mesh. Returns the refined geometry as result.
    void convertToBSpline( gsTensorBSpline<d,T,gsCompactKnotVector<T> >& result );

    /// Refines the whole domain to the finest level present in the mesh. Returns the refined geometry as result.
    void convertToBSpline( gsTensorBSpline<d,T>& result );

    /// Increase multiplicity of knot-value \a knotValue in level \a lvl and direction \a dir by \a mult
    void increaseMultiplicity(index_t lvl, int dir, T knotValue, int mult = 1);


private:

    ///get B-spline control points on a given box of a certain level by refining eveywhere
    void getBsplinePatchGlobal(gsVector<unsigned> b1, gsVector<unsigned> b2, unsigned l, gsTensorBSpline<2> geo) const;
    ///function for getBsplinePatchGlobal
    void globalRefinement(int level)const;

    ///initialization of cmatrix
    void initialize_cmatrix(int col, int c_level) const;

    ///convert the coefficient matrix mat in the given direction to a column of the control points matrix
    void return_cp_1D(const gsMatrix<T> & mat, int direction, gsMatrix<T>& cp)const;

public:

    /// Constucts an isoparametric slice of this THBSpline by fixing
    /// \a par in direction \a dir_fixed. The resulting THBSpline has
    /// one less dimension and is given back in \a result.
    void slice(index_t dir_fixed,T par,BoundaryGeometryType & result) const;

}; // class gsTHBSpline


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTHBSpline.hpp)
#endif
