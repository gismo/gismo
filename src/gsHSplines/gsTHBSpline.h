/** @file gsTHBSpline.h

    @brief Provides declaration of THBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
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

template<short_t d, class T>
class gsTHBSpline : public gsGeoTraits<d,T>::GeometryBase
{
public:
    typedef gsTHBSplineBasis<d,T> Basis;

    typedef typename Basis::tensorBasis tensorBasis;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    /// Shared pointer for gsTHBSpline
    typedef memory::shared_ptr< gsTHBSpline > Ptr;

    /// Unique pointer for gsTHBSpline
    typedef memory::unique_ptr< gsTHBSpline > uPtr;

    typedef typename
    util::conditional<d==1, gsConstantFunction<T>, gsTHBSpline<static_cast<short_t>(d-1),T>
                      >::type BoundaryGeometryType;

    typedef typename gsTHBSplineBasis<d,T>::BoundaryBasisType BoundaryBasisType;

public:

    /// Default empty constructor
    gsTHBSpline() { }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs ) 
    { }

    /// Construct THB-Spline by basis functions and coefficient matrix
    gsTHBSpline( const Basis& basis, gsMatrix<T> & coefs ) :
        Base( basis, coefs )
    { }

    /// Construct B-Spline from a Tensor B-Spline
    explicit gsTHBSpline( const gsTensorBSpline<d,T> & tbsp )
    {
        this->m_basis = new Basis( tbsp.basis() );
        this->m_coefs = tbsp.coefs();
    }

    GISMO_CLONE_FUNCTION(gsTHBSpline)

    GISMO_BASIS_ACCESSORS

public:


// ************************************************
//  Other member functions
// ************************************************
public:
    /** \brief Return the list of B-spline patches to represent a THB-spline geometry

    \param[out] b1 bottom left corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
    \param[out] b2 top right corners of the box (vector of indices with respect to the gsCompactKnotVector of the highest possible level)
    \param[out] level levels of the boxes (level[i]: level of the i-th box,)
    */
    void getBsplinePatches(gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level) const;
    //\param[out] bpatches list of B-spline patches associated with the boxes
    //void getBsplinePatches(gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level, std::vector< gsTensorBSpline<2> > & bpatches) const;

    /// Refines the whole domain to the finest level present in the mesh. Returns the refined geometry as result.
//    void convertToBSpline( gsTensorBSpline<d,T,gsCompactKnotVector<T> >& result );

    /// Refines the whole domain to the finest level present in the mesh. Returns the refined geometry as result.
    void convertToBSpline( gsTensorBSpline<d,T>& result );

    /// Increase multiplicity of knot-value \a knotValue in level \a lvl and direction \a dir by \a mult
    void increaseMultiplicity(index_t lvl, int dir, T knotValue, int mult = 1);


private:

    ///get B-spline control points on a given box of a certain level by refining eveywhere
    void getBsplinePatchGlobal(gsVector<index_t> b1, gsVector<index_t> b2, unsigned l, gsTensorBSpline<2> geo) const;
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


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTHBSpline.hpp)
#endif
