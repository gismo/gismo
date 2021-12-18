/** @file gsHBSpline.h

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
#include <gsHSplines/gsHBSplineBasis.h>

namespace gismo
{

/** \brief
    A hierarchical B-Spline function, in \em d dimensions.

    This is the geometry type associated with gsHBSplineBasis.

    \tparam T is the coefficient type

    \ingroup geometry
    \ingroup HSplines
*/

template<short_t d, class T>
class gsHBSpline : public gsGeoTraits<d,T>::GeometryBase
{
public:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsHBSplineBasis<d,T> Basis;

    /// Shared pointer for gsHBSpline
    typedef memory::shared_ptr< gsHBSpline > Ptr;

    /// Unique pointer for gsHBSpline
    typedef memory::unique_ptr< gsHBSpline > uPtr;

    typedef typename
    util::conditional<d==1, gsConstantFunction<T>, gsHBSpline<static_cast<short_t>(d-1),T>
                      >::type BoundaryGeometryType;

    typedef typename gsHBSplineBasis<d,T>::BoundaryBasisType BoundaryBasisType;

public:

    /// Default empty constructor
    gsHBSpline() { }

    /// Construct HB-Spline by basis functions and coefficient matrix
    gsHBSpline( const Basis * basis, const gsMatrix<T> * coefs ) :
    Base( basis, coefs ) { }

    /// Construct HB-Spline by basis functions and coefficient matrix
    gsHBSpline( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs ) { }

    /// Construct B-Spline from a Tensor B-Spline
    explicit gsHBSpline( const gsTensorBSpline<d,T> & tbsp )
    {
        this->m_basis = new Basis( tbsp.basis() );
        this->m_coefs = tbsp.coefs();
    }

    GISMO_CLONE_FUNCTION(gsHBSpline)

    GISMO_BASIS_ACCESSORS

public:

    /// Constucts an isoparametric slice of this HBSpline by fixing
    /// \a par in direction \a dir_fixed. The resulting HBSpline has
    /// one less dimension and is given back in \a result.
    void slice(index_t dir_fixed,T par,BoundaryGeometryType & result) const
    {
        GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
        GISMO_ASSERT(dir_fixed>=0 && static_cast<unsigned>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
        const BoundaryBasisType * bBasis = this->basis().basisSlice(dir_fixed,par);
        if(d==1)
        {
            gsMatrix<T> val(1,1),point;
            val(0,0)=par;
            this->eval_into(val,point);
            result=BoundaryGeometryType(*bBasis,point);
        }
        else
        {
            gsMatrix<T> vals,anchorsSlice,anchorsInGeom;
            bBasis->anchors_into(anchorsSlice);
            anchorsInGeom.resize(anchorsSlice.rows()+1,anchorsSlice.cols());
            anchorsInGeom.topRows(dir_fixed)=anchorsSlice.topRows(dir_fixed);
            anchorsInGeom.row(dir_fixed)=gsVector<T>::Constant(anchorsSlice.cols(),par);
            anchorsInGeom.bottomRows(anchorsSlice.rows()-dir_fixed)=anchorsSlice.bottomRows(anchorsSlice.rows()-dir_fixed);
            this->eval_into(anchorsInGeom,vals);
            BoundaryGeometryType* geom =
                    dynamic_cast<BoundaryGeometryType *>(bBasis->interpolateAtAnchors(vals).release()); //todo make it better
            GISMO_ASSERT(geom!=NULL,"bBasis should have BoundaryGeometryType.");
            result = *geom;
            delete geom;
        }
        delete bBasis;
    }
}; // class gsHBSpline


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHBSpline.hpp)
#endif
