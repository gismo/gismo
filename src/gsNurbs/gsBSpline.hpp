/** @file gsBSpline.hpp

    @brief Implementation of a B-spline curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once 

#include <gsNurbs/gsBSplineAlgorithms.h>
#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{


template<class T>
void gsBSpline<T>::merge( gsGeometry<T> * otherG )
{
    // See also gsNurbs::merge().
    // check geometric dimension
    GISMO_ASSERT(this->geoDim()==otherG->geoDim(),
                 "gsNurbs: cannot merge curves in different spaces ( R^"
                 << this->geoDim() << ", R^" << otherG->geoDim() << " ).");

    // check if the type of other is BSpline
    gsBSpline *  other = dynamic_cast<gsBSpline *>( otherG );
    GISMO_ASSERT( other!=NULL, "Can only merge with B-spline curves.");
    other= other->clone().release();

    GISMO_ASSERT( this ->basis().isPeriodic() == false &&
                  other->basis().isPeriodic() == false,
                  "Cannot merge a closed curve with anything." );

    // check degree
    const int mDeg = this ->basis().degree();
    const int oDeg = other->basis().degree();
    const int deg  = math::max(mDeg,oDeg);

    other->gsBSpline::degreeElevate( deg - oDeg ); // degreeElevate(0) does nothing (and very quickly)
    this ->gsBSpline::degreeElevate( deg - mDeg );

    // check whether the resulting curve will be continuous
    // TODO: ideally, the tolerance should be a parameter of the function
    T tol = 1e-8;
    gsMatrix<T> mValue = this ->eval(this ->support().col(1));
    gsMatrix<T> oValue = other->eval(other->support().col(0));
    bool continuous = gsAllCloseAbsolute(mValue,oValue,tol);

    // merge knot vectors.
    KnotVectorType& mKnots = this ->basis().knots();
    KnotVectorType& oKnots = other->basis().knots();
    T lastKnot = mKnots.last();
    if (continuous) // reduce knot multiplicity
    {
        // TODO check for clamped knot vectors otherwise
        // we should do knot insertion beforehands
        mKnots.remove(lastKnot);
    }// else there is a knot of multiplicity deg + 1 at the discontinuity

    oKnots.addConstant(lastKnot-oKnots.first());
    mKnots.append( oKnots.begin()+deg+1, oKnots.end());

    // merge coefficients
    int n= this->coefsSize();
    int skip = continuous ? 1 : 0;
    this->m_coefs.conservativeResize( n + other->coefsSize() -skip, Eigen::NoChange ) ;

    this->m_coefs.block( n,0,other->coefsSize()-skip,other->geoDim() ) =
        other->m_coefs.block( 1,0,other->coefsSize()-skip,other->geoDim() ) ;

    delete other;
}

template<class T>    
bool gsBSpline<T>::isOn(gsMatrix<T> const &u, T tol) const
{
    GISMO_ASSERT( u.cols() == 1, "Expecting single point.");
    gsBSplineSolver<T> slv;
    gsMatrix<T> e;
    
    for ( index_t k = 0; k != u.rows(); ++k)
    {        
        std::vector<T> roots;
        slv.allRoots(*this, roots, k, u(k,0) ) ;

        if( roots.size()!=0 )
        {
            gsAsMatrix<T> xx(roots); 
            this->eval_into( xx, e ); 
            
            //if( math::abs( e(!k,0)-u(!k,0) ) < tol )
            for(index_t j=0; j!=e.cols(); j++)
                if( ( e.col(j)-u ).norm() < tol )
                    return true;
        }
    }
    return false;
}

template<class T>    
bool gsBSpline<T>::isPatchCorner(gsMatrix<T> const &v, T tol) const
{
    return (( v - m_coefs.row(0)       ).squaredNorm() < tol ||
            ( v - m_coefs.bottomRows(1)).squaredNorm() < tol );
}

template<class T>    
void gsBSpline<T>::findCorner(const gsMatrix<T> & v, 
                                             gsVector<index_t,1> & curr,
                                             T tol)
{
    if ((v - m_coefs.row(0)).squaredNorm() < tol)
        curr[0] = 0;
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < tol)
        curr[0] = m_coefs.rows()-1;
    else
    {
        curr[0] = m_coefs.rows(); // invalidate result
        gsWarn<<"Point "<< v <<" is not an corner of the patch. (Call isPatchCorner() first!).\n";
    }
}


template<class T>    
void gsBSpline<T>::setOriginCorner(gsMatrix<T> const &v)
{
    if ((v - m_coefs.row(0)).squaredNorm() < 1e-3)
        return;
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < 1e-3)
        this->reverse();
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template<class T>    
void gsBSpline<T>::setFurthestCorner(gsMatrix<T> const &v)
{
    if ((v - m_coefs.bottomRows(1)).squaredNorm() < 1e-3)
        return;
    else if ((v - m_coefs.row(0)).squaredNorm() < 1e-3)
        this->reverse();
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template <class T>
void gsBSpline<T>::swapDirections(const unsigned i, const unsigned j)
{
    GISMO_UNUSED(i);
    GISMO_UNUSED(j);
    GISMO_ASSERT( static_cast<int>(i) == 0 && static_cast<int>(j) == 0,
                  "Invalid basis components "<<i<<" and "<<j<<" requested" );
}


template<class T>
void gsBSpline<T>::degreeElevate(short_t const i, short_t const dir)
{
    GISMO_UNUSED(dir);
    GISMO_ASSERT( (dir == -1) || (dir == 0),
                  "Invalid basis component "<< dir <<" requested for degree elevation" );
    
    bspline::degreeElevateBSpline(this->basis(), this->m_coefs, i);
}

namespace internal
{

/// Get a BSpline from XML data
template<class T>
class gsXml< gsBSpline<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBSpline<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "BSpline"; }

    static gsBSpline<T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsBSpline<T> >(node);
    }

    static gsXmlNode * put (const gsBSpline<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsBSpline<T> >(obj,data);
    }
};


}// namespace internal

} // namespace gismo
