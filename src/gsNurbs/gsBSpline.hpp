/** @file gsBSpline.hpp

    @brief Implementation of a B-spline curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once 

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

template<class T, class KnotVectorType >    
bool gsBSpline<T,KnotVectorType>::isOn(gsMatrix<T> const &u, T tol) const
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

template<class T, class KnotVectorType >    
bool gsBSpline<T,KnotVectorType>::isPatchCorner(gsMatrix<T> const &v, T tol) const
{
    return (( v - m_coefs.row(0)       ).squaredNorm() < tol ||
            ( v - m_coefs.bottomRows(1)).squaredNorm() < tol );
}

template<class T, class KnotVectorType >    
void gsBSpline<T,KnotVectorType>::findCorner(const gsMatrix<T> & v, 
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


template<class T, class KnotVectorType >    
void gsBSpline<T,KnotVectorType>::setOriginCorner(gsMatrix<T> const &v)
{
    if ((v - m_coefs.row(0)).squaredNorm() < 1e-3)
        return;
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < 1e-3)
        this->reverse();
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template<class T, class KnotVectorType >    
void gsBSpline<T,KnotVectorType>::setFurthestCorner(gsMatrix<T> const &v)
{
    if ((v - m_coefs.row(0)).squaredNorm() < 1e-3)
        this->reverse();
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < 1e-3)
        return;
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template <class T, class KnotVectorType>
void gsBSpline<T,KnotVectorType>::swapDirections(const unsigned i, const unsigned j)
{
    GISMO_ASSERT( static_cast<int>(i) == 0 && static_cast<int>(j) == 0,
                  "Invalid basis components "<<i<<" and "<<j<<" requested" );
}


namespace internal
{

/// Get a BSpline from XML data
template<class T>
class gsXml< gsBSpline<T, gsKnotVector<T> > > // TO DO: KnotVectorType
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBSpline<TMPLA2(T,gsKnotVector<T>)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "BSpline"; }

    static gsBSpline<T, gsKnotVector<T> > * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsBSpline<T, gsKnotVector<T> > >(node);
    }

    static gsXmlNode * put (const gsBSpline<T, gsKnotVector<T> > & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsBSpline<T, gsKnotVector<T> > >(obj,data);
    }
};


}// namespace internal

} // namespace gismo
