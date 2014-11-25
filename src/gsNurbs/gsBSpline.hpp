
#pragma once 

#include <gsIO/gsXmlUtils.h>

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


namespace internal
{

/// Get a BSpline from XML data
template<class T>
class gsXml< gsBSpline<T, gsKnotVector<T> > > // TO DO: KnotVectorType
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBSpline<A2(T,gsKnotVector<T>)>);
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
