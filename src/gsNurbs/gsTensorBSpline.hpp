
#pragma once 

#include <gsIO/gsXmlUtils.h>

#include <gsNurbs/gsBSplineBasis.h>

namespace gismo
{

namespace internal
{

/// @brief Get a Tensor BSpline from XML data
///
/// \ingroup Nurbs
template<unsigned d, class T, class KnotVectorType>
class gsXml< gsTensorBSpline<d,T, KnotVectorType> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorBSpline<A3(d,T,KnotVectorType)>);
    static std::string tag ()  { return "Geometry"; }
    static std::string type () { return "TensorBSpline" +  to_string(d); }

    static gsTensorBSpline<d,T,KnotVectorType> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTensorBSpline<d,T,KnotVectorType> >( node );
    }
    
    static gsXmlNode * put (const gsTensorBSpline<d,T,KnotVectorType> & obj,
                            gsXmlTree & data)
    {
        return putGeometryToXml(obj,data);
    }
};



}// namespace internal

} // namespace gismo
