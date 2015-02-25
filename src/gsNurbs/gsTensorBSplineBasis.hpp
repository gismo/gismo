
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsNurbs/gsKnotVector.h>

#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlUtils.hpp> // !


namespace gismo
{


namespace internal
{

/// @brief Get a TensorBSplineBasis from XML data
template<unsigned d, class T, class KnotVectorType>
class gsXml< gsTensorBSplineBasis<d,T,KnotVectorType> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorBSplineBasis<TMPLA3(d,T,KnotVectorType)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "TensorBSplineBasis"+to_string(d); }

    static gsTensorBSplineBasis<d,T,KnotVectorType> * get (gsXmlNode * node)
    {
        return getTensorBasisFromXml<gsTensorBSplineBasis<d,T,KnotVectorType> >( node );
    }
    
    static gsXmlNode * put (const gsTensorBSplineBasis<d,T,KnotVectorType> & obj, 
                            gsXmlTree & data )
    {
        return putTensorBasisToXml<gsTensorBSplineBasis<d,T,KnotVectorType> >(obj,data);
    }
};

} // internal

} // gismo

