
namespace gismo
{

namespace internal
{


/// Get a HBSpline from XML data
template<class T, unsigned d>
class gsXml< gsHBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHBSpline<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "HBSpline"+to_string(d); }

    static gsHBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsHBSpline<d,T> >(node);
    }
    
    static gsXmlNode * put (const gsHBSpline<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsHBSpline<d,T> >(obj,data);
    }
};


} // end namespace internal

} // end namespace gismo
