
#include <gsIO/gsXmlUtils.h>

namespace gismo
{

namespace internal
{

/// Get a KnotVector from XML data
///
/// \ingroup Nurbs
template<class T>
class gsXml< gsCompactKnotVector<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsCompactKnotVector<T>);
    static std::string tag () { return "KnotVector"; }
    static std::string type() { return ""; } // "Compact" ?
    
    static gsCompactKnotVector<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ! strcmp( node->name(), "KnotVector"), "Invalid tag" );
        // && node->first_attribute("type")->value() is "Plain" or "Compact");
        
        int p = atoi(node->first_attribute("degree")->value() );
        gsCompactKnotVector<T> * kv = new gsCompactKnotVector<T>(p);        
        
        std::istringstream str;
        str.str( node->value() );
        for (T knot; str >> knot;) 
            kv->push_back(knot);
        
        return kv;
    }
    
    static gsXmlNode * put (const gsCompactKnotVector<T> & obj, gsXmlTree & data)
    {
        // Write the knot values (for now WITH multiplicities)            
        std::ostringstream str;
        str << std::setprecision(FILE_PRECISION);

        for ( typename gsCompactKnotVector<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            str << *it <<" ";
        }
        
        // Make a new XML KnotVector node 
        gsXmlNode * tmp = internal::makeNode("KnotVector", str.str(), data);
        // Append the degree attribure
        str.str(std::string());// clean the ostream
        str<< obj.degree();
        tmp->append_attribute( makeAttribute("degree", str.str(),data) );
        
        // todo : append "Compact" and write out as compact ?
        
        return tmp;
    }
};

}// namespace internal

}// namespace gismo
 
