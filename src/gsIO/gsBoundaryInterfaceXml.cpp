/** @file gsBoundaryInterfaceXml.cpp

    @brief XML serialization for boundaryInterface (put/tag)

    Minimal gsXml specialization for boundaryInterface providing
    tag() and put() so gsFileData can write interfaces into XML.

*/

#include <sstream>
#include <gsCore/gsBoundary.h>
#include <gsIO/gsXml.h>

namespace gismo {
namespace internal {

// Provide XML tag for boundaryInterface
template<>
std::string gsXml<boundaryInterface>::tag()
{
    return std::string("interface");
}

// Minimal type string (left empty)
template<>
std::string gsXml<boundaryInterface>::type()
{
    return std::string("");
}

// Serialize a boundaryInterface into an XML node
template<>
gsXmlNode * gsXml<boundaryInterface>::put(const boundaryInterface & obj, gsXmlTree & data)
{
    std::ostringstream oss;
    // write patch/side pairs
    oss << obj.first().patch << " " << obj.first().index() << " "
        << obj.second().patch << " " << obj.second().index() << " ";

    // write direction map
    const gsVector<index_t> & dm = obj.dirMap();
    for (index_t i = 0; i < dm.rows(); ++i)
        oss << dm(i) << " ";

    // write orientation flags
    const gsVector<bool> & orient = obj.dirOrientation();
    for (index_t i = 0; i < orient.rows(); ++i)
        oss << (orient(i) ? 1 : 0) << " ";

    std::string content = oss.str();
    gsXmlNode * node = makeNode(tag(), content, data);
    if (!obj.label().empty())
        node->append_attribute(makeAttribute("name", obj.label(), data));
    return node;
}

} // namespace internal
} // namespace gismo
