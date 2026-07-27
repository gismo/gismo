#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
// gsXml<gsTrimSurface> still lives centrally in gsXmlUtils.hpp
// (relocated to this module in step A6)
#include <gsIO/gsXmlUtils.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{

// historically writable through gsGeometry but not readable
GISMO_XML_REGISTER_PUT(gsGeometry<real_t>, gsTrimSurface<real_t>, 210)

}

#ifdef GISMO_PRIMARY_INSTANCE
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsModeling(void) { }
#endif
