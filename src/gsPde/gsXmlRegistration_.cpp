#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
// gsXml<gsPoissonPde>/<gsSurfacePoissonPde> still live centrally in
// gsXmlUtils.hpp (relocated to this module in step A6)
#include <gsIO/gsXmlUtils.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{

// Pde reading through the abstract base (writing was never supported)
GISMO_XML_REGISTER_GET(gsPde<real_t>, gsPoissonPde<real_t>)
GISMO_XML_REGISTER_GET(gsPde<real_t>, gsSurfacePoissonPde<real_t>)

}

#ifdef GISMO_PRIMARY_INSTANCE
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsPde(void) { }
#endif
