#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsPde/gsPdeXml.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{

// Pde reading through the abstract base (writing was never supported)
GISMO_XML_REGISTER_GET(gsPde<real_t>, gsPoissonPde<real_t>)


// instantiations of the specs defined in gsPdeXml.hpp
namespace internal {
CLASS_TEMPLATE_INST gsXml< gsPoissonPde<real_t> >;
}

}

#ifdef GISMO_PRIMARY_INSTANCE
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsPde(void) { }
#endif
