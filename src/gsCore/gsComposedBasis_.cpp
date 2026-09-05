
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedBasis.h>
#include <gsCore/gsComposedBasis.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedBasis<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedBasis<real_t> >;


// XML dispatch registration: last in the historical put chain
GISMO_XML_REGISTER(gsBasis<real_t>, gsComposedBasis<real_t>, 900)

}
