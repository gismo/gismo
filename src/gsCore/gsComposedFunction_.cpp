
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedFunction.h>
#include <gsCore/gsComposedFunction.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedFunction<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedFunction<real_t> >;


// XML dispatch registration (priorities: see gsCoreXmlRegistration.h)
GISMO_XML_REGISTER(gsFunction<real_t>, gsComposedFunction<real_t>, 120)

}
