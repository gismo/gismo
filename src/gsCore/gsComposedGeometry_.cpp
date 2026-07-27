
#include <gsCore/gsTemplateTools.h>
#include <gsIO/gsXmlRegistry.h>

#include <gsCore/gsComposedGeometry.h>
#include <gsCore/gsComposedGeometry.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedGeometry<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedGeometry<real_t> >;


// XML dispatch registration: last in the historical put chain
GISMO_XML_REGISTER(gsGeometry<real_t>, gsComposedGeometry<real_t>, 900)

}
