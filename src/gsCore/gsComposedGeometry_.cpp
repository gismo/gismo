
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedGeometry.h>
#include <gsCore/gsComposedGeometry.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedGeometry<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedGeometry<real_t> >;

}
