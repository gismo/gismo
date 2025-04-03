
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedFunction.h>
#include <gsCore/gsComposedFunction.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedFunction<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedFunction<real_t> >;

}
