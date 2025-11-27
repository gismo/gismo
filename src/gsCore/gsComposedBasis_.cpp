
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedBasis.h>
#include <gsCore/gsComposedBasis.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsComposedBasis<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsComposedBasis<real_t> >;

}
