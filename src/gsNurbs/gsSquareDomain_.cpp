
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsSquareDomain.h>
#include <gsNurbs/gsSquareDomain.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsSquareDomain<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsSquareDomain<real_t> >;

}
