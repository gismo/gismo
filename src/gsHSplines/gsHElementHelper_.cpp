
#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsHElementHelper.h>
#include <gsHSplines/gsHElementHelper.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHElementHelper<2,real_t>;
    CLASS_TEMPLATE_INST gsHElementHelper<3,real_t>;
}