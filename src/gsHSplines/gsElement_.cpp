
#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsElement.h>
#include <gsHSplines/gsElement.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsElement<1,real_t>;
    CLASS_TEMPLATE_INST gsElement<2,real_t>;
    CLASS_TEMPLATE_INST gsElement<3,real_t>;
}