
#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsHElementMarker.h>
#include <gsHSplines/gsHElementMarker.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHElementMarker<2,real_t>;
    CLASS_TEMPLATE_INST gsHElementMarker<3,real_t>;
}