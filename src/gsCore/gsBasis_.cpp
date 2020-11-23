
#include <gsCore/gsDebug.h> // to purge warning on MSVC
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsBasis.h>
#include <gsCore/gsBasis.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBasis<real_t>;
}
