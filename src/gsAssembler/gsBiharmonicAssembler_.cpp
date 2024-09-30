
#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsBiharmonicAssembler.h>
#include <gsAssembler/gsBiharmonicAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBiharmonicAssembler<real_t, gsVisitorBiharmonic<real_t> >;
}
