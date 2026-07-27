#include <gsCore/gsTemplateTools.h>

#include <gsOptimizer/gsOptimizer.h>
#include <gsCore/gsRegistry.h>
#include <gsCore/gsRegistry.hpp>

namespace gismo
{

// The one authoritative optimizer registry instance lives in this
// translation unit (compiled into libgismo); runtime-loaded module
// libraries resolve the exported gsRegistry::get() from here.
CLASS_TEMPLATE_INST gsRegistry<gsOptimizer<real_t> >;

}
