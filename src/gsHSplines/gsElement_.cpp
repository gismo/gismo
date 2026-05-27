
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsElement.h>
#include <gsHSplines/gsElement.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsElement<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST
}