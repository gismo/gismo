
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsHElementHelper.h>
#include <gsHSplines/gsHElementHelper.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsHElementHelper<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST
}