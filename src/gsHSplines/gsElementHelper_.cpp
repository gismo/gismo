
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsElementHelper.h>
#include <gsHSplines/gsElementHelper.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsElementHelper<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST
}