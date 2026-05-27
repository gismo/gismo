
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsHElementMarker.h>
#include <gsHSplines/gsHElementMarker.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsHElementMarker<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST
}