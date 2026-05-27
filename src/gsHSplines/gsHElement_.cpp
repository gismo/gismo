
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElement.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST gsHElement<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST
}