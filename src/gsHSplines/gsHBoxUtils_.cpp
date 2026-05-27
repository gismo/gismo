#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsCore/gsForwardDeclarations.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>

#include <gsHSplines/gsHBoxUtils.h>
#include <gsHSplines/gsHBoxUtils.hpp>

namespace gismo
{

#define INST(D) STRUCT_TEMPLATE_INST gsHBoxUtils<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) STRUCT_TEMPLATE_INST gsHBoxCompare<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) STRUCT_TEMPLATE_INST gsHBoxEqual<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) STRUCT_TEMPLATE_INST gsHBoxContains<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

#define INST(D) STRUCT_TEMPLATE_INST gsHBoxIsContained<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

}
