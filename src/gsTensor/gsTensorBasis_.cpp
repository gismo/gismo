
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorBasis.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST gsTensorBasis<D, real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

}
