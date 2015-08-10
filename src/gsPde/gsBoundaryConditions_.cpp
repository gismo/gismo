
#include <gsCore/gsTemplateTools.h>

#include <gsPde/gsBoundaryConditions.h>
#include <gsPde/gsBoundaryConditions.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBoundaryConditions<real_t>;
    CLASS_TEMPLATE_INST internal::gsXml< gsBoundaryConditions<real_t> >;
}

