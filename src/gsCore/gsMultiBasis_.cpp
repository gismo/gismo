
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiBasis.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsMultiBasis<real_t>;

TEMPLATE_INST bool gsMultiBasis<real_t>::
repairInterfaceFindElements<2>(const boundaryInterface & bi,
                               std::vector<unsigned> & refEltsFirst,
                               std::vector<unsigned> & refEltsSecond );

TEMPLATE_INST bool gsMultiBasis<real_t>::
repairInterfaceFindElements<3>(const boundaryInterface & bi,
                               std::vector<unsigned> & refEltsFirst,
                               std::vector<unsigned> & refEltsSecond );
}
