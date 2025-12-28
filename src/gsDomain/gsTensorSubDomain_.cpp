#include <gsCore/gsTemplateTools.h>

#include <gsDomain/gsTensorSubDomain.h>
#include <gsDomain/gsTensorSubDomain.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTensorSubDomain<1,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<2,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<3,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<4,real_t>;

    CLASS_TEMPLATE_INST gsTensorSubDomainIterator<real_t,1>;
    CLASS_TEMPLATE_INST gsTensorSubDomainIterator<real_t,2>;
    CLASS_TEMPLATE_INST gsTensorSubDomainIterator<real_t,3>;
    CLASS_TEMPLATE_INST gsTensorSubDomainIterator<real_t,4>;
}