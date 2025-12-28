#include <gsCore/gsTemplateTools.h>

#include <gsDomain/gsHDomain.h>
#include <gsDomain/gsHDomain.hpp>
#include <gsDomain/gsHDomainIterator.hpp>
#include <gsDomain/gsHDomainBoundaryIterator.hpp>
#include <gsDomain/gsCompositeDomain.h>
#include <gsDomain/gsHSubDomain.h>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHDomain<1,real_t,index_t>;
    CLASS_TEMPLATE_INST gsHDomain<2,real_t,index_t>;
    CLASS_TEMPLATE_INST gsHDomain<3,real_t,index_t>;
    CLASS_TEMPLATE_INST gsHDomain<4,real_t,index_t>;
    CLASS_TEMPLATE_INST gsHDomainIterator<real_t,1,index_t>;
    CLASS_TEMPLATE_INST gsHDomainIterator<real_t,2,index_t>;
    CLASS_TEMPLATE_INST gsHDomainIterator<real_t,3,index_t>;
    CLASS_TEMPLATE_INST gsHDomainIterator<real_t,4,index_t>;
    CLASS_TEMPLATE_INST gsHDomainBoundaryIterator<real_t,1,index_t>;
    CLASS_TEMPLATE_INST gsHDomainBoundaryIterator<real_t,2,index_t>;
    CLASS_TEMPLATE_INST gsHDomainBoundaryIterator<real_t,3,index_t>;
    CLASS_TEMPLATE_INST gsHDomainBoundaryIterator<real_t,4,index_t>;
}
