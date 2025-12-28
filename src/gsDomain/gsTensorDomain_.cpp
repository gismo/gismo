#include <gsCore/gsTemplateTools.h>

#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsTensorDomain.hpp>
#include <gsDomain/gsTensorDomainIterator.hpp>
#include <gsDomain/gsTensorDomainBoundaryIterator.hpp>

#include <gsDomain/gsTensorSubDomain.h>
#include <gsDomain/gsTensorSubDomain.hpp>

#include <vector>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTensorDomain<1,real_t>;
    CLASS_TEMPLATE_INST gsTensorDomain<2,real_t>;
    CLASS_TEMPLATE_INST gsTensorDomain<3,real_t>;
    CLASS_TEMPLATE_INST gsTensorDomain<4,real_t>;
    CLASS_TEMPLATE_INST gsTensorDomain<-1,real_t>;

    CLASS_TEMPLATE_INST gsTensorSubDomain<1,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<2,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<3,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<4,real_t>;
    CLASS_TEMPLATE_INST gsTensorSubDomain<-1,real_t>;

    CLASS_TEMPLATE_INST gsTensorDomainIterator<real_t,1>;
    CLASS_TEMPLATE_INST gsTensorDomainIterator<real_t,2>;
    CLASS_TEMPLATE_INST gsTensorDomainIterator<real_t,3>;
    CLASS_TEMPLATE_INST gsTensorDomainIterator<real_t,4>;
    CLASS_TEMPLATE_INST gsTensorDomainIterator<real_t,-1>;

    CLASS_TEMPLATE_INST gsTensorDomainBoundaryIterator<real_t,1,typename gsKnotVector<real_t>::const_uiterator>;
    CLASS_TEMPLATE_INST gsTensorDomainBoundaryIterator<real_t,2,typename gsKnotVector<real_t>::const_uiterator>;
    CLASS_TEMPLATE_INST gsTensorDomainBoundaryIterator<real_t,3,typename gsKnotVector<real_t>::const_uiterator>;
    CLASS_TEMPLATE_INST gsTensorDomainBoundaryIterator<real_t,4,typename gsKnotVector<real_t>::const_uiterator>;
    // This is needed for gsRemapInterface
    CLASS_TEMPLATE_INST gsTensorDomainBoundaryIterator<real_t,-1,std::vector<real_t>::const_iterator>;

}