#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsHPatchRule.h>
#include <gsAssembler/gsHPatchRule.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>

namespace gismo
{

CLASS_TEMPLATE_INST gsHPatchRule<real_t>;
CLASS_TEMPLATE_INST gsHPatchRuleCache<real_t>;

// Explicit instantiation of the dimension-templated makeIntegrationBasis
TEMPLATE_INST gsTHBSplineBasis<1,real_t,true>
gsHPatchRule<real_t>::makeIntegrationBasis<1>(const gsTHBSplineBasis<1,real_t,true> &);
TEMPLATE_INST gsTHBSplineBasis<2,real_t,true>
gsHPatchRule<real_t>::makeIntegrationBasis<2>(const gsTHBSplineBasis<2,real_t,true> &);
TEMPLATE_INST gsTHBSplineBasis<3,real_t,true>
gsHPatchRule<real_t>::makeIntegrationBasis<3>(const gsTHBSplineBasis<3,real_t,true> &);
TEMPLATE_INST gsTHBSplineBasis<4,real_t,true>
gsHPatchRule<real_t>::makeIntegrationBasis<4>(const gsTHBSplineBasis<4,real_t,true> &);

}
