
#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsFieldDistance.h>
#include <gsAssembler/gsFieldDistance.hpp>

namespace gismo
{

TEMPLATE_INST
real_t distanceL2(const gsField<real_t> &, const gsField<real_t> &);
TEMPLATE_INST
real_t distanceL2(const gsField<real_t> &, const gsFunctionSet<real_t> &, bool);
TEMPLATE_INST
real_t distanceL2(const gsField<real_t> &, const gsFunctionSet<real_t> &,
                  const gsMultiBasis<real_t> &, bool);
TEMPLATE_INST
real_t distanceH1(const gsField<real_t> &, const gsFunctionSet<real_t> &, bool);
TEMPLATE_INST
real_t distanceH1(const gsField<real_t> &, const gsFunctionSet<real_t> &,
                  const gsMultiBasis<real_t> &, bool);
TEMPLATE_INST
real_t distanceH2(const gsField<real_t> &, const gsFunctionSet<real_t> &, bool);
TEMPLATE_INST
real_t distanceDG(const gsField<real_t> &, const gsFunctionSet<real_t> &, bool);

} // namespace gismo
