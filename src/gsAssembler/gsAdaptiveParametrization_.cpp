
#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveParametrization.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsOptSigma<real_t>;
TEMPLATE_INST real_t gsCheckSigmaGradient<real_t>(gsOptProblem<real_t> &, gsSquareDomain<real_t> &, real_t, real_t);
CLASS_TEMPLATE_INST gsOptMesh<real_t,MonitorMode::ValueBased>;
CLASS_TEMPLATE_INST gsOptMesh<real_t,MonitorMode::GradientBased>;
CLASS_TEMPLATE_INST gsOptFit<real_t>;
CLASS_TEMPLATE_INST gsOptL2<real_t>;

CLASS_TEMPLATE_INST gsAdaptiveParametrizationBase<real_t>;

CLASS_TEMPLATE_INST gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>;
CLASS_TEMPLATE_INST gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>;

// Member templates of an explicitly-instantiated class template are NOT
// instantiated by CLASS_TEMPLATE_INST above; every <d> a caller outside this
// translation unit uses (unittests, examples) needs its own line here, or
// the call is an undefined reference at link time (GISMO_BUILD_LIB hides the
// .hpp from every other TU -- see the #ifndef guard at the bottom of
// gsAdaptiveParametrization.h). Only d=2 is used anywhere in the tree today.
TEMPLATE_INST gsTensorBSplineBasis<2,real_t> gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::makeIntegrationBasis<2>(const gsTensorBSplineBasis<2,real_t> &, const gsTensorBSplineBasis<2,real_t> &);
TEMPLATE_INST gsTensorBSplineBasis<2,real_t> gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>::makeIntegrationBasis<2>(const gsTensorBSplineBasis<2,real_t> &, const gsTensorBSplineBasis<2,real_t> &);
TEMPLATE_INST gsBasis<real_t>::uPtr gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::makeIntegrationBasis<2>(const gsTHBSplineBasis<2,real_t> &, const gsTensorBSplineBasis<2,real_t> &);
TEMPLATE_INST gsBasis<real_t>::uPtr gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>::makeIntegrationBasis<2>(const gsTHBSplineBasis<2,real_t> &, const gsTensorBSplineBasis<2,real_t> &);
}
