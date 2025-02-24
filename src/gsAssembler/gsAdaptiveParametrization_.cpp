
#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveParametrization.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsOptMesh<real_t,MonitorMode::ValueBased>;
CLASS_TEMPLATE_INST gsOptMesh<real_t,MonitorMode::GradientBased>;

CLASS_TEMPLATE_INST gsAdaptiveParametrizationBase<real_t>;

CLASS_TEMPLATE_INST gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>;
CLASS_TEMPLATE_INST gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>;
}
