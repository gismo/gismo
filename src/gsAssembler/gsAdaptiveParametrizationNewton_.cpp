
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsComposedFunction.h>
#include <gsCore/gsMultiPatch.h>

#include <gsAssembler/gsAdaptiveParametrizationNewton.h>
#include <gsAssembler/gsAdaptiveParametrizationNewton.hpp>

namespace gismo
{
CLASS_TEMPLATE_INST gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased>;
CLASS_TEMPLATE_INST gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased>;
}
