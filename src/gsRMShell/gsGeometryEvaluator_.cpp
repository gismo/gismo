
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsBasis.h>
#include <gsRMShell/gsGeometryEvaluator.h>
#include <gsRMShell/gsGeometryEvaluator.hpp>


namespace gismo
{
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 1, 0>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 1, 1>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 1, 2>;

	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 2, 0>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 2, 1>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 2, -1>;

	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 3, 0>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 3, 1>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 3, -2>;

	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 4, 0>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 4, 1>;
	CLASS_TEMPLATE_INST gsGenericGeometryEvaluator<real_t, 4, -3>;
}
