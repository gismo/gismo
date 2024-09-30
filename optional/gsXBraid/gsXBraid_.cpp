
#include <gsCore/gsTemplateTools.h>
#include <gsXBraid/gsXBraid.h>
#include <gsXBraid/gsXBraid.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsXBraid<real_t>;
  CLASS_TEMPLATE_INST gsXBraid< gsMatrix<real_t> >;
  CLASS_TEMPLATE_INST gsXBraid< gsVector<real_t> >;

}
