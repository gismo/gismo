#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsJob.h>

namespace gismo
{

  CLASS_TEMPLATE_INST gsJob<gsJob_enum::Python> ;
  CLASS_TEMPLATE_INST gsJob<gsJob_enum::CXX> ;

}
