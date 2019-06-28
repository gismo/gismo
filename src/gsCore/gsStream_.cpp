#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsStream.h>

namespace gismo
{

  CLASS_TEMPLATE_INST gsStream<gsJob_enum::Python> ;
  CLASS_TEMPLATE_INST gsStream<gsJob_enum::CXX> ;

  static gsStream<gsJob_enum::Python> _stream_python;
  static gsStream<gsJob_enum::CXX> _stream_cxx;

}
