#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsConstantFunction.hpp>


namespace gismo
{

CLASS_TEMPLATE_INST gsConstantFunction<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsConstantFunction<real_t> >;

}
