#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsConstantFunction.hpp>
#include <gsIO/gsXmlRegistry.h>


namespace gismo
{

CLASS_TEMPLATE_INST gsConstantFunction<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsConstantFunction<real_t> >;


// XML dispatch registration (priorities: see gsCoreXmlRegistration.h)
GISMO_XML_REGISTER(gsFunction<real_t>, gsConstantFunction<real_t>, 100)

}
