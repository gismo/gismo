#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsFunctionExpr.hpp>


namespace gismo
{

CLASS_TEMPLATE_INST gsFunctionExpr<real_t> ;
CLASS_TEMPLATE_INST internal::gsXml< gsFunctionExpr<real_t> >;

}
