#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBSpline<real_t>;

CLASS_TEMPLATE_INST internal::gsXml<gsBSpline<real_t> >;

}
