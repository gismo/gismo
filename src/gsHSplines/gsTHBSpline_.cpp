#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.hpp>

#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsTHBSplineBasis <1, real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <2,real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <3,real_t>;
CLASS_TEMPLATE_INST gsTHBSplineBasis <4,real_t>;

CLASS_TEMPLATE_INST gsTHBSpline      <1,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <2,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <3,real_t>;
CLASS_TEMPLATE_INST gsTHBSpline      <4,real_t>;

namespace internal
{
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsTHBSpline<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTHBSpline<4,real_t> >;
}

}
