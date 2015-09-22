#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXml.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>

#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsHBSplineBasis.hpp>

#include <gsHSplines/gsHBSpline.h>
#include <gsHSplines/gsHBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsHTensorBasis <1,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <2,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <3,real_t>;
CLASS_TEMPLATE_INST gsHTensorBasis <4,real_t>;

CLASS_TEMPLATE_INST gsHBSplineBasis <1,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <2,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <3,real_t>;
CLASS_TEMPLATE_INST gsHBSplineBasis <4,real_t>;

namespace internal
{

CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<4,real_t> >;

CLASS_TEMPLATE_INST gsXml< gsHBSpline<1,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsHBSpline<4,real_t> >;
}

}
