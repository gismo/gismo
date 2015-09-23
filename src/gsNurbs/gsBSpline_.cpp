#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsCompactKnotVector.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSpline.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBSpline<real_t, gsKnotVector<real_t>        >;
CLASS_TEMPLATE_INST gsBSpline<real_t, gsCompactKnotVector<real_t> >;

CLASS_TEMPLATE_INST internal::gsXml< gsBSpline<real_t, gsKnotVector<real_t>        > >;
// to do
//CLASS_TEMPLATE_INST internal::gsXml< gsBSpline<real_t, gsCompactKnotVector<real_t> > >;

}
