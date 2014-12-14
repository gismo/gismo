#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsCompactKnotVector.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsBSplineBasis<real_t, gsKnotVector<real_t>        >;
    CLASS_TEMPLATE_INST gsBSplineBasis<real_t, gsCompactKnotVector<real_t> >;

    CLASS_TEMPLATE_INST internal::gsXml< gsBSplineBasis<real_t, gsKnotVector<real_t>        > >;
    CLASS_TEMPLATE_INST internal::gsXml< gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >;

}
