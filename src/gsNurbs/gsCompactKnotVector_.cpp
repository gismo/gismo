#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsCompactKnotVector.hpp>

namespace gismo
{
    //CLASS_TEMPLATE_INST gsCompactKnotVector<real_t>;

    CLASS_TEMPLATE_INST internal::gsXml< gsCompactKnotVector<real_t> >;
}
