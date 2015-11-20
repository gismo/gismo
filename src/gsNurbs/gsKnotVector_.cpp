#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp>

namespace gismo
{
    //CLASS_TEMPLATE_INST gsKnotVector<real_t>;
    template class gsKnotVector<real_t>;

    CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;
}
