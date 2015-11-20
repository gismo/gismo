#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp>

namespace gismo
{
    //CLASS_TEMPLATE_INST gsKnotVector<real_t>;

    #ifdef _MSC_VER
    #pragma warning(suppress : 4231)
    #endif
    template class gsKnotVector<real_t>;

    CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;
}
