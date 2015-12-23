#include <gsCore/gsTemplateTools.h>

#ifdef _MSC_VER
//non-standard extension: extern before explicit inst.
#pragma warning(push)
#pragma warning(disable : 4231)
#endif

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsKnotVector<real_t>;

    CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
