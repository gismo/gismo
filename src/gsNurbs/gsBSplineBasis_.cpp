
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.hpp>

//#include <misc/gsInstanceTools.h>
//#include <gsNurbs/gsBSplineBasis.tpl>

namespace gismo
{

CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t, gsKnotVector<real_t> >;

CLASS_TEMPLATE_INST gsBSplineBasis<real_t, gsKnotVector<real_t>        >;

CLASS_TEMPLATE_INST internal::gsXml< gsBSplineBasis<real_t, gsKnotVector<real_t>        > >;

}
