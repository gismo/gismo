/* Symbol export for G+Smo shared object */

#define gsBSplineBasis_EXPORT

#include <gsIO/gsXml.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.hpp> //dependency (otherwise already included)

namespace gismo
{

// CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t>;
// CLASS_TEMPLATE_INST gsBSplineBasis<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsBSplineBasis<real_t> >;

}
