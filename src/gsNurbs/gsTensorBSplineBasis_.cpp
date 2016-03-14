/* Symbol export for G+Smo shared object */

//#define gsTensorBSplineBasis_EXPORT

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsTensorBSplineBasis<2,real_t>;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<3,real_t>;
CLASS_TEMPLATE_INST gsTensorBSplineBasis<4,real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t> >;

}
