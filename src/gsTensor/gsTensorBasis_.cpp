
#include <gsCore/gsTemplateTools.h>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorBasis.hpp>

namespace gismo
{

//TEMPLATE_INST gsXmlNode * putTensorBasisToXml 

CLASS_TEMPLATE_INST gsTensorBasis<1, real_t  >;
CLASS_TEMPLATE_INST gsTensorBasis<2, real_t  >;
CLASS_TEMPLATE_INST gsTensorBasis<3, real_t  >;
CLASS_TEMPLATE_INST gsTensorBasis<4, real_t  >;

}
