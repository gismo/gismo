#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsModeling/gsModelingXml.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{

// historically writable through gsGeometry but not readable
GISMO_XML_REGISTER_PUT(gsGeometry<real_t>, gsTrimSurface<real_t>, 210)


// instantiations of the specs defined in gsModelingXml.hpp
namespace internal {
CLASS_TEMPLATE_INST gsXml< gsSolid<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTrimSurface<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsPlanarDomain<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsCurveLoop<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsCurveFitting<real_t> >;
}

}

#ifdef GISMO_PRIMARY_INSTANCE
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsModeling(void) { }
#endif
