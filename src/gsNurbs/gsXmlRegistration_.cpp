#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsNurbs/gsNurbsXml.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{

// XML dispatch registration (priorities: see gsNurbsXmlRegistration.h)
GISMO_XML_REGISTER(gsGeometry<real_t>, gsNurbs<real_t>, 110)
GISMO_XML_REGISTER(gsCurve<real_t>,    gsNurbs<real_t>, 110)
GISMO_XML_REGISTER(gsGeometry<real_t>, gsTensorNurbs<TMPLA2(2,real_t)>, 150)
GISMO_XML_REGISTER(gsGeometry<real_t>, gsTensorNurbs<TMPLA2(3,real_t)>, 160)
GISMO_XML_REGISTER(gsGeometry<real_t>, gsTensorNurbs<TMPLA2(4,real_t)>, 170)
GISMO_XML_REGISTER(gsSurface<real_t>,  gsTensorNurbs<TMPLA2(2,real_t)>, 110)

// Basis family (specs also central until step A6)
GISMO_XML_REGISTER(gsBasis<real_t>, gsNurbsBasis<real_t>, 110)
GISMO_XML_REGISTER(gsBasis<real_t>, gsTensorNurbsBasis<TMPLA2(2,real_t)>, 150)
GISMO_XML_REGISTER(gsBasis<real_t>, gsTensorNurbsBasis<TMPLA2(3,real_t)>, 160)
GISMO_XML_REGISTER(gsBasis<real_t>, gsTensorNurbsBasis<TMPLA2(4,real_t)>, 170)


// instantiations of the specs defined in gsNurbsXml.hpp
namespace internal {
CLASS_TEMPLATE_INST gsXml< gsNurbs<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<4,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsNurbsBasis<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<2,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<3,real_t> >;
CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<4,real_t> >;
}

}

#ifdef GISMO_PRIMARY_INSTANCE
// Anchor against linker dead-stripping in static-archive consumers
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsNurbs(void) { }
#endif
