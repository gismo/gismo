#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
// gsXml<gsNurbs>/<gsTensorNurbs> still live centrally in gsXmlUtils.hpp
// (relocated to this module in step A6); that header includes no class
// .hpp files, so it is safe to include outside the instantiation units.
#include <gsIO/gsXmlUtils.hpp>
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

}

#ifdef GISMO_PRIMARY_INSTANCE
// Anchor against linker dead-stripping in static-archive consumers
extern "C" GISMO_EXPORT void gismo_xml_anchor_gsNurbs(void) { }
#endif
