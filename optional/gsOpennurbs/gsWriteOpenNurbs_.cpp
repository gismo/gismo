
#include <gsCore/gsTemplateTools.h>

#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.hpp>

namespace gismo {

namespace extensions {

TEMPLATE_INST bool writeON_PlanarDomain( const gsPlanarDomain<real_t> &, const std::string &);

TEMPLATE_INST bool writeON_MultiPatch  ( const gsMultiPatch<real_t> &, const std::string &);

TEMPLATE_INST bool writeON_NurbsSurface( const gsSurface<real_t> & srf, const std::string & name);

TEMPLATE_INST bool writeON_NurbsCurve( const gsCurve<real_t> & curve, const std::string & name);

TEMPLATE_INST bool writeON_Mesh ( const gsMesh<real_t> &, const std::string &);

}// namespace extensions

}// namespace gismo
