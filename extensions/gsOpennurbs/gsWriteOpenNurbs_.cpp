
#include <gsCore/gsTemplateTools.h>

#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.hpp>

namespace gismo {

namespace extensions {

TEMPLATE_INST bool writeON_PlanarDomain( const gsPlanarDomain<real_t> &);

TEMPLATE_INST bool writeON_MultiPatch  ( const gsMultiPatch<real_t> &);

TEMPLATE_INST bool writeON_Mesh ( const gsMesh<real_t> &, const std::string &);

}// namespace extensions

}// namespace gismo
