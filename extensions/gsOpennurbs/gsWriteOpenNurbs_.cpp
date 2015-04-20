
#include <gsCore/gsTemplateTools.h>

#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.hpp>

namespace gismo {

namespace extensions {

TEMPLATE_INST bool writeON_PlanarDomain( const gsPlanarDomain<real_t> & pd);

TEMPLATE_INST bool writeON_MultiPatch  ( const gsMultiPatch<real_t> & patches);

}// namespace extensions

}// namespace gismo
