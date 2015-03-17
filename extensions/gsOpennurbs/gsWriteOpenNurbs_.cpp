
#include <gsCore/gsTemplateTools.h>

#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.hpp>

namespace gismo {

namespace extensions {

TEMPLATE_INST bool writeON_PlanarDomain( const gsPlanarDomain<real_t> & pd);

}// namespace extensions

}// namespace gismo
