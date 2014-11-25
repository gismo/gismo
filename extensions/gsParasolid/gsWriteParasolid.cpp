
#include <gsCore/gsTemplateTools.h>

#include <gsParasolid/gsWriteParasolid.h>
#include <gsParasolid/gsWriteParasolid.hpp>

namespace gismo {
namespace extensions {

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsGeometry<real_t> & ggeo, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsMultiPatch<real_t> & mp, std::string const & filename );


}// extensions

}// gismo
