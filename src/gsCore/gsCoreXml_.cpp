/* gsCore container XML instantiations (moved from gsIO/gsXmlInstance_.cpp,
   modularization stream S3 step A6) */
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsCore/gsCoreXml.hpp>

namespace gismo {
namespace internal {

CLASS_TEMPLATE_INST gsXml< gsMultiPatch<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsMultiBasis<real_t> >;

} // namespace internal
} // namespace gismo
