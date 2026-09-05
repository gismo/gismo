/* gsMesh XML instantiations (moved from gsIO/gsXmlInstance_.cpp,
   modularization stream S3 step A6) */
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsUtils/gsMesh/gsMeshXml.hpp>

namespace gismo {
namespace internal {

CLASS_TEMPLATE_INST gsXml< gsMesh<real_t> >;

} // namespace internal
} // namespace gismo
