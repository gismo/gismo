/* gsMatrix XML instantiations (moved from gsIO/gsXmlInstance_.cpp,
   modularization stream S3 step A6) */
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsXml.h>
#include <gsMatrix/gsMatrixXml.hpp>

namespace gismo {
namespace internal {

CLASS_TEMPLATE_INST gsXml< gsMatrix<real_t> >;
CLASS_TEMPLATE_INST gsXml< gsMatrix<index_t> >;
CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<real_t, RowMajor, index_t> >;
CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<real_t, ColMajor, index_t> >;

} // namespace internal
} // namespace gismo
