#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlUtils.hpp>

namespace gismo
{

// Explicit instantiation

namespace internal
{
    // CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<index_t> >;
    // CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<bool> >;

    CLASS_TEMPLATE_INST gsXml< gsFunction<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsGeometry<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsCurve<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSurface<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsBasis<real_t> >;




    //CLASS_TEMPLATE_INST gsXml< gsBezier<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsPde<real_t>        >;
//    CLASS_TEMPLATE_INST gsXml< gsSurfacePoissonPde<real_t> >;



} // end namespace internal


} // end namespace gismo
