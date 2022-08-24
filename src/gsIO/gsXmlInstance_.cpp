#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlUtils.hpp>

namespace gismo
{

// Explicit instantiation

namespace internal
{
    CLASS_TEMPLATE_INST gsXml< gsMatrix<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMatrix<index_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<real_t> >;
    // CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<index_t> >;
    // CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<bool> >;

    CLASS_TEMPLATE_INST gsXml< gsGeometry<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsCurve<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSurface<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsBasis<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMultiPatch<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMultiBasis<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSolid<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTrimSurface<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsNurbs<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsNurbsBasis<real_t> >;
        
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<4,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<4,real_t> >;

    //CLASS_TEMPLATE_INST gsXml< gsBezier<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMesh<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsCurveFitting<real_t> >;
    
    CLASS_TEMPLATE_INST gsXml< gsPde<real_t>        >;
    CLASS_TEMPLATE_INST gsXml< gsPoissonPde<real_t> >;
//    CLASS_TEMPLATE_INST gsXml< gsSurfacePoissonPde<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsCurveLoop<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsPlanarDomain<real_t> >;

    
} // end namespace internal


} // end namespace gismo
