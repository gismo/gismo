#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlUtils.hpp>

namespace gismo
{

// Explicit instantiation

namespace internal
{
    CLASS_TEMPLATE_INST gsXml< gsMatrix<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSparseMatrix<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsGeometry<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsCurve<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSurface<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsBasis<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMultiPatch<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsSolid<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTrimSurface<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsBSpline<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsBSplineBasis<real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsNurbs<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsNurbsBasis<real_t> >;
        
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbs<4,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTensorNurbsBasis<4,real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHBSplineBasis<4,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHBSpline<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHBSpline<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHBSpline<4,real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<1,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSplineBasis<4,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSpline<1,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSpline<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSpline<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsTHBSpline<4,real_t> >;

    CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<2,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<3,real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsHTensorBasis<4,real_t> >;


    //CLASS_TEMPLATE_INST gsXml< gsBezier<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsMesh<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsCurveFitting<real_t> >;
    
    CLASS_TEMPLATE_INST gsXml< gsPde<real_t>        >;
    CLASS_TEMPLATE_INST gsXml< gsPoissonPde<real_t> >;
//    CLASS_TEMPLATE_INST gsXml< gsSurfacePoissonPde<real_t> >;
    CLASS_TEMPLATE_INST gsXml< gsBVProblem<real_t>  >;

    CLASS_TEMPLATE_INST gsXml< gsPlanarDomain<real_t> >;
    
    TEMPLATE_INST
    gsXmlNode * makeNode<real_t>( const std::string & name, 
                             const gsMatrix<real_t> & value, 
                             gsXmlTree & data,
                             bool transposed );
    TEMPLATE_INST
    void getMatrixFromXml<real_t> ( gsXmlNode * node, 
                               unsigned const & rows, 
                               unsigned const & cols, gsMatrix<real_t> & result );

    TEMPLATE_INST
    void getSparseEntriesFromXml ( gsXmlNode * node, 
                                   gsSparseEntries<real_t> & result );

    TEMPLATE_INST
    void getFunctionFromXml<real_t> ( gsXmlNode * node, 
                                      gsMFunctionExpr<real_t> & result );
    
    TEMPLATE_INST
    gsXmlNode * putMatrixToXml<real_t> ( gsMatrix<real_t> const & mat, 
                                         gsXmlTree & data, 
                                         std::string name );

    TEMPLATE_INST
    gsXmlNode * putSparseMatrixToXml ( gsSparseMatrix<real_t> const & mat, 
                                       gsXmlTree & data, std::string name);

} // end namespace internal


} // end namespace gismo
