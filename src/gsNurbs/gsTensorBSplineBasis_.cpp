
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>

#include <gsNurbs/gsCompactKnotVector.h>

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorBasis.hpp>


#include <gsIO/gsXmlUtils.hpp>

namespace gismo
{

//TEMPLATE_INST gsXmlNode * putTensorBasisToXml 


// instantiate tensor B-spline bases
CLASS_TEMPLATE_INST gsTensorBasis<1, gsBSplineBasis<real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<2, gsBSplineBasis<real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<3, gsBSplineBasis<real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<4, gsBSplineBasis<real_t, gsKnotVector<real_t> > >;

// instantiate tensor B-spline bases with compact knot vector
CLASS_TEMPLATE_INST gsTensorBasis<1, gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<2, gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<3, gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST gsTensorBasis<4, gsBSplineBasis<real_t, gsCompactKnotVector<real_t> > >; 

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t, gsKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t, gsKnotVector<real_t> > >;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<1,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<2,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<3,real_t, gsCompactKnotVector<real_t> > >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<4,real_t, gsCompactKnotVector<real_t> > >;

}
