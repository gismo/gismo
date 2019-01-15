
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSpline.hpp>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{
TEMPLATE_INST
void constructCoefsForSlice<1, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 1> & sizes,
                                       gsMatrix<real_t> & result
                                      );

TEMPLATE_INST
void constructCoefsForSlice<2, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 2> & sizes,
                                       gsMatrix<real_t> & result
                                      );

TEMPLATE_INST
void constructCoefsForSlice<3, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 3> & sizes,
                                       gsMatrix<real_t> & result
                                      );
TEMPLATE_INST
void constructCoefsForSlice<4, real_t>(index_t dir_fixed,
                                       const index_t index,
                                       const gsMatrix<real_t> & fullCoefs,
                                       const gsVector<index_t, 4> & sizes,
                                       gsMatrix<real_t> & result
                                      );


CLASS_TEMPLATE_INST gsTensorBSpline<1,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<2,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<3,real_t>;
CLASS_TEMPLATE_INST gsTensorBSpline<4,real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<1,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<2,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<3,real_t> >;
CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSpline<4,real_t> >;

}
