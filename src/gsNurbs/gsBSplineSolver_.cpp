
#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsBSplineSolver.h>
#include <gsNurbs/gsBSplineSolver.hpp>

#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

  CLASS_TEMPLATE_INST gsBSplineSolver<real_t>;
  TEMPLATE_INST       unsigned findHyperPlaneIntersections<real_t> (
         const gsBSpline<real_t>    &curve, const gsVector<real_t>     &normal,
        real_t reference, real_t tolerance, std::vector<Root<real_t> > &roots  );

}
