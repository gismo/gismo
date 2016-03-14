/* Symbol export for G+Smo shared object */

#define gsKnotVector_EXPORT

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp> //dependancy

namespace gismo
{

//CLASS_TEMPLATE_INST gsKnotVector<real_t>;
CLASS_TEMPLATE_INST internal::gsXml< gsKnotVector<real_t> >;

}
