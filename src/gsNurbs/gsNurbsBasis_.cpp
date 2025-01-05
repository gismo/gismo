/* Symbol export for G+Smo shared object */

#define gsNurbsBasis_EXPORT

#include <gsIO/gsXml.h>
#include <gsCore/gsBasisFun.h>
#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsNurbsBasis.hpp>

namespace gismo
{

// CLASS_TEMPLATE_INST gsTensorBSplineBasis<1, real_t>;
// CLASS_TEMPLATE_INST gsNurbsBasis<real_t>;

CLASS_TEMPLATE_INST internal::gsXml< gsNurbsBasis<real_t> >;

/*
#ifdef GISMO_WITH_PYBIND11

//namespace py = pybind11;

#endif
*/

}
