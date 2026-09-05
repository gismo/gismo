/* Symbol export for G+Smo shared object */

#define gsNurbsBasis_EXPORT
#include <gsNurbs/gsNurbsBasis.h>

namespace gismo
{

CLASS_TEMPLATE_INST internal::gsXml< gsNurbsBasis<real_t> >;

/*
#ifdef GISMO_WITH_PYBIND11

//namespace py = pybind11;

#endif
*/

}
