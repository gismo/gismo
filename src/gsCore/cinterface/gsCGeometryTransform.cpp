
#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometryTransform.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsCore/cinterface/gsCGeometry.h>
#include <gsCore/cinterface/gsCGeometryTransform.h>
#include <gsDomain/cinterface/gsCKnotVector.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCGeometryTransform * gsGeometryTransform_create(gsCGeometry* g, gsCMatrix * m,
                                                                   gsCVector * v)
    {
        GISMO_CAPI_BEGIN
        auto * g_ptr = RICAST_G(g);
        auto * mm = RICAST_M(m);
        auto * vv = RICAST_V(v);
        return RICAST_CG(new gsGeometryTransform<double>(g_ptr,*mm, *vv));
            GISMO_CAPI_END(NULL)
    }


#ifdef __cplusplus
}
#endif
