#include <gsCore/gsLinearAlgebra.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsMatrix/cinterface/gsCVector.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCVector * gsVector_create(void)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVector*>(new gsVector<double>());     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCVector * gsVector_create_r(int rows)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVector*>(new gsVector<double>(rows));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCVector * gsVector_create_rd(int rows, double * data)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVector*>(new gsVector<double>(gsAsVector<double>(data,rows)));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsVector_delete(gsCVector * m)
{
    GISMO_CAPI_BEGIN delete RICAST_V(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsVector_print(gsCVector * m)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_V(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT double * gsVector_data(gsCVector * m)
{
    GISMO_CAPI_BEGIN return RICAST_V(m)->data();     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsVector_transposeInPlace(gsCVector * m)
{
    GISMO_CAPI_BEGIN return RICAST_V(m)->transposeInPlace();     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsVector_rows(gsCVector * m)
{
    GISMO_CAPI_BEGIN return RICAST_V(m)->rows();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsVector_cols(gsCVector * m)
{
    GISMO_CAPI_BEGIN return RICAST_V(m)->cols();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT void gsVector_setZero(gsCVector * m)
{
    GISMO_CAPI_BEGIN RICAST_V(m)->setZero();     GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif