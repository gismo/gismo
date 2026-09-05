#include <gsCore/gsLinearAlgebra.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCVectorInt.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCVectorInt * gsVectorInt_create(void)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVectorInt*>(new gsVector<int>());     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCVectorInt * gsVectorInt_create_r(int rows)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVectorInt*>(new gsVector<int>(rows));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCVectorInt * gsVectorInt_create_rd(int rows, int * data)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCVectorInt*>(new gsVector<int>(gsAsVector<int>(data,rows)));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsVectorInt_delete(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN delete RICAST_Vi(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsVectorInt_print(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_Vi(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int * gsVectorInt_data(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Vi(m)->data();     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsVectorInt_transposeInPlace(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Vi(m)->transposeInPlace();     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsVectorInt_rows(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Vi(m)->rows();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsVectorInt_cols(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Vi(m)->cols();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT void gsVectorInt_setZero(gsCVectorInt * m)
{
    GISMO_CAPI_BEGIN RICAST_Vi(m)->setZero();     GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif