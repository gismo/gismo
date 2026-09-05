#include <gsCore/gsLinearAlgebra.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCMatrix * gsMatrix_create(void)
{
    GISMO_CAPI_BEGIN return RICAST_CM(new gsMatrix<double>());     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMatrix * gsMatrix_create_rc(int rows, int cols)
{
    GISMO_CAPI_BEGIN return RICAST_CM(new gsMatrix<double>(rows,cols));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMatrix * gsMatrix_create_rcd(int rows, int cols, double * data)
{
    GISMO_CAPI_BEGIN return RICAST_CM(new gsMatrix<double>(gsAsMatrix<double>(data,rows,cols)));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMatrix_delete(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN delete RICAST_M(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsMatrix_print(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_M(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT double * gsMatrix_data(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_M(m)->data();     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMatrix_transposeInPlace(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_M(m)->transposeInPlace();     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsMatrix_rows(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_M(m)->rows();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsMatrix_cols(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_M(m)->cols();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT void gsMatrix_setZero(gsCMatrix * m)
{
    GISMO_CAPI_BEGIN RICAST_M(m)->setZero();     GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif