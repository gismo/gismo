#include <gsCore/gsLinearAlgebra.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsMatrix/cinterface/gsCMatrixInt.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create(void)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCMatrixInt*>(new gsMatrix<int>());     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create_rc(int rows, int cols)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCMatrixInt*>(new gsMatrix<int>(rows,cols));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create_rcd(int rows, int cols, int * data)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCMatrixInt*>(new gsMatrix<int>(gsAsMatrix<int>(data,rows,cols)));     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMatrixInt_delete(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN delete RICAST_Mi(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsMatrixInt_print(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_Mi(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int * gsMatrixInt_data(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Mi(m)->data();     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMatrixInt_transposeInPlace(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Mi(m)->transposeInPlace();     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsMatrixInt_rows(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Mi(m)->rows();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsMatrixInt_cols(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN return RICAST_Mi(m)->cols();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT void gsMatrixInt_setZero(gsCMatrixInt * m)
{
    GISMO_CAPI_BEGIN RICAST_Mi(m)->setZero();     GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif