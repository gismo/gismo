#include <gsCore/gsFunctionSet.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsCore/cinterface/gsCFunctionSet.h>
#include <gsMatrix/cinterface/gsCMemory.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT void gsFunctionSet_print(gsCFunctionSet * fs)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_F(fs)<<"\n";     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFunctionSet_delete(gsCFunctionSet * ptr)
{
    GISMO_CAPI_BEGIN delete RICAST_F(ptr);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsFunctionSet_domainDim(gsCFunctionSet * fs)
{
    GISMO_CAPI_BEGIN return RICAST_F(fs)->domainDim();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsFunctionSet_targetDim(gsCFunctionSet * fs)
{
    GISMO_CAPI_BEGIN return RICAST_F(fs)->targetDim();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsFunctionSet_nPieces(gsCFunctionSet * fs)
{
    GISMO_CAPI_BEGIN return RICAST_F(fs)->nPieces();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT gsCMatrix* gsFunctionSet_support(gsCFunctionSet * fs)
{
    GISMO_CAPI_BEGIN return RICAST_CM( new gsMatrix<double>(RICAST_F(fs)->support()) );     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsFunctionSet_eval_into(gsCFunctionSet * fs,
                            gsCMatrix * u,
                            gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_F(fs)->eval_into(*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFunctionSet_deriv_into(gsCFunctionSet * fs,
                             gsCMatrix * u,
                             gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_F(fs)->deriv_into(*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFunctionSet_deriv2_into(gsCFunctionSet * fs,
                                            gsCMatrix * u,
                                            gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_F(fs)->deriv2_into(*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFunctionSet_evalAllDers1_into(gsCFunctionSet * fs,
                                                  gsCMatrix * u,
                                                  gsCMatrix * values,
                                                  gsCMatrix * deriv)
{
    GISMO_CAPI_BEGIN
    std::vector<gsMatrix<> > result(2);
    RICAST_F(fs)->evalAllDers_into(*RICAST_M(u), 1, result, false);
    *RICAST_M(values) = give(result[0]);
    *RICAST_M(deriv) = give(result[1]);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFunctionSet_evalAllDers2_into(gsCFunctionSet * fs,
                                                  gsCMatrix * u,
                                                  gsCMatrix * values,
                                                  gsCMatrix * deriv,
                                                  gsCMatrix * deriv2)
{
    GISMO_CAPI_BEGIN
    std::vector<gsMatrix<> > result(3);
    RICAST_F(fs)->evalAllDers_into(*RICAST_M(u), 2, result, false);
    *RICAST_M(values) = give(result[0]);
    *RICAST_M(deriv) = give(result[1]);
    *RICAST_M(deriv2) = give(result[2]);
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif
