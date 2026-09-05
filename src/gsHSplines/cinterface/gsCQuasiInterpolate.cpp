
// Header-only builds: the B-spline basis header must come first, it
// completes gsTensorBSplineBasis before gsSurfMesh.hpp needs it
#include <gsNurbs/gsBSplineBasis.h>
#include <gsHSplines/gsQuasiInterpolate.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsAssembler/cinterface/gsCQuadRule.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif


GISMO_EXPORT void gsQuasiInterpolate_localIntpl_into( gsCBasis * basis,
                                                      gsCFunctionSet * fun,
                                                      gsCMatrix * result)
{
    GISMO_CAPI_BEGIN
    auto basis_ptr = RICAST_B(basis);
    auto fun_ptr = RICAST_F(fun);
    auto result_ptr = RICAST_M(result);
    gsQuasiInterpolate<double>::localIntpl(*basis_ptr, fun_ptr->function(0), *result_ptr);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsQuasiInterpolate_Taylor_into( gsCBasis * basis,
                                                  gsCFunctionSet * fun,
                                                  int deg,
                                                  gsCMatrix * result)
{
    GISMO_CAPI_BEGIN
    auto basis_ptr = RICAST_B(basis);
    auto fun_ptr = RICAST_F(fun);
    auto result_ptr = RICAST_M(result);
    gsQuasiInterpolate<double>::Taylor(*basis_ptr, fun_ptr->function(0), deg, *result_ptr);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsQuasiInterpolate_Schoenberg_into( gsCBasis * basis,
                                                      gsCFunctionSet * fun,
                                                      gsCMatrix * result)
{
    GISMO_CAPI_BEGIN
    auto basis_ptr = RICAST_B(basis);
    auto fun_ptr = RICAST_F(fun);
    auto result_ptr = RICAST_M(result);
    gsQuasiInterpolate<double>::Schoenberg(*basis_ptr, fun_ptr->function(0), *result_ptr);
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif