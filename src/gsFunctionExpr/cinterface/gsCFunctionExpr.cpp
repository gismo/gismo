#include <gsFunctionExpr/gsFunctionExpr.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMemory.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCFunctionExpr * gsFunctionExpr1_create(const char * expression_string,
                                                      short_t ddim)
{
    GISMO_CAPI_BEGIN
    return RICAST_CF(new gsFunctionExpr<double>(expression_string,
                                                       ddim));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCFunctionExpr * gsFunctionExpr2_create(const char * expression_string1,
                                                      const char * expression_string2,
                                                      short_t ddim)
{
    GISMO_CAPI_BEGIN
    return RICAST_CF(new gsFunctionExpr<double>(expression_string1,
                                                       expression_string2, ddim));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCFunctionExpr * gsFunctionExpr3_create(const char * expression_string1,
                                                      const char * expression_string2,
                                                      const char * expression_string3,
                                                      short_t ddim)
{
    GISMO_CAPI_BEGIN
    return RICAST_CF(new gsFunctionExpr<double>(expression_string1,
                                                       expression_string2,
                                                       expression_string3,
                                                       ddim));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCFunctionExpr * gsFunctionExpr4_create(const char * expression_string1,
                                                      const char * expression_string2,
                                                      const char * expression_string3,
                                                      const char * expression_string4,
                                                      short_t ddim)
{
    GISMO_CAPI_BEGIN
    return RICAST_CF(new gsFunctionExpr<double>(expression_string1,
                                                       expression_string2,
                                                       expression_string3,
                                                       expression_string4,
                                                       ddim));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCFunctionExpr * gsFunctionExpr9_create(const char * expression_string1,
                                                      const char * expression_string2,
                                                      const char * expression_string3,
                                                      const char * expression_string4,
                                                      const char * expression_string5,
                                                      const char * expression_string6,
                                                      const char * expression_string7,
                                                      const char * expression_string8,
                                                      const char * expression_string9,
                                                      short_t ddim)
{
    GISMO_CAPI_BEGIN
    return RICAST_CF(new gsFunctionExpr<double>(expression_string1,
                                                       expression_string2,
                                                       expression_string3,
                                                       expression_string4,
                                                       expression_string5,
                                                       expression_string6,
                                                       expression_string7,
                                                       expression_string8,
                                                       expression_string9,
                                                       ddim));
    GISMO_CAPI_END(NULL)
}


#ifdef __cplusplus
}
#endif
