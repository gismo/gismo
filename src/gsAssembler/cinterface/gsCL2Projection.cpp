
#include <gsAssembler/gsProjection.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsAssembler/cinterface/gsCQuadRule.h>
#include <gsCInterface/gsCError.h>
#include <limits>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT double gsL2Projection_into( gsCFunctionSet * projectionBasis,
                                         gsCMultiBasis * integrationBasis,
                                         gsCMultiPatch * geometryMap,
                                         gsCFunctionSet * sourceFunction,
                                         gsCMatrix * coefs,
                                         gsCOptionList * options)
{
    GISMO_CAPI_BEGIN
    auto * projBasis_ptr = RICAST_F(projectionBasis);
    auto * intBasis_ptr = RICAST_MB(integrationBasis);
    auto * geomMap_ptr = RICAST_MP(geometryMap);
    auto * sourceFunc_ptr = RICAST_F(sourceFunction);
    auto * coefs_ptr = RICAST_M(coefs);
    auto * options_ptr = reinterpret_cast<gsOptionList*>(options);
    double error= gsL2Projection<double>::project(*projBasis_ptr, *intBasis_ptr, *geomMap_ptr, *sourceFunc_ptr, *coefs_ptr, gsBoundaryConditions<double>(), *options_ptr);
    return error;
    GISMO_CAPI_END(std::numeric_limits<double>::quiet_NaN())
}

#ifdef __cplusplus
}
#endif