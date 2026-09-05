
#include <gsModeling/gsFitting.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsModeling/cinterface/gsCFitting.h>
#include <gsCInterface/gsCError.h>
#include <limits>

#ifdef __cplusplus
extern "C"
{
#endif

using namespace gismo;

GISMO_EXPORT gsCFitting * gsFitting_create(gsCMatrix * param_values, gsCMatrix * points, gsCBasis * basis)
{
    GISMO_CAPI_BEGIN
    auto * param_values_ptr = RICAST_M(param_values);
    auto * points_ptr = RICAST_M(points);
    auto * basis_ptr = RICAST_B(basis);
    return reinterpret_cast<gsCFitting*>(new gsFitting<double>(*param_values_ptr, *points_ptr, *basis_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsFitting_delete(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    delete reinterpret_cast<gsFitting<double>*>(fitter);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFitting_compute(gsCFitting * fitter, double lambda)
{
    GISMO_CAPI_BEGIN
    reinterpret_cast<gsFitting<double>*>(fitter)->compute(lambda);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFitting_parameterCorrection(gsCFitting * fitter, double accuracy, int maxIter, double tolOrth)
{
    GISMO_CAPI_BEGIN
    reinterpret_cast<gsFitting<double>*>(fitter)->parameterCorrection(accuracy, maxIter, tolOrth);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsFitting_computeErrors(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    reinterpret_cast<gsFitting<double>*>(fitter)->computeErrors();
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT double gsFitting_minPointError(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsFitting<double>*>(fitter)->minPointError();
    GISMO_CAPI_END(std::numeric_limits<double>::quiet_NaN())
}

GISMO_EXPORT double gsFitting_maxPointError(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsFitting<double>*>(fitter)->maxPointError();
    GISMO_CAPI_END(std::numeric_limits<double>::quiet_NaN())
}

GISMO_EXPORT double* gsFitting_pointWiseErrors(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    const double * errors = reinterpret_cast< const double* >(reinterpret_cast<gsFitting<double>*>(fitter)->pointWiseErrors().data());
    return const_cast<double *>(errors);
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT int gsFitting_numPointsBelow(gsCFitting * fitter, double threshold)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsFitting<double>*>(fitter)->numPointsBelow(threshold);
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT gsCGeometry* gsFitting_result(gsCFitting * fitter)
{
    GISMO_CAPI_BEGIN
    auto * result = reinterpret_cast<gsFitting<double>*>(fitter)->result();
    return RICAST_CG(result->clone().release());
    GISMO_CAPI_END(NULL)
}

#ifdef __cplusplus
}
#endif