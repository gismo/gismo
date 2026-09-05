
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsNurbsCreator.h>
#include <gsContainers/gsMultiPatch.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsCore/cinterface/gsCGeometry.h>
#include <gsDomain/cinterface/gsCKnotVector.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineUnitInterval(int deg)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineUnitInterval(deg).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineRectangle(double low_x, double low_y, double upp_x, double upp_y, double turndeg)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineRectangle(low_x, low_y, upp_x, upp_y, turndeg).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineTrapezium(double Lbot, double Ltop, double H, double d, double turndeg)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineTrapezium(Lbot, Ltop, H, d, turndeg).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineSquare(double r, double x, double y)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineSquare(r, x, y).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineSquareGrid(int n, int m, double r, double lx, double ly)
{
    GISMO_CAPI_BEGIN
    return RICAST_CMP(new gsMultiPatch<double>(gsNurbsCreator<double>::BSplineSquareGrid(n, m, r, lx, ly)));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineCube(double r, double x, double y, double z)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineCube(r, x, y, z).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineCubeGrid(int n, int m, int p, double r, double lx, double ly, double lz)
{
    GISMO_CAPI_BEGIN
    return RICAST_CMP(new gsMultiPatch<double>(gsNurbsCreator<double>::BSplineCubeGrid(n, m, p, r, lx, ly, lz)));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsQuarterAnnulus(double r0, double r1)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::NurbsQuarterAnnulus(r0, r1).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsAnnulus(double r0, double r1)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::NurbsAnnulus(r0, r1).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineSaddle()
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineSaddle().release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsSphere(double r, double x, double y, double z)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::NurbsSphere(r, x, y, z).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsCircle(double r, double x, double y)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::NurbsCircle(r, x, y).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineTriangle(double H, double W)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(gsNurbsCreator<double>::BSplineTriangle(H, W).release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineStar(int N, double R0, double R1)
{
    GISMO_CAPI_BEGIN
    return RICAST_CMP(new gsMultiPatch<double>(gsNurbsCreator<double>::BSplineStar(N, R0, R1)));
    GISMO_CAPI_END(NULL)
}

#ifdef __cplusplus
}
#endif
