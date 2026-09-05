
#include <gsCore/gsGeometry.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsFuncData.h>
#include <gsIO/gsFileData.h>
#include <gsTensor/gsGridIterator.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsCore/cinterface/gsCGeometry.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>
#include <limits>

#ifdef __cplusplus
extern "C"
{
#endif

using namespace gismo;

GISMO_EXPORT gsCGeometry * gsGeometry_read(char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data(filename);
    if (data.hasAny< gsGeometry<> >())
    {
        gsGeometry<>::uPtr ptr = data.getAnyFirst< gsGeometry<> >();
        return RICAST_CG(ptr.release());
    }
    else
    {
        gsWarn<<"[G+Smo] No gsGeometry found in file "<<filename<<"\n";
        return NULL;
    }
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsGeometry_write(gsCGeometry * obj, char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data;
    data.add(*RICAST_G(obj));
    data.save(filename);
    GISMO_CAPI_END_VOID
}


GISMO_EXPORT gsCGeometry* gsGeometry_clone(gsCGeometry * g)
{
    GISMO_CAPI_BEGIN
    return RICAST_CG(RICAST_G(g)->clone().release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsGeometry_basis(gsCGeometry * g)
{
    GISMO_CAPI_BEGIN
    return RICAST_CB(RICAST_G(g)->basis().clone().release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsGeometry_coefs_into(gsCGeometry * g, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    *RICAST_M(coefs) = RICAST_G(g)->coefs();
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_setCoefs(gsCGeometry * g, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    RICAST_G(g)->coefs() = *RICAST_M(coefs);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_uniformRefine(gsCGeometry * g, int numKnots, int mul, int dir)
{
    GISMO_CAPI_BEGIN RICAST_G(g)->uniformRefine(numKnots, mul, dir);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_refineElements(gsCGeometry * g, int * boxData, int boxSize)
{
    GISMO_CAPI_BEGIN
    std::vector<int> boxes(boxData,boxData+boxSize);
    RICAST_G(g)->refineElements(boxes);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_refine(gsCGeometry * g, gsCMatrix * boxes, int refExt)
{
    GISMO_CAPI_BEGIN
    std::vector<int> boxData = RICAST_G(g)->basis().asElements(*RICAST_M(boxes),refExt);
    RICAST_G(g)->refineElements(boxData);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_degreeElevate(gsCGeometry * g, int i, int dir)
{
    GISMO_CAPI_BEGIN RICAST_G(g)->degreeElevate(i, dir);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_recoverPoints(gsCGeometry * g, gsCMatrix * uv, gsCMatrix * xyz,
                                            int k, double accuracy)
{
    GISMO_CAPI_BEGIN RICAST_G(g)->recoverPoints(*RICAST_M(xyz), *RICAST_M(uv), k, accuracy);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_recoverPointGrid(gsCGeometry * g, gsCVector * a, gsCVector * b,
                                               gsCVectorInt * sz,  gsCMatrix * xyz,
                                               gsCMatrix * uv,int c, double accuracy)
{
    GISMO_CAPI_BEGIN
    gismo::gsGridIterator<real_t,gismo::CUBE> git(*RICAST_V(a),*RICAST_V(b), *RICAST_Vi(sz));
    RICAST_G(g)->recoverPointGrid(git,*RICAST_M(xyz), *RICAST_M(uv), c, accuracy);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_invertPoints(gsCGeometry * g, gsCMatrix * points, gsCMatrix * result,
                                double accuracy)
{
    GISMO_CAPI_BEGIN RICAST_G(g)->invertPoints(*RICAST_M(points), *RICAST_M(result), accuracy,false);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_invertPointGrid(gsCGeometry * g, gsCVector * a, gsCVector * b,
                                   gsCVectorInt * sz, gsCMatrix * result,
                                   double accuracy)
{
    GISMO_CAPI_BEGIN
    gismo::gsGridIterator<real_t,gismo::CUBE> git(*RICAST_V(a),*RICAST_V(b), *RICAST_Vi(sz));
    RICAST_G(g)->invertPointGrid(git, *RICAST_M(result), accuracy,false);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsGeometry_normal_into(gsCGeometry * fs,
                             gsCMatrix * u,
                             gsCMatrix * result)
{
    GISMO_CAPI_BEGIN
    gismo::gsMapData<> mapData;
    mapData.flags = gismo::NEED_NORMAL ;
    mapData.points = *RICAST_M(u);
    RICAST_G(fs)->computeMap(mapData);
    *RICAST_M(result) = mapData.normals;
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT double gsGeometry_closestPointTo(gsCGeometry * fs,
                                 gsCMatrix * pt,
                                 gsCMatrix * result,
                                 double accuracy)
{
    GISMO_CAPI_BEGIN
    gismo::gsVector<> tmp;
    double dist = RICAST_G(fs)->closestPointTo(*RICAST_V(pt),tmp,accuracy,false);
    *RICAST_M(result) = tmp;
    // gsDebugVar(*RICAST_M(result));
    // gsDebugVar(*RICAST_V(result));
    return dist;
    GISMO_CAPI_END(std::numeric_limits<double>::quiet_NaN())
}

#ifdef __cplusplus
}
#endif
