#include <gsContainers/gsMultiPatch.h>
#include <gsIO/gsFileData.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsContainers/cinterface/gsCMultiPatch.h>
#include <gsMatrix/cinterface/gsCMemory.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCMultiPatch* gsMultiPatch_create()
{
    GISMO_CAPI_BEGIN
    return RICAST_CMP(new gsMultiPatch<double>());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiPatch* gsMultiPatch_create_geometry(gsCGeometry * geo)
{
    GISMO_CAPI_BEGIN
    auto * geo_ptr = RICAST_G(geo);
    return RICAST_CMP(new gsMultiPatch<double>(*geo_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMultiPatch_addPatch(gsCMultiPatch* mp, gsCGeometry* geo)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    auto * geo_ptr = RICAST_G(geo);
    mp_ptr->addPatch(*geo_ptr);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT gsCBasis * gsMultiPatch_basis(gsCMultiPatch * mp, int i)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    return RICAST_CB(&mp_ptr->basis(i));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry * gsMultiPatch_patch(gsCMultiPatch * mp, int i)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    return RICAST_CG(&mp_ptr->patch(i));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMultiPatch_computeTopology(gsCMultiPatch * mp)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    mp_ptr->computeTopology();
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsMultiPatch_embed(gsCMultiPatch * mp, int dim)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    mp_ptr->embed(dim);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsMultiPatch_uniformRefine(gsCMultiPatch * mp, int numKnots, int mul, int dir)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    mp_ptr->uniformRefine(numKnots, mul, dir);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsMultiPatch_degreeElevate(gsCMultiPatch * mp, int i, int dir)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    mp_ptr->degreeElevate(i, dir);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT gsCMultiPatch * gsMultiPatch_read(char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data(filename);
    if (data.hasAny< gsMultiPatch<> >())
    {
        gsMultiPatch<>::uPtr ptr = data.getAnyFirst< gsMultiPatch<> >();
        return RICAST_CMP(ptr.release());
    }
    else
    {
        gsWarn<<"[G+Smo] No gsMultiPatch found in file "<<filename<<"\n";
        return NULL;
    }
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMultiPatch_write(gsCMultiPatch * obj, char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data;
    data.add(*RICAST_MP(obj));
    data.save(filename);
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif
