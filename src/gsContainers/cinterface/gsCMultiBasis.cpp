#include <gsContainers/gsMultiBasis.h>
#include <gsIO/gsFileData.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMatrix.h>
#include <gsContainers/cinterface/gsCMultiBasis.h>
#include <gsMatrix/cinterface/gsCMemory.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create()
{
    GISMO_CAPI_BEGIN
    return RICAST_CMP(new gsMultiBasis<double>());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create_basis(gsCBasis * basis)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = RICAST_B(basis);
    return RICAST_CMP(new gsMultiBasis<double>(*basis_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create_patches(gsCMultiPatch * mp)
{
    GISMO_CAPI_BEGIN
    auto * mp_ptr = RICAST_MP(mp);
    return RICAST_CMP(new gsMultiBasis<double>(*mp_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMultiBasis_addBasis(gsCMultiBasis* mb, gsCBasis* basis)
{
    GISMO_CAPI_BEGIN
    auto * mb_ptr = RICAST_MB(mb);
    auto * basis_ptr = RICAST_B(basis);
    mb_ptr->addBasis(basis_ptr);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT gsCBasis * gsMultiBasis_basis(gsCMultiBasis * mb, int i)
{
    GISMO_CAPI_BEGIN
    auto * mb_ptr = RICAST_MB(mb);
    return RICAST_CB(&mb_ptr->basis(i));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCMultiBasis * gsMultiBasis_read(char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data(filename);
    if (data.hasAny< gsMultiBasis<> >())
    {
        gsMultiBasis<>::uPtr ptr = data.getAnyFirst< gsMultiBasis<> >();
        return RICAST_CMB(ptr.release());
    }
    else
    {
        gsWarn<<"[G+Smo] No gsMultiBasis found in file "<<filename<<"\n";
        return NULL;
    }
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsMultiBasis_write(gsCMultiBasis * obj, char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data;
    data.add(*RICAST_MB(obj));
    data.save(filename);
    GISMO_CAPI_END_VOID
}


#ifdef __cplusplus
}
#endif
