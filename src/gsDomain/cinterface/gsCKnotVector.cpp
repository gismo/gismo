
#include <gsDomain/gsKnotVector.h>
#include <gsCInterface/gsCTypes.h>
#include <gsDomain/cinterface/gsCKnotVector.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCKnotVector * gsKnotVector_create(double* knots, int size)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsCKnotVector*>(new gsKnotVector<double>(-1,knots,knots+size));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsKnotVector_delete(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    delete reinterpret_cast<gsKnotVector<double>*>(kv);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsKnotVector_print(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    reinterpret_cast<gsKnotVector<double>*>(kv)->print(gsInfo);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT double* gsKnotVector_data(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    return const_cast<double *>(reinterpret_cast<gsKnotVector<double>*>(kv)->get().data());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT int gsKnotVector_size(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsKnotVector<double>*>(kv)->size();
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsKnotVector_uSize(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsKnotVector<double>*>(kv)->uSize();
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsKnotVector_numElements(gsCKnotVector * kv)
{
    GISMO_CAPI_BEGIN
    return reinterpret_cast<gsKnotVector<double>*>(kv)->numElements();
    GISMO_CAPI_END(-1)
}

// GISMO_EXPORT double* gsKnotVector_unique(gsCKnotVector * kv)
// {
//     return reinterpret_cast<gsKnotVector<double>*>(kv)->unique().data();
// }

//GISMO_EXPORT double * greville(void * kv)
//{
//    gsKnotVector<double>* _kv = static_cast<gsKnotVector<double>*>(kv);
//}

#ifdef __cplusplus
}
#endif