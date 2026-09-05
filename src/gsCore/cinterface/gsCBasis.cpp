
#include <gsCore/gsBasis.h>
#include <gsDomain/gsDomain.h>
#include <gsDomain/gsDomainIterator.h>
#include <gsIO/gsFileData.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;



#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCBasis * gsBasis_read(char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data(filename);
    if (data.hasAny< gsBasis<> >())
    {
        gsBasis<>::uPtr ptr = data.getAnyFirst< gsBasis<> >();
        return RICAST_CB(ptr.release());
    }
    else
    {
        gsWarn<<"[G+Smo] No gsBasis found in file "<<filename<<"\n";
        return NULL;
    }
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsBasis_write(gsCBasis * obj, char* filename)
{
    GISMO_CAPI_BEGIN
    gsFileData<> data;
    data.add(*RICAST_B(obj));
    data.save(filename);
    GISMO_CAPI_END_VOID
}


//
// Methods, gsBasis
//

GISMO_EXPORT gsCBasis* gsBasis_clone(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    return RICAST_CB(RICAST_B(b)->clone().release());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsBasis_active_into(gsCBasis * b,
                              gsCMatrix * u,
                              gsCMatrixInt * result)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->active_into(*RICAST_M(u), *RICAST_Mi(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_evalSingle_into(gsCBasis * b,
                                          int i,
                                          gsCMatrix * u,
                                          gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->evalSingle_into(i,*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_derivSingle_into(gsCBasis * b,
                                           int i,
                                           gsCMatrix * u,
                                           gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->derivSingle_into(i,*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_deriv2Single_into(gsCBasis * b,
                                            int i,
                                            gsCMatrix * u,
                                            gsCMatrix * result)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->deriv2Single_into(i,*RICAST_M(u), *RICAST_M(result) );     GISMO_CAPI_END_VOID
}

GISMO_EXPORT gsCBasis * gsBasis_component(gsCBasis * b, int dir)
{
    GISMO_CAPI_BEGIN
    gsBasis<double> * c = & RICAST_B(b)->component(dir);
    return reinterpret_cast<gsCBasis*>(c);
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT int gsBasis_degree(gsCBasis * b, int dir)
{
    GISMO_CAPI_BEGIN return RICAST_B(b)->component(dir).degree(dir);     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsBasis_numElements(gsCBasis * b)
{
    GISMO_CAPI_BEGIN return RICAST_B(b)->numElements();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsBasis_dim(gsCBasis * b)
{
    GISMO_CAPI_BEGIN return RICAST_B(b)->dim();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsBasis_size(gsCBasis * b)
{
    GISMO_CAPI_BEGIN return RICAST_B(b)->size();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT gsCMatrix* gsBasis_support(gsCBasis * b, int i)
{
    GISMO_CAPI_BEGIN return reinterpret_cast<gsCMatrix*>( new gsMatrix<double>(RICAST_B(b)->support(i)) );     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsBasis_uniformRefine(gsCBasis * b, int numKnots, int mul, int dir)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->uniformRefine(numKnots, mul, dir);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_refineElements(gsCBasis * b, int * boxData, int boxSize)
{
    GISMO_CAPI_BEGIN
    std::vector<int> boxes(boxData,boxData+boxSize);
    RICAST_B(b)->refineElements(boxes);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_refine(gsCBasis * b, gsCMatrix * boxes, int refExt)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->refine(*RICAST_M(boxes),refExt);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_degreeElevate(gsCBasis * b, int i, int dir)
{
    GISMO_CAPI_BEGIN RICAST_B(b)->degreeElevate(i,dir);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_boundary_into(gsCBasis * b, int side, gsCMatrixInt * result)
{
    GISMO_CAPI_BEGIN *RICAST_Mi(result) = RICAST_B(b)->boundary(side);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_boundaryOffset_into(gsCBasis * b, int side, int offset, gsCMatrixInt * result)
{
    GISMO_CAPI_BEGIN *RICAST_Mi(result) = RICAST_B(b)->boundaryOffset(side,offset);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_elements_into(gsCBasis * b, gsCMatrix* elements)
{
    GISMO_CAPI_BEGIN
    auto * el = RICAST_M(elements);
    el->resize(RICAST_B(b)->domainDim(),2*RICAST_B(b)->numElements());
    auto domain = RICAST_B(b)->domain();
    auto domIt  = domain->beginAll();
    auto domEnd = domain->endAll();
    int id=0;
    for (; domIt<domEnd; ++domIt, ++id)
    {
        el->col(2*id) = domIt.lowerCorner();
        el->col(2*id+1) = domIt.upperCorner();
    }
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBasis_elementsBdr_into(gsCBasis * b, int side, gsCMatrix* elements)
{
    GISMO_CAPI_BEGIN
    auto * el = RICAST_M(elements);
    el->resize(RICAST_B(b)->domainDim(),2*RICAST_B(b)->numElements());
    auto domain = RICAST_B(b)->domain();
    auto domIt  = domain->beginBdr(side);
    auto domEnd = domain->endBdr(side);
    int id=0;
    for (; domIt<domEnd; ++domIt, ++id)
    {
        el->col(2*id) = domIt.lowerCorner();
        el->col(2*id+1) = domIt.upperCorner();
    }
    GISMO_CAPI_END_VOID
}

//
// Methods, Other
//



#ifdef __cplusplus
}
#endif