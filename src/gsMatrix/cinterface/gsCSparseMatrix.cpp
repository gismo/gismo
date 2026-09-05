#include <gsCore/gsLinearAlgebra.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsMacros.h>
#include <gsMatrix/cinterface/gsCSparseMatrix.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCSparseMatrix * gsSparseMatrix_create(void)
{
    GISMO_CAPI_BEGIN return RICAST_CSM(new gsSparseMatrix<double>());     GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsSparseMatrix_delete(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN delete RICAST_SM(m);     GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsSparseMatrix_print(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN gsInfo<<*RICAST_SM(m);     GISMO_CAPI_END_VOID
}

// GISMO_EXPORT double * gsSparseMatrix_data(gsCSparseMatrix * m)
// { return RICAST_SM(m)->data(); }

GISMO_EXPORT double* gsSparseMatrix_valuePtr(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->valuePtr();     GISMO_CAPI_END(NULL)
} // get pointer to matrix values
GISMO_EXPORT int*    gsSparseMatrix_innerIndexPtr(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->innerIndexPtr();     GISMO_CAPI_END(NULL)
} // get pointer to matrix rows
GISMO_EXPORT int*    gsSparseMatrix_outerIndexPtr(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->outerIndexPtr();     GISMO_CAPI_END(NULL)
} // get pointer to matrix columns

GISMO_EXPORT int gsSparseMatrix_rows(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->rows();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsSparseMatrix_cols(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->cols();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsSparseMatrix_nnz(gsCSparseMatrix * m)
{
    GISMO_CAPI_BEGIN return RICAST_SM(m)->nonZeros();     GISMO_CAPI_END(-1)
}

GISMO_EXPORT void gsSparseMatrix_setFromTriplets(gsCSparseMatrix * m, gsCVectorInt * rows, gsCVectorInt * cols, gsCVector * values)
{
    GISMO_CAPI_BEGIN
    auto * R = RICAST_Vi(rows);
    auto * C = RICAST_Vi(cols);
    auto * V = RICAST_V(values);

    GISMO_ENSURE(R->size() == C->size() && R->size() == V->size(), "Input vectors must have the same size.");

    gsSparseEntries<double> entries;
    entries.reserve(R->size());
    for (int i = 0; i < R->size(); i++)
        entries.add((*R)[i], (*C)[i], (*V)[i]);

    RICAST_SM(m)->resize(R->size(), C->size());
    RICAST_SM(m)->setFrom(entries);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsSparseMatrix_intoTriplets(gsCSparseMatrix * m, gsCVectorInt * rows, gsCVectorInt * cols, gsCVector * vals)
{
    GISMO_CAPI_BEGIN
    auto * sm = RICAST_SM(m);
    auto * R = RICAST_Vi(rows);
    auto * C = RICAST_Vi(cols);
    auto * V = RICAST_V(vals);

    R->resize(sm->nonZeros());
    C->resize(sm->nonZeros());
    V->resize(sm->nonZeros());
    auto R_it = R->begin();
    auto C_it = C->begin();
    auto V_it = V->begin();

    for (int i = 0; i!=sm->outerSize(); i++)
    {
        for (typename gsSparseMatrix<double>::InnerIterator it(*sm,i); it;
                    ++it, ++R_it, ++C_it, ++V_it)
        {
            *R_it = it.row();
            *C_it = it.col();
            *V_it = it.value();
        }
    }
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif