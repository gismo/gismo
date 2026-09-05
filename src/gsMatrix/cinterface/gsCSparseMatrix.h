
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCSparseMatrix * gsSparseMatrix_create(void);                     // create matrix object

    GISMO_EXPORT void gsSparseMatrix_delete(gsCSparseMatrix * m);
    GISMO_EXPORT void gsSparseMatrix_print(gsCSparseMatrix * m);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    // GISMO_EXPORT double* gsSparseMatrix_data(gsCSparseMatrix * m);   // get pointer to matrix data

    GISMO_EXPORT double* gsSparseMatrix_valuePtr(gsCSparseMatrix * m); // get pointer to matrix values
    GISMO_EXPORT int*    gsSparseMatrix_innerIndexPtr(gsCSparseMatrix * m);    // get pointer to matrix rows
    GISMO_EXPORT int*    gsSparseMatrix_outerIndexPtr(gsCSparseMatrix * m);    // get pointer to matrix cols
    GISMO_EXPORT int     gsSparseMatrix_rows(gsCSparseMatrix * m);
    GISMO_EXPORT int     gsSparseMatrix_cols(gsCSparseMatrix * m);
    GISMO_EXPORT int     gsSparseMatrix_nnz(gsCSparseMatrix * m);


    GISMO_EXPORT void gsSparseMatrix_setFromTriplets(gsCSparseMatrix * m, gsCVectorInt * rows, gsCVectorInt * cols, gsCVector * values);
    GISMO_EXPORT void gsSparseMatrix_intoTriplets(gsCSparseMatrix * m, gsCVectorInt * rows, gsCVectorInt * cols, gsCVector * vals);


#ifdef __cplusplus
}
#endif
