
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCVector * gsVector_create(void);
    GISMO_EXPORT gsCVector * gsVector_create_r (int rows);
    GISMO_EXPORT gsCVector * gsVector_create_rd(int rows, double * data);
    GISMO_EXPORT void gsVector_delete(gsCVector * m);
    GISMO_EXPORT void gsVector_print(gsCVector * m);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    GISMO_EXPORT double * gsVector_data(gsCVector * m);

    GISMO_EXPORT void gsVector_transposeInPlace(gsCVector * m);
    GISMO_EXPORT int gsVector_rows(gsCVector * m);
    GISMO_EXPORT int gsVector_cols(gsCVector * m);
    GISMO_EXPORT void gsVector_setZero(gsCVector * m);

#ifdef __cplusplus
}
#endif
