
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCVectorInt * gsVectorInt_create(void);
    GISMO_EXPORT gsCVectorInt * gsVectorInt_create_r (int rows);
    GISMO_EXPORT gsCVectorInt * gsVectorInt_create_rd(int rows, int * data);

    GISMO_EXPORT void gsVectorInt_delete(gsCVectorInt * m);
    GISMO_EXPORT void gsVectorInt_print(gsCVectorInt * m);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    GISMO_EXPORT int * gsVectorInt_data(gsCVectorInt * m);

    GISMO_EXPORT void gsVectorInt_transposeInPlace(gsCVectorInt * m);
    GISMO_EXPORT int gsVectorInt_rows(gsCVectorInt * m);
    GISMO_EXPORT int gsVectorInt_cols(gsCVectorInt * m);
    GISMO_EXPORT void gsVectorInt_setZero(gsCVectorInt * m);
    
#ifdef __cplusplus
}
#endif
