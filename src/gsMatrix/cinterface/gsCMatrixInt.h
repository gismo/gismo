
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create(void);
    GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create_rc (int rows, int cols);
    GISMO_EXPORT gsCMatrixInt * gsMatrixInt_create_rcd(int rows, int cols, int * data);

//    typedef struct NoArg NoArg;
//#define gsMatrixInt_create(...) gsMatrixInt_create_(__VA_ARGS__ __VA_OPT__(,) (NoArg*)0, ~ )(__VA_ARGS__)
//#define gsMatrixInt_create_(X, ...) _Generic((X),
//                                          NoArg*: gsMatrixInt_create_,
//        int: _Generic((FIRST(__VA_ARGS__,)), 
//       int: foo_float_int))(_1, __VA_ARGS__)
//#define FIRST(A, ...) A
    GISMO_EXPORT void gsMatrixInt_delete(gsCMatrixInt * m);
    GISMO_EXPORT void gsMatrixInt_print(gsCMatrixInt * m);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    GISMO_EXPORT int* gsMatrixInt_data(gsCMatrixInt * m);

    GISMO_EXPORT void gsMatrixInt_transposeInPlace(gsCMatrixInt * m);
    GISMO_EXPORT int gsMatrixInt_rows(gsCMatrixInt * m);
    GISMO_EXPORT int gsMatrixInt_cols(gsCMatrixInt * m);
    GISMO_EXPORT void gsMatrixInt_setZero(gsCMatrixInt * m);
    
#ifdef __cplusplus
}
#endif
