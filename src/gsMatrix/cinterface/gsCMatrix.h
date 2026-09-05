
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCMatrix * gsMatrix_create(void);                     // create matrix object
    GISMO_EXPORT gsCMatrix * gsMatrix_create_rc (int rows, int cols);   // create matrix object with rows+cols
    GISMO_EXPORT gsCMatrix * gsMatrix_create_rcd(int rows, int cols, double * data); // copy data into matrix

//    typedef struct NoArg NoArg;
//#define gsMatrix_create(...) gsMatrix_create_(__VA_ARGS__ __VA_OPT__(,) (NoArg*)0, ~ )(__VA_ARGS__)
//#define gsMatrix_create_(X, ...) _Generic((X),                     
//                                          NoArg*: gsMatrix_create_,  
//        int: _Generic((FIRST(__VA_ARGS__,)), 
//       int: foo_float_int))(_1, __VA_ARGS__)
//#define FIRST(A, ...) A
    GISMO_EXPORT void gsMatrix_delete(gsCMatrix * m);
    GISMO_EXPORT void gsMatrix_print(gsCMatrix * m);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    GISMO_EXPORT double* gsMatrix_data(gsCMatrix * m);   // get pointer to matrix data

    GISMO_EXPORT void gsMatrix_transposeInPlace(gsCMatrix * m);
    GISMO_EXPORT int gsMatrix_rows(gsCMatrix * m);
    GISMO_EXPORT int gsMatrix_cols(gsCMatrix * m);
    GISMO_EXPORT void gsMatrix_setZero(gsCMatrix * m);
    
#ifdef __cplusplus
}
#endif
