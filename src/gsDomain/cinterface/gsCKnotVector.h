
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCKnotVector * gsKnotVector_create(double* knots, int size);
    GISMO_EXPORT void gsKnotVector_delete(gsCKnotVector * kv);
    GISMO_EXPORT void gsKnotVector_print(gsCKnotVector * kv);
    /* LIFETIME: the returned pointer aliases internal storage of the
       object; it is invalidated by any modifying call on the same object
       and by its *_delete. Do not free it. */
    GISMO_EXPORT double* gsKnotVector_data(gsCKnotVector * kv);
    GISMO_EXPORT int gsKnotVector_size(gsCKnotVector * kv);
    GISMO_EXPORT int gsKnotVector_uSize(gsCKnotVector * kv);
    GISMO_EXPORT int gsKnotVector_numElements(gsCKnotVector * kv);
    // GISMO_EXPORT double* gsKnotVector_unique(gsCKnotVector * kv);


#ifdef __cplusplus
}
#endif
