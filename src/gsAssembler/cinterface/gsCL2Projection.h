
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT double gsL2Projection_into( gsCFunctionSet * projectionBasis,
                                             gsCMultiBasis * integrationBasis,
                                             gsCMultiPatch * geometryMap,
                                             gsCFunctionSet * sourceFunction,
                                             gsCMatrix * coefs,
                                             gsCOptionList * options);

#ifdef __cplusplus
}
#endif
