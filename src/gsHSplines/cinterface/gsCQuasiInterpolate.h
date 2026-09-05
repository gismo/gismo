
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT void gsQuasiInterpolate_localIntpl_into( gsCBasis * basis,
                                                          gsCFunctionSet * fun,
                                                          gsCMatrix * result);

    GISMO_EXPORT void gsQuasiInterpolate_Taylor_into( gsCBasis * basis,
                                                      gsCFunctionSet * fun,
                                                      int deg,
                                                      gsCMatrix * result);

    GISMO_EXPORT void gsQuasiInterpolate_Schoenberg_into( gsCBasis * basis,
                                                          gsCFunctionSet * fun,
                                                          gsCMatrix * result);

#ifdef __cplusplus
}
#endif
