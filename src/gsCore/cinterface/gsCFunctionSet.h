
#ifdef __cplusplus
extern "C"
{
#endif
    GISMO_EXPORT void gsFunctionSet_delete(gsCFunctionSet * ptr);

    GISMO_EXPORT void gsFunctionSet_print(gsCFunctionSet * fs);
    GISMO_EXPORT int gsFunctionSet_domainDim(gsCFunctionSet * fs);
    GISMO_EXPORT int gsFunctionSet_targetDim(gsCFunctionSet * fs);

    GISMO_EXPORT int gsFunctionSet_nPieces(gsCFunctionSet * fs);

    GISMO_EXPORT gsCMatrix* gsFunctionSet_support(gsCFunctionSet * fs);

    GISMO_EXPORT void gsFunctionSet_eval_into(gsCFunctionSet * fs,
                                              gsCMatrix * u,
                                              gsCMatrix * result);

    GISMO_EXPORT void gsFunctionSet_deriv_into(gsCFunctionSet * fs,
                                               gsCMatrix * u,
                                               gsCMatrix * result);

    GISMO_EXPORT void gsFunctionSet_deriv2_into(gsCFunctionSet * fs,
                                                gsCMatrix * u,
                                                gsCMatrix * result);

    // gsFunctionSet::evalAllDers_into(u, 1, result,false);
    GISMO_EXPORT void gsFunctionSet_evalAllDers1_into(gsCFunctionSet * fs,
                                                      gsCMatrix * u,
                                                      gsCMatrix * values,
                                                      gsCMatrix * deriv);

    // gsFunctionSet::evalAllDers_into(u, 2, result,false);
    GISMO_EXPORT void gsFunctionSet_evalAllDers2_into(gsCFunctionSet * fs,
                                                      gsCMatrix * u,
                                                      gsCMatrix * values,
                                                      gsCMatrix * deriv,
                                                      gsCMatrix * deriv2);

#ifdef __cplusplus
}
#endif
