
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCFitting* gsFitting_create(gsCMatrix * param_values, gsCMatrix * points, gsCBasis * basis);
    GISMO_EXPORT void gsFitting_delete(gsCFitting * fitter);

    GISMO_EXPORT void gsFitting_compute(gsCFitting* fitter, double lambda);
    GISMO_EXPORT void gsFitting_parameterCorrection(gsCFitting* fitter, double accuracy, int maxIter, double tolOrth);

    GISMO_EXPORT void gsFitting_computeErrors(gsCFitting* fitter);

    GISMO_EXPORT double  gsFitting_minPointError(gsCFitting* fitter);
    GISMO_EXPORT double  gsFitting_maxPointError(gsCFitting* fitter);
    GISMO_EXPORT double* gsFitting_pointWiseErrors(gsCFitting* fitter);
    GISMO_EXPORT int     gsFitting_numPointsBelow(gsCFitting* fitter, double threshold);

    GISMO_EXPORT gsCGeometry* gsFitting_result(gsCFitting* fitter);

#ifdef __cplusplus
}
#endif
