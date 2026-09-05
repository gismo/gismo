
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCGeometry * gsGeometry_read(char* filename);
    GISMO_EXPORT void          gsGeometry_write(gsCGeometry * obj, char* filename);

#   define gsGeometry_print gsFunctionSet_print
#   define gsGeometry_delete gsFunctionSet_delete
#   define gsBSpline_delete gsFunctionSet_delete






    GISMO_EXPORT gsCGeometry* gsGeometry_clone(gsCGeometry * g);

    GISMO_EXPORT gsCBasis* gsGeometry_basis(gsCGeometry * g);

    GISMO_EXPORT void gsGeometry_coefs_into(gsCGeometry * g, gsCMatrix * coef);
    GISMO_EXPORT void gsGeometry_setCoefs(gsCGeometry * g, gsCMatrix * coef);


    GISMO_EXPORT void gsGeometry_uniformRefine(gsCGeometry * b, int numKnots, int mul, int dir);
    GISMO_EXPORT void gsGeometry_refineElements(gsCGeometry * b, int * boxData, int boxSize);
    GISMO_EXPORT void gsGeometry_refine(gsCGeometry * b, gsCMatrix * boxes, int refExt);

    GISMO_EXPORT void gsGeometry_degreeElevate(gsCGeometry * b, int i, int dir);

    GISMO_EXPORT void gsGeometry_recoverPoints(gsCGeometry * g, gsCMatrix * uv, gsCMatrix * xyz,
                                                int k, double accuracy);
    GISMO_EXPORT void gsGeometry_recoverPointGrid(gsCGeometry * g, gsCVector * a, gsCVector * b,
                                                    gsCVectorInt * sz, gsCMatrix * xyz,
                                                    gsCMatrix * uv, int c, double accuracy);

    GISMO_EXPORT void gsGeometry_invertPointGrid(gsCGeometry * g, gsCVector * a, gsCVector * b,
                                       gsCVectorInt * sz,  gsCMatrix * result, double accuracy);
    GISMO_EXPORT void gsGeometry_invertPoints(gsCGeometry * g, gsCMatrix * points, gsCMatrix * result,
                                    double accuracy);

    GISMO_EXPORT void gsGeometry_normal_into(gsCGeometry * fs,
                                 gsCMatrix * u,
                                 gsCMatrix * result);


    GISMO_EXPORT double gsGeometry_closestPointTo(gsCGeometry * fs,
                                       gsCMatrix * pt,
                                       gsCMatrix * result,
                                       double accuracy);


#ifdef __cplusplus
}
#endif
