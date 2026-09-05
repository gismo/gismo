
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineUnitInterval(int);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineRectangle(double, double, double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineTrapezium(double, double, double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineSquare(double, double, double);

    GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineSquareGrid(int, int, double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineCube(double, double, double, double);

    GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineCubeGrid(int, int, int, double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsQuarterAnnulus(double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsAnnulus(double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineSaddle();

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsSphere(double, double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_NurbsCircle(double, double, double);

    GISMO_EXPORT gsCGeometry* gsNurbsCreator_BSplineTriangle(double, double);

    GISMO_EXPORT gsCMultiPatch* gsNurbsCreator_BSplineStar(int, double, double);

#ifdef __cplusplus
}
#endif
