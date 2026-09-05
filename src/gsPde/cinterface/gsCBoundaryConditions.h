
#ifdef __cplusplus
extern "C"
{
#endif

    #include<stdbool.h>

    GISMO_EXPORT gsCBoundaryConditions * gsBoundaryConditions_create();

    GISMO_EXPORT void gsBoundaryConditions_addCondition(gsCBoundaryConditions * bc,
                                                        int patch,
                                                        int side,
                                                        int type,
                                                        gsCFunctionSet * fun,
                                                        int unknown,
                                                        int component,
                                                        bool parametric);

    GISMO_EXPORT void gsBoundaryConditions_addCornerValue(gsCBoundaryConditions * bc,
                                                          int corner,
                                                          double value,
                                                          int patch,
                                                          int unknown,
                                                          int component);

    GISMO_EXPORT void gsBoundaryConditions_setGeoMap(gsCBoundaryConditions * bc, gsCFunctionSet * gm);

    GISMO_EXPORT void gsBoundaryConditions_print(gsCBoundaryConditions * bc);

    GISMO_EXPORT void gsBoundaryConditions_delete(gsCBoundaryConditions * bc);


#ifdef __cplusplus
}
#endif
