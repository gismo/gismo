
#ifdef __cplusplus
extern "C"
{
#endif

    #include <stdbool.h>

    GISMO_EXPORT gsCBasis * gsBasis_read(char* filename);
    GISMO_EXPORT void       gsBasis_write(gsCBasis * obj, char* filename);

#   define gsBasis_print gsFunctionSet_print
#   define gsBasis_delete gsFunctionSet_delete





    //
    // Methods, gsBasis
    //
    GISMO_EXPORT gsCBasis* gsBasis_clone(gsCBasis * b);

    GISMO_EXPORT void gsBasis_active_into(gsCBasis * b, gsCMatrix * u, gsCMatrixInt * result);

    GISMO_EXPORT void gsBasis_evalSingle_into(gsCBasis * b, int i, gsCMatrix * u, gsCMatrix * result);
    GISMO_EXPORT void gsBasis_derivSingle_into(gsCBasis * b, int i, gsCMatrix * u, gsCMatrix * result);
    GISMO_EXPORT void gsBasis_deriv2Single_into(gsCBasis * b, int i, gsCMatrix * u, gsCMatrix * result);

    GISMO_EXPORT gsCBasis * gsBasis_component(gsCBasis * b, int dir);
    GISMO_EXPORT int gsBasis_degree(gsCBasis * b, int dir);
    GISMO_EXPORT int gsBasis_numElements(gsCBasis * b);
    GISMO_EXPORT int gsBasis_dim(gsCBasis * b);
    GISMO_EXPORT int gsBasis_size(gsCBasis * b);
    GISMO_EXPORT gsCMatrix* gsBasis_support(gsCBasis * b, int i);

    GISMO_EXPORT void gsBasis_uniformRefine(gsCBasis * b, int numKnots, int mul, int dir);
    GISMO_EXPORT void gsBasis_refineElements(gsCBasis * b, int * boxData, int boxSize);
    GISMO_EXPORT void gsBasis_refine(gsCBasis * b, gsCMatrix * boxes, int refExt);

    GISMO_EXPORT void gsBasis_degreeElevate(gsCBasis * b, int i, int dir);

    GISMO_EXPORT void gsBasis_boundary_into(gsCBasis * b, int side, gsCMatrixInt * result);
    GISMO_EXPORT void gsBasis_boundaryOffset_into(gsCBasis * b, int side, int offset, gsCMatrixInt * result);

    GISMO_EXPORT void gsBasis_elements_into(gsCBasis * b, gsCMatrix*);
    GISMO_EXPORT void gsBasis_elementsBdr_into(gsCBasis * b, int, gsCMatrix*);


    // TODO:
    // - DegreeElevate

    //
    // Methods, Other
    //




#ifdef __cplusplus
}
#endif
