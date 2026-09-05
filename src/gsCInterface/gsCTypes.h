// Opaque types that we use as a handles
#ifdef __cplusplus
extern "C"
{
#endif
    struct gsCMatrix;
    typedef struct gsCMatrix gsCMatrix;

    struct gsCMatrixInt;
    typedef struct gsCMatrixInt gsCMatrixInt;

    struct gsCVector;
    typedef struct gsCVector gsCVector;

    struct gsCVectorInt;
    typedef struct gsCVectorInt gsCVectorInt;

    struct gsCSparseMatrix;
    typedef struct gsCSparseMatrix gsCSparseMatrix;

    struct gsCFunctionSet; // Opaque type that we use as a handle
    typedef struct gsCFunctionSet gsCFunctionSet;

    typedef struct gsCFunctionSet gsCMultiPatch;
    typedef struct gsCFunctionSet gsCMultiBasis;
    typedef struct gsCFunctionSet gsCBasis;
    typedef struct gsCFunctionSet gsCGeometry;
    typedef struct gsCFunctionSet gsCGeometryTransform;
    typedef struct gsCFunctionSet gsCFunctionExpr;

    struct gsCKnotVector;
    typedef struct gsCKnotVector gsCKnotVector;

    struct gsCBoundaryConditions;
    typedef struct gsCBoundaryConditions gsCBoundaryConditions;

    struct gsCQuadRule;
    typedef struct gsCQuadRule gsCQuadRule;

    struct gsCOptionList;
    typedef struct gsCOptionList gsCOptionList;

    struct gsCFitting;
    typedef struct gsCFitting gsCFitting;

    //struct gsCGridIteratorCube;
    //typedef struct gsCGridIteratorCube gsCGridIteratorCube;

#ifdef __cplusplus
}
#endif
