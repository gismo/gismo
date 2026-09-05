
#ifdef __cplusplus
extern "C"
{
#endif

#   define gsGeometryTransform_print gsFunctionSet_print
#   define gsGeometryTransform_delete gsFunctionSet_delete
    
    GISMO_EXPORT gsCGeometryTransform * gsGeometryTransform_create(gsCBasis* b, gsCMatrix * m,
                                                                   gsCVector * v);
    
#ifdef __cplusplus
}
#endif
