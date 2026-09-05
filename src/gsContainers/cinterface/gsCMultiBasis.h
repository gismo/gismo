
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create();
    GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create_basis(gsCBasis * basis);
    GISMO_EXPORT gsCMultiBasis * gsMultiBasis_create_patches(gsCMultiPatch* mp);

    GISMO_EXPORT void gsMultiBasis_addBasis(gsCMultiBasis* mp, gsCGeometry* geo);

    GISMO_EXPORT gsCBasis * gsMultiBasis_basis(gsCMultiBasis * mp, int i);

    GISMO_EXPORT gsCMultiBasis * gsMultiBasis_read(char* filename);
    GISMO_EXPORT void            gsMultiBasis_write(gsCMultiBasis * obj, char* filename);

#ifdef __cplusplus
}
#endif
