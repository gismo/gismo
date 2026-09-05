
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCQuadRule* gsQuadrature_get(gsCBasis * basis, gsCOptionList * options, int fixDir);
    GISMO_EXPORT gsCQuadRule* gsGaussRule_create(int d, int* numNodes, int digits);
    GISMO_EXPORT gsCQuadRule* gsLobattoRule_create(int d, int* numNodes, int digits);
    GISMO_EXPORT void gsQuadRule_delete(gsCQuadRule * qr);

    GISMO_EXPORT int gsQuadRule_numNodes(gsCQuadRule * qr);
    GISMO_EXPORT int gsQuadRule_dim(gsCQuadRule * qr);

    GISMO_EXPORT void gsQuadRule_mapTo(gsCQuadRule* qr, gsCVector* lower, gsCVector* upper, gsCMatrix* nodes, gsCVector* weights);
    GISMO_EXPORT void gsQuadRule_mapToScalar(gsCQuadRule* qr, double startVal, double endVal, gsCMatrix* nodes, gsCVector* weights);

#ifdef __cplusplus
}
#endif
