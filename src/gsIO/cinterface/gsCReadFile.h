#include <gsCInterface/gsCTypes.h>

#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT void* gsCReadFile (char* filename);
    GISMO_EXPORT void  gsCWriteFile(gsCFunctionSet * obj, char * filename);

#ifdef __cplusplus
}
#endif
