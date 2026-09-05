
#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCOptionList* gsOptionList_create();
    GISMO_EXPORT void gsOptionList_delete(gsCOptionList * ol);

    GISMO_EXPORT void gsOptionList_print(gsCOptionList * ol);

    GISMO_EXPORT void gsOptionList_addString(gsCOptionList * ol, const char * label, const char * description, const char * value);
    GISMO_EXPORT void gsOptionList_addInt(gsCOptionList * ol, const char * label, const char * description, int value);
    GISMO_EXPORT void gsOptionList_addReal(gsCOptionList * ol, const char * label, const char * description, double value);
    GISMO_EXPORT void gsOptionList_addSwitch(gsCOptionList * ol, const char * label, const char * description, int value);

    GISMO_EXPORT void gsOptionList_setString(gsCOptionList * ol, const char * label, const char * value);
    GISMO_EXPORT void gsOptionList_setInt(gsCOptionList * ol, const char * label, int value);
    GISMO_EXPORT void gsOptionList_setReal(gsCOptionList * ol, const char * label, double value);
    GISMO_EXPORT void gsOptionList_setSwitch(gsCOptionList * ol, const char * label, int value);

    GISMO_EXPORT char*  gsOptionList_getString(gsCOptionList * ol, const char * label);
    GISMO_EXPORT int    gsOptionList_getInt(gsCOptionList * ol, const char * label);
    GISMO_EXPORT double gsOptionList_getReal(gsCOptionList * ol, const char * label);
    GISMO_EXPORT int    gsOptionList_getSwitch(gsCOptionList * ol, const char * label);

#ifdef __cplusplus
}
#endif
