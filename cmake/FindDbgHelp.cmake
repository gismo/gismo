INCLUDE(FindPackageHandleStandardArgs)
#INCLUDE(HandleLibraryTypes)

FIND_PATH(DBGHELP_INCLUDE_DIR dbghelp.h
  PATHS ${DbgHelp_DIR} $ENV{DbgHelp_DIR}
  PATH_SUFFIXES include inc Include Inc
)
FIND_LIBRARY(DBGHELP_LIBRARY dbghelp.lib
  PATHS ${DbgHelp_DIR} $ENV{DXSDK_DIR}
  PATH_SUFFIXES lib Lib lib/x86 Lib/x86
  NO_DEFAULT_PATH # Or else CMake 2.6.0 will find the dll in system32 on windows
)
FIND_LIBRARY(DBGHELP_LIBRARY dbghelp.lib # search again if the library wasn't found in the dependencies package
  PATHS ${DbgHelp_DIR} $ENV{DXSDK_DIR}
  #PATH_SUFFIXES lib Lib lib/x64 Lib/x64
)

# Handle the REQUIRED argument and set DBGHELP_FOUND
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DbgHelp DEFAULT_MSG
  DBGHELP_LIBRARY
  DBGHELP_INCLUDE_DIR
)

MARK_AS_ADVANCED(
  DBGHELP_INCLUDE_DIR
  DBGHELP_LIBRARY
)
