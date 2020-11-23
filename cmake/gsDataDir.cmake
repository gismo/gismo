######################################################################
## gsDataDir.cmake
## 
## set GISMO_DATA_DIR based on compiler-toolchain
##
## Author: Harald Weiner
######################################################################

# Data directory
if (DEFINED EMSCRIPTEN)
	set(GISMO_DATA_DIR "filedata/" )
	message("GISMO_DATA_DIR set to relative path '${GISMO_DATA_DIR}'")
	set(GISMO_ORIG_DATA_DIR "${PROJECT_SOURCE_DIR}/filedata/")
else()
	set (GISMO_DATA_DIR   "${PROJECT_SOURCE_DIR}/filedata/")
endif()

# Config-files directory
set(GISMO_CONFIG_DIR "${CMAKE_BINARY_DIR}/config/")

# Set default search paths
if (NOT DEFINED GISMO_SEARCH_PATHS)
	set(GISMO_SEARCH_PATHS "${GISMO_DATA_DIR}" CACHE STRING
		"Define paths where files should be searched; seperated by semicolon.")
endif()
