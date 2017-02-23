### gsDataDir.cmake ---
## 
## set GISMO_DATA_DIR based on compiler-toolchain
## 
## Author: Harald Weiner
## Copyright (C) 2015 - RICAM-Linz.
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

