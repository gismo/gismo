######################################################################
## CMakeLists.txt --- gzstream
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
######################################################################

project(gzstream)

set(GZSTREAM_SOURCES ${gismo_externals}/gzstream/gzstream.cpp)
set(GZSTREAM_HEADERS${gismo_externals}/gzstream/gzstream.h)

# Apply same configuration as G+Smo
include( gsConfig)

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_NONSTDC_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_library(${PROJECT_NAME} OBJECT ${GZSTREAM_SOURCES} ${GZSTREAM_HEADERS})

set_target_properties(${PROJECT_NAME} 
PROPERTIES POSITION_INDEPENDENT_CODE ON
FOLDER "G+Smo extensions")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:gzstream>
    CACHE INTERNAL "Gismo extensions to be included")
