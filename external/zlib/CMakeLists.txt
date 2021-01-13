######################################################################
## CMakeLists.txt --- zlib
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
######################################################################

project(zlib)

#set(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} -H")

## #################################################################
## Collect files
## #################################################################

set(ZLIB_SOURCES
${gismo_externals}/zlib/adler32.c
${gismo_externals}/zlib/compress.c
${gismo_externals}/zlib/crc32.c
${gismo_externals}/zlib/deflate.c
${gismo_externals}/zlib/gzclose.c
${gismo_externals}/zlib/gzlib.c
${gismo_externals}/zlib/gzread.c
${gismo_externals}/zlib/gzwrite.c
${gismo_externals}/zlib/infback.c
${gismo_externals}/zlib/inffast.c
${gismo_externals}/zlib/inflate.c
${gismo_externals}/zlib/inftrees.c
${gismo_externals}/zlib/trees.c
${gismo_externals}/zlib/uncompr.c
${gismo_externals}/zlib/zutil.c
)

set(ZLIB_HEADERS
${gismo_externals}/zlib/crc32.h
${gismo_externals}/zlib/deflate.h
${gismo_externals}/zlib/gzguts.h
${gismo_externals}/zlib/inffast.h
${gismo_externals}/zlib/inffixed.h
${gismo_externals}/zlib/inflate.h
${gismo_externals}/zlib/inftrees.h
${gismo_externals}/zlib/trees.h
${gismo_externals}/zlib/zlib.h
${gismo_externals}/zlib/zconf.h
${gismo_externals}/zlib/zutil.h
)

# Apply same configuration as G+Smo
include( gsConfig)

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_NONSTDC_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

#if(CMAKE_C_COMPILER_ID MATCHES "Clang")
#    STRING(REPLACE "-Wimplicit-function-declaration" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
#endif()

add_definitions(-DZ_PREFIX)   # use prefixed symbols
if(GISMO_WITH_ONURBS)
  add_definitions(-DMY_ZCALLOC) # use oncalloc() and onfree()
endif()

add_library(${PROJECT_NAME} OBJECT ${ZLIB_SOURCES} ${ZLIB_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES
COMPILE_DEFINITIONS "Z_PREFIX"
POSITION_INDEPENDENT_CODE ON
LINKER_LANGUAGE CXX
FOLDER "G+Smo extensions")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:zlib>
    CACHE INTERNAL "Gismo extensions to be included")
