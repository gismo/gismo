######################################################################
## CMakeLists.txt --- turbob64
## This file is part of the G+Smo library. 
## 
## Author: Matthias Moller
######################################################################

project(turbob64)

set(TURBOB64_SOURCES
${gismo_externals}/turbob64/turbob64c.c
${gismo_externals}/turbob64/turbob64d.c
${gismo_externals}/turbob64/turbob64sse.c
${gismo_externals}/turbob64/turbob64avx2.c
)

set(TURBOB64_HEADERS
${gismo_externals}/turbob64/conf.h
${gismo_externals}/turbob64/turbob64.h
${gismo_externals}/turbob64/turbob64_.h
)

# Apply same configuration as G+Smo
include( gsConfig)

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_NONSTDC_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)64le")
  add_definitions(-D__SSSE3__)
endif()

add_library(${PROJECT_NAME} OBJECT ${TURBOB64_SOURCES} ${TURBOB64_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES
POSITION_INDEPENDENT_CODE ON
FOLDER "G+Smo extensions")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:turbob64>
    CACHE INTERNAL "Gismo extensions to be included")
