######################################################################
### CMakeLists.txt --- gsOpennurbs
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

## OpenNurbs extension
project(gsHLBFGS)

set(HLBFGS_SOURCES
${gismo_externals}/HLBFGS/HLBFGS.cpp
${gismo_externals}/HLBFGS/HLBFGS_BLAS.cpp
${gismo_externals}/HLBFGS/ICFS.cpp
${gismo_externals}/HLBFGS/LineSearch.cpp
${gismo_externals}/HLBFGS/Lite_Sparse_Matrix.cpp
)

SET(HLBFGS_HEADERS
${gismo_externals}/HLBFGS/HLBFGS.h
${gismo_externals}/HLBFGS/HLBFGS_BLAS.h
${gismo_externals}/HLBFGS/ICFS.h
${gismo_externals}/HLBFGS/LineSearch.h
${gismo_externals}/HLBFGS/Lite_Sparse_Matrix.h
)

# Apply same configuration as G+Smo
include(gsConfig)
set(CMAKE_CXX_VISIBILITY_PRESET default)

add_library(${PROJECT_NAME} OBJECT ${HLBFGS_SOURCES}) #STATIC
#target_link_libraries(${PROJECT_NAME} Rpcrt4) # Note: linked with gismolib instead

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/W0")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-w")
endif()
#In CMake 2.8.12 the above can be:
#target_compile_options(${PROJECT_NAME} PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-w> $<$<CXX_COMPILER_ID:MSVC>: /W0>)

set_target_properties(${PROJECT_NAME}
PROPERTIES POSITION_INDEPENDENT_CODE ON
COMPILE_DEFINITIONS gismo_EXPORTS
FOLDER "G+Smo extensions")

#SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION "2010.04.95")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")
# if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" )
# set(gismo_LINKER ${gismo_LINKER} Rpcrt4
#     CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
# endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "*.h")
