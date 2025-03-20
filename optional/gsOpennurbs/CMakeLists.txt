######################################################################
### CMakeLists.txt --- gsOpennurbs
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

## OpenNurbs extension
project(gsOpennurbs LANGUAGES CXX)

# Find the opennurbs library package
# find_package(opennurbs REQUIRED)

include(${CMAKE_SOURCE_DIR}/external/gsOpenNurbs.cmake)

# get_target_property(opennurbs_INCLUDE_DIR opennurbs INTERFACE_INCLUDE_DIRECTORIES)
# include_directories(${opennurbs_INCLUDE_DIR})

set(OPENNURBS_SOURCES
gsReadOpenNurbs.cpp
gsWriteOpenNurbs_.cpp
)

#if (BUILD_SHARED_LIBS)
#set(OPENNURBS_SOURCES ${OPENNURBS_SOURCES} ${gismo_externals}/onurbs/opennurbs_dll.cpp)
#endif()
#set(OPENNURBS_SOURCES_UNUSED
#${gismo_externals}/onurbs/opennurbs_gl.cpp
#  )

SET(OPENNURBS_HEADERS
gsReadOpenNurbs.h
gsWriteOpenNurbs.h
)

## TODO: CLB not sure what to do with this
set(OPENNURBS_UNUSED
${gismo_externals}/onurbs/opennurbs_gl.h
${gismo_externals}/onurbs/opennurbs_precompiledheader.cpp
  )

SET(ON_DLL_HEADERS
  ${gismo_externals}/onurbs/opennurbs_dll_resource.h
  )
########################

INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES(emmintrin.h HAVE_EMMINTRIN_H)
CHECK_INCLUDE_FILES(ieeefp.h HAVE_IEEEFP_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
if(HAVE_EMMINTRIN_H)
  add_definitions(-DHAVE_EMMINTRIN_H)
endif(HAVE_EMMINTRIN_H)
if(HAVE_IEEEFP_H)
  add_definitions(-DHAVE_IEEEFP_H)
endif(HAVE_IEEEFP_H)
if(HAVE_SYS_STAT_H)
  add_definitions(-DHAVE_SYS_STAT_H)
endif(HAVE_SYS_STAT_H)

# Apply same configuration as G+Smo
include(gsConfig)

#if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
#  add_definitions(
#    -DON_COMPILING_OPENNURBS
#    -DON_DLL_EXPORTS
#    -DON_DLL_IMPORTS
#    )
#  SET(OPENNURBS_HEADERS "${OPENNURBS_HEADERS};${ON_DLL_HEADERS}")
#else ()
add_definitions(-DON_COMPILING_OPENNURBS)
#endif()


if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  add_definitions(-DUNICODE)
  add_definitions(
    -wd4458 # Suppress 'declaration of '..' hides class member'	~2400 warns
    -wd4456 # declaration of '..' hides previous local declaration ~37 warns
    )
else()
  add_definitions(-DON_NO_WINDOWS) #for Mingw
endif()

if (GISMO_ZLIB_STATIC)
  #message("Using prefixed Zlib")
  add_definitions(-DZ_PREFIX)   # use prefixed zlib
  add_definitions(-DMY_ZCALLOC) # zlib uses oncalloc() and onfree()
endif()

add_library(${PROJECT_NAME} OBJECT ${OPENNURBS_SOURCES}) #STATIC

target_link_libraries(${PROJECT_NAME} PUBLIC opennurbs)

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
if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" )
set(gismo_LINKER ${gismo_LINKER} Rpcrt4
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "*.h")

if(TARGET gmp)
  add_dependencies(${PROJECT_NAME} gmp)
endif()
if (TARGET mpfr)
  add_dependencies(${PROJECT_NAME} mpfr)
endif()
