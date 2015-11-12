######################################################################
## CMakeLists.txt
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012-2015 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 OLD)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

project(gismo)

if(NOT GISMO_INSOURCE_BUILD AND "${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. \n How to fix this:\n 1. Please make a separate folder for building, otherwise type \nmake\nthat will create a ./build folder and will compile in that folder. Cmake has created garbage files/dirs (you may manually remove them):\n CMakeCache.txt CMakeFiles \n 2. If your source folder contains these two files please remove them first!")
endif()

## #################################################################
## Cmake Inclusions and configuration options
## #################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CTest)
enable_testing()

include( gsVersion )
include( gismoUse  )
include( gsConfig )

message( "   GGGGGGGGG      GGGG      GGGGGGGGG  GGGGGG   GGGGGG  GGGGGGGGGG" )
message( "  GGGG            GGGG     GGGG        GGGGGG  GGGGGG  GGGG   GGGG" )
message( " GGGG         GGGGGGGGGGGG GGGGGGGGG   G GGGG  G GGGG GGGG    GGGG" )
message( " GGGG GGGGGG GGGGGGGGGGGGG GGGGGGGGGG GG GGGG GG GGGG GGGG   GGGGG" )
message( "GGGGG  GGGGG GGGGGGGGGGGG  GGGGGGGGG  GG GGGGGG GGGG  GGGG   GGGG " )
message( "GGGG   GGGG      GGGG           GGGG  GG  GGGG  GGGG  GGGG   GGGG " )
message( " GGGGGGGGGG      GGGG     GGGGGGGGG  GG   GGG   GGGG  GGGGGGGGGG  " )
#message( "                          version ${${PROJECT_NAME}_VERSION}" )
message( "                  version ${${PROJECT_NAME}_VERSION} (revision ${GISMO_WC_REVISION})" )

## #################################################################
## Project options
## #################################################################

include( gsOptions )

## #################################################################
## Preprocessor definitions
## #################################################################

# Data directory
include ( gsDataDir )

# Configure
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsConfig.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsConfig.h" )
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsExport.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsExport.h" )
file(COPY "${PROJECT_SOURCE_DIR}/CTestConfig.cmake"
	  DESTINATION "${PROJECT_BINARY_DIR}")

include_directories(${PROJECT_BINARY_DIR})

## #################################################################
## Collect headers and sources
## #################################################################

# Clean previously set variables in the cache
unset(${PROJECT_NAME}_HEADERS    CACHE)
unset(${PROJECT_NAME}_SOURCES    CACHE)
unset(${PROJECT_NAME}_EXTENSIONS CACHE)
unset(${PROJECT_NAME}_MODULES    CACHE)
unset(${PROJECT_NAME}_optionals  CACHE)
unset(gismo_TMPL_HEADERS         CACHE)
set(${PROJECT_NAME}_optionals ""
    CACHE INTERNAL "${PROJECT_NAME} list of optionals")

set (GISMO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src;${PROJECT_SOURCE_DIR}/external;${PROJECT_SOURCE_DIR}/extensions"
     CACHE INTERNAL "${PROJECT_NAME} include directories")

# Add main sources directory
add_subdirectory(src)

# Add externals directory
add_subdirectory(external)

# Inclusion paths
include_directories(${GISMO_INCLUDE_DIRS})

## #################################################################
## Extensions
## #################################################################

if(GISMO_WITH_MPFR)
  add_definitions(-DGISMO_WITH_MPFR="${GISMO_WITH_MPFR}")
endif(GISMO_WITH_MPFR)

if(GISMO_WITH_PSOLID)
  add_subdirectory(extensions/gsParasolid)
endif(GISMO_WITH_PSOLID)

if(GISMO_WITH_ONURBS)
  add_subdirectory(extensions/gsOpennurbs)
endif(GISMO_WITH_ONURBS)

#if (GISMO_WITH_VTK)
#  add_subdirectory(extensions/gsVtk)
#  add_definitions(-DGISMO_WITH_VTK="${GISMO_WITH_VTK}")
#endif(GISMO_WITH_VTK)

if(GISMO_WITH_IPOPT)
  add_subdirectory(extensions/gsIpopt)
  include_directories(${IPOPT_INCLUDE_DIR})
endif(GISMO_WITH_IPOPT)

if(GISMO_WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
endif(GISMO_WITH_SUPERLU)

## #################################################################
## Build G+Smo library instantized for GISMO_COEFF_TYPE number type
## #################################################################

if(GISMO_BUILD_LIB)

  if (GISMO_EXTRA_DEBUG)
    set(gismo_SOURCES ${gismo_SOURCES} ${PROJECT_SOURCE_DIR}/src/misc/gsStackWalker.cpp)
  endif()

endif()

include(gsLibrary)

## #################################################################
## Documentation
## #################################################################

add_subdirectory(doc)

## #################################################################
## Plugins
## #################################################################

if(GISMO_BUILD_AXL)
  add_subdirectory(plugins/gsAxel)
endif(GISMO_BUILD_AXL)

if(GISMO_BUILD_MEX)
  add_subdirectory(plugins/gsMex)
endif(GISMO_BUILD_MEX)


## #################################################################
## Install
## #################################################################

include( gsInstall )

## #################################################################
## Executables
## #################################################################

# Note: if the shared libs are created, examples are linked against the
#  shared lib, otherwise targets are not linked against anything, only
#  templated code is used.
if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(GISMO_BUILD_EXAMPLES)

add_subdirectory(optional)

if (NOT gismo_optionals STREQUAL "")
  string (REPLACE ";" ", " gismo_optionals "${gismo_optionals}")
  message(STATUS "Optional compile list: ${gismo_optionals}")
endif()

