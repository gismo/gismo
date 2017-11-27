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

if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
   # Set default build type to RelWithDebInfo
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
   "Type of build (None Debug Release RelWithDebInfo MinSizeRel)" FORCE)
   if(NOT CMAKE_CONFIGURATION_TYPES)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
     "RelWithDebInfo" "MinSizeRel")
   endif()
endif()

if (GISMO_BUILD_PCH)
project(gismo C CXX CXXPCH)
else()
project(gismo C CXX)
endif()

if(NOT GISMO_INSOURCE_BUILD AND "${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. \n How to fix this:\n 1. Type \nmake cleancmake\nto remove garbage CMake files. \n\n 2. Please make a separate folder for building, otherwise typing \nmake\nthat will create a ./build folder and will compile in that folder. ")
endif()

## #################################################################
## Cmake Inclusions and configuration options
## #################################################################

include(CTest)
enable_testing()

#Use solution folders for Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Remove NDEBUG flag from RelWithDebInfo builds
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

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
message( "                          version ${${PROJECT_NAME}_VERSION}" )

## #################################################################
## Project options
## #################################################################

include( gsOptions )

## #################################################################
## Preprocessor definitions
## #################################################################

# Data directory
include ( gsDataDir )

## #################################################################
## JIT compiler
## #################################################################

include(gsJITConfigXml)
gsJITConfigXml(${CMAKE_CURRENT_SOURCE_DIR}/cmake/jit.xml.in
  ${CMAKE_BINARY_DIR}/config/jit.xml)

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
unset(${PROJECT_NAME}_LINKER     CACHE)
unset(gismo_TMPL_HEADERS         CACHE)
set(${PROJECT_NAME}_optionals ""
    CACHE INTERNAL "${PROJECT_NAME} list of optionals")

set (GISMO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src"
           "${PROJECT_SOURCE_DIR}/external"
      "${PROJECT_SOURCE_DIR}/extensions"
     CACHE INTERNAL "${PROJECT_NAME} include directories" )

include_directories(${GISMO_INCLUDE_DIRS})

## #################################################################
## Extensions
## #################################################################

# Add externals directory
add_subdirectory(external)

if(GISMO_WITH_CODIPACK)
  add_subdirectory(extensions/gsCoDiPack)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${CODIPACK_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  include_directories(${CODIPACK_INCLUDE_DIR})
endif(GISMO_WITH_CODIPACK)

#include_directories(${METIS_INCLUDE_DIR})
#set(gismo_LINKER ${gismo_LINKER} ${METIS_LIBRARIES}
#    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

if(GISMO_WITH_PSOLID)
  add_subdirectory(extensions/gsParasolid)
endif(GISMO_WITH_PSOLID)

if(GISMO_WITH_ONURBS)
  add_subdirectory(extensions/gsOpennurbs)
endif(GISMO_WITH_ONURBS)

if(GISMO_WITH_FDBB)
  add_subdirectory(extensions/gsFdbb)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${FDBB_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  include_directories(${FDBB_INCLUDE_DIR})
endif(GISMO_WITH_FDBB)

#if (GISMO_WITH_VTK)
#  add_subdirectory(extensions/gsVtk)
#  add_definitions(-DGISMO_WITH_VTK="${GISMO_WITH_VTK}")
#endif(GISMO_WITH_VTK)

if(GISMO_WITH_IPOPT)
  add_subdirectory(extensions/gsIpopt)
  #set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIR}
  #CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${IPOPT_INCLUDE_DIR})
endif(GISMO_WITH_IPOPT)

if(GISMO_WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${SUPERLU_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${SUPERLU_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(SYSTEM ${SUPERLU_INCLUDES})
endif(GISMO_WITH_SUPERLU)

if(GISMO_WITH_TAUCS)
  find_package(Taucs REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TAUCS_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${TAUCS_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(SYSTEM ${SUPERLU_INCLUDES})
endif(GISMO_WITH_TAUCS)

if(GISMO_WITH_UMFPACK)
  find_package(Umfpack REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UMFPACK_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${UMFPACK_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(SYSTEM ${UMFPACK_INCLUDES})
endif(GISMO_WITH_UMFPACK)

#if(GISMO_WITH_METIS)

if(GISMO_WITH_MPI)
  find_package(MPI REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPI_INCLUDE_PATH}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${MPI_CXX_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(GISMO_WITH_MPI)

if(GISMO_WITH_TRILINOS)
  add_subdirectory(extensions/gsTrilinos)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TRILINOS_INCLUDE_DIR} 
    CACHE INTERNAL "${PROJECT_NAME} include directories") 
  include_directories(${TRILINOS_INCLUDE_DIR})
endif(GISMO_WITH_TRILINOS)

## #################################################################
## Build G+Smo library instantized for GISMO_COEFF_TYPE number type
## #################################################################

# external inclusion paths
#get_property(EXTERNAL_INCLUDE_DIRS GLOBAL PROPERTY EXTERNAL_INCLUDE_DIRS)
#include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external)

# Add main sources directory
add_subdirectory(src)

if(GISMO_BUILD_LIB)

  if (GISMO_EXTRA_DEBUG)
    set(gismo_SOURCES ${gismo_SOURCES} ${PROJECT_SOURCE_DIR}/src/misc/gsStackWalker.cpp)
  endif()

endif()

include(gsLibrary)

## #################################################################
## Documentation
## #################################################################

add_subdirectory(doc/snippets EXCLUDE_FROM_ALL)
add_subdirectory(doc)

## #################################################################
## Plugins
## #################################################################

add_custom_target(plugins)

if(GISMO_BUILD_AXL)
  add_subdirectory(plugins/gsAxel)
endif(GISMO_BUILD_AXL)

if(GISMO_BUILD_RHINOPLUGINS)
  add_subdirectory(plugins/RhinoTHBSplines)
endif(GISMO_BUILD_RHINOPLUGINS)

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
else()
  add_subdirectory(examples EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

## #################################################################
## Unittests
## #################################################################

if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(unittests)
else()
  add_subdirectory(unittests EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_UNITTESTS)


add_subdirectory(optional)

if (NOT gismo_optionals STREQUAL "")
  string (REPLACE ";" ", " gismo_optionals "${gismo_optionals}")
  message(STATUS "Optional compile list: ${gismo_optionals}")
endif()
