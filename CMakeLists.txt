#####################################################################
## CMakeLists.txt
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012-2018 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 NEW)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

if(POLICY CMP0074)# CMake 3.12
cmake_policy(SET CMP0074 NEW)
endif()

if(POLICY CMP0077)# CMake 3.13
cmake_policy(SET CMP0077 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
   # Set default build type to Release
   set(CMAKE_BUILD_TYPE Release CACHE STRING
   "Type of build (None Debug Release RelWithDebInfo MinSizeRel)" FORCE)
   if(NOT CMAKE_CONFIGURATION_TYPES)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
     "RelWithDebInfo" "MinSizeRel")
   endif()
endif()

set(gismo_VERSION_MAJOR "0")
set(gismo_VERSION_MINOR "8")
set(gismo_VERSION_PATCH "5")
set(gismo_VERSION
  "${gismo_VERSION_MAJOR}.${gismo_VERSION_MINOR}.${gismo_VERSION_PATCH}" CACHE INTERNAL "gismo version number")

if(POLICY CMP0048)# CMake 3.0
  project(gismo LANGUAGES C CXX  VERSION ${gismo_VERSION})
else()
  project(gismo C CXX)
endif()

if (GISMO_BUILD_PCH)
  enable_language(CXXPCH)
endif()

if(NOT GISMO_INSOURCE_BUILD AND "${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. \n How to fix this:\n 1. Type \nmake cleancmake\nto remove garbage CMake files. \n\n 2. Please make a separate folder for building, otherwise typing \nmake\nthat will create a ./build folder and will compile in that folder. ")
endif()

## #################################################################
## Cmake Inclusions and configuration options
## #################################################################

set(CMAKE_INSTALL_MESSAGE NEVER)

#Use solution folders for Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Remove NDEBUG flag from RelWithDebInfo builds
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

include( gismoUse )
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
include(CTest)
# Note: clean cache to reset to OFF
set(BUILD_TESTING "Enable testing" ${GISMO_BUILD_UNITTESTS})

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

set (GISMO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/external"
     "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/extensions"
     CACHE INTERNAL "${PROJECT_NAME} include directories" )

# Offer the user the choice of overriding the installation directories
set(LIB_INSTALL_DIR     lib     CACHE STRING "Installation directory for libraries")
set(BIN_INSTALL_DIR     bin     CACHE STRING "Installation directory for executables")
set(INCLUDE_INSTALL_DIR include CACHE STRING "Installation directory for header files")

## #################################################################
## Extensions
## #################################################################

if (EIGEN_USE_MKL_ALL)
  # See http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html
  find_package(MKL REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MKL_INCLUDE_DIR}
    CACHE INTERNAL "${PROJECT_NAME} include directories")
endif()

#first time
include_directories(${GISMO_INCLUDE_DIRS})

include(gsFetch)

# Add externals directory
add_subdirectory(external)

if(GISMO_BUILD_UNITTESTS) # Unittests
  add_subdirectory(unittests)
else()
  add_custom_target(unittests COMMAND "" COMMENT "Set CMake argument GISMO_BUILD_UNITTESTS=ON to enable unittests")
endif(GISMO_BUILD_UNITTESTS)

if(GISMO_WITH_CODIPACK)
  add_subdirectory(extensions/gsCoDiPack)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${CODIPACK_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${CODIPACK_INCLUDE_DIR})
endif(GISMO_WITH_CODIPACK)

if(GISMO_WITH_IPOPT)
  add_subdirectory(extensions/gsIpopt)
  #set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIR}
  #CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${IPOPT_INCLUDE_DIR})
endif(GISMO_WITH_IPOPT)

#if(GISMO_WITH_METIS)
#include_directories(${METIS_INCLUDE_DIR})
#set(gismo_LINKER ${gismo_LINKER} ${METIS_LIBRARIES}
#    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
#endif(GISMO_WITH_METIS)

if(GISMO_WITH_MPI)
  find_package(MPI REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPI_INCLUDE_PATH}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${MPI_CXX_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(GISMO_WITH_MPI)

if(GISMO_WITH_MPFR OR GISMO_WITH_MPQ)
  find_package(GMP QUIET)
  find_package(MPFR QUIET)

  if (NOT GMP_FOUND OR NOT MPFR_FOUND)
    add_subdirectory(extensions/gsMultiPrecision)
  endif()

  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GMP_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${GMP_LIBRARY} ${GMPXX_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${GMP_INCLUDE_DIR})

  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPFR_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${MPFR_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${MPFR_INCLUDE_DIR})
endif()

if(GISMO_WITH_ONURBS)
  add_subdirectory(extensions/gsOpennurbs)
endif(GISMO_WITH_ONURBS)

#if(GISMO_WITH_PARDISO)
#endif(GISMO_WITH_PARDISO)

#if(GISMO_WITH_PASTIX)
#endif(GISMO_WITH_PASTIX)

if(GISMO_WITH_PSOLID)
  add_subdirectory(extensions/gsParasolid)
endif(GISMO_WITH_PSOLID)

if(GISMO_WITH_SPECTRA)
  add_subdirectory(extensions/gsSpectra)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${SPECTRA_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${SPECTRA_INCLUDE_DIR})
endif(GISMO_WITH_SPECTRA)

if(GISMO_WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${SUPERLU_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${SUPERLU_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${SUPERLU_INCLUDES})
endif(GISMO_WITH_SUPERLU)

if(GISMO_WITH_TAUCS)
  find_package(Taucs REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TAUCS_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${TAUCS_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${SUPERLU_INCLUDES})
endif(GISMO_WITH_TAUCS)

if(GISMO_WITH_UMFPACK)
  find_package(Umfpack REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UMFPACK_INCLUDES}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_LINKER ${gismo_LINKER} ${UMFPACK_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  #include_directories(SYSTEM ${UMFPACK_INCLUDES})
endif(GISMO_WITH_UMFPACK)

if(GISMO_WITH_UNUM)
  add_subdirectory(extensions/gsUnum)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UNUM_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${UNUM_INCLUDE_DIR})
endif(GISMO_WITH_UNUM)

if(GISMO_WITH_OCC)
  add_subdirectory(extensions/gsOpenCascade)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${OCC_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${OCC_INCLUDE_DIR})
endif(GISMO_WITH_OCC)

if(GISMO_WITH_SMESH)
  include(external/surface_mesh.cmake)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${SMESH_INCLUDE_DIRS}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${SMESH_INCLUDE_DIR})
endif()

#second time
include_directories(${GISMO_INCLUDE_DIRS})

if(GISMO_WITH_TRILINOS)
  add_subdirectory(extensions/gsTrilinos)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TRILINOS_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
  #include_directories(${TRILINOS_INCLUDE_DIR})
endif(GISMO_WITH_TRILINOS)

#third time
include_directories(${GISMO_INCLUDE_DIRS})

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

## #################################################################
## Submodules
## #################################################################

option(GISMO_COMPFLOW "Enable compflow module" OFF)
if(GISMO_COMPFLOW)
  gismo_fetch_module(gsCompFlow)
endif()

option(GISMO_ELASTICITY "Enable elasticity module" OFF)
if(GISMO_ELASTICITY)
  gismo_fetch_module(gsElasticity)
endif()

option(GISMO_EXASTENCILS "Enable exastenciles module" OFF)
if(GISMO_EXASTENCILS)
  gismo_fetch_module(gsExaStencils)
endif()

option(GISMO_UNSUPPORTED "Enable unsupported module" OFF)
option(GISMO_MOTOR "Enable motor module" OFF)
if(GISMO_UNSUPPORTED OR GISMO_MOTOR)
  gismo_fetch_module(unsupported)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GISMO_DEV_INCLUDE_DIRS}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
endif()

if(GISMO_MOTOR)
  gismo_fetch_module(motor)
endif()

## #################################################################
## Shared and static library
## #################################################################

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

if(GISMO_BUILD_MEX)
  add_subdirectory(plugins/gsMex)
endif(GISMO_BUILD_MEX)

#if(GISMO_BUILD_PVIEW)
#  add_subdirectory(plugins/gsParaview)
#endif(GISMO_BUILD_PVIEW)

if(GISMO_BUILD_RHINOPLUGINS)
  add_subdirectory(plugins/RhinoTHBSplines)
endif(GISMO_BUILD_RHINOPLUGINS)

## #################################################################
## Install
## #################################################################

include(gsInstall)

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
## Misc
## #################################################################

add_subdirectory(optional)

if (NOT gismo_optionals STREQUAL "")
  string (REPLACE ";" ", " gismo_optionals "${gismo_optionals}")
  message(STATUS "Optional compile list: ${gismo_optionals}")
endif()

unset(GISMO_UNAME CACHE)
unset(GISMO_PASS CACHE)

#message("GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS}")
