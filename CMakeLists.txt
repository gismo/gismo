#####################################################################
## CMakeLists.txt --- root
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

cmake_minimum_required(VERSION 2.8.12)

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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(NOT CMAKE_BUILD_TYPE)
   # Set default build type to Release
   set(CMAKE_BUILD_TYPE Release CACHE STRING
   "Type of build (None Debug Release RelWithDebInfo MinSizeRel)" FORCE)
   if(NOT CMAKE_CONFIGURATION_TYPES)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
     "RelWithDebInfo" "MinSizeRel")
   endif()
endif()

set(gismo_VERSION_MAJOR 20) #year
set(gismo_VERSION_MINOR 12) #month
set(gismo_VERSION_PATCH 0 ) #patch
set(gismo_VERSION
  "${gismo_VERSION_MAJOR}.${gismo_VERSION_MINOR}.${gismo_VERSION_PATCH}" CACHE INTERNAL "gismo version number")
set(GISMO_VERSION "${gismo_VERSION}" CACHE INTERNAL "G+Smo version number")

if(CMAKE_QUIET)
  function(message)
    list(GET ARGV 0 MessageType)
    if(MessageType STREQUAL FATAL_ERROR OR
	MessageType STREQUAL SEND_ERROR OR
	MessageType STREQUAL WARNING OR
	MessageType STREQUAL AUTHOR_WARNING)
      list(REMOVE_AT ARGV 0)
      _message(${MessageType} "${ARGV}")
    endif()
  endfunction()
  unset(CMAKE_QUIET CACHE)
endif()

if(POLICY CMP0048)# CMake 3.0
  project(gismo LANGUAGES C CXX  VERSION ${gismo_VERSION})
else()
  project(gismo C CXX)
endif()

if (GISMO_BUILD_PCH)
  enable_language(CXXPCH)
endif()

if(NOT GISMO_INSOURCE_BUILD AND "${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. \nHow to fix this:\n1. Type \nmake cleancmake\nto remove garbage CMake files.\n2. Please make a separate folder for building, otherwise typing \nmake\nthat will create a ./build folder and will compile in that folder. ")
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
file(WRITE "${PROJECT_BINARY_DIR}/CTestCustom.cmake" "set(CTEST_SOURCE_DIRECTORY ${PROJECT_SOURCE_DIR})\nset(CTEST_BINARY_DIRECTORY   ${PROJECT_BINARY_DIR})\nset(CTEST_SITE         ${SITE})\nset(CTEST_BUILD_NAME   ${BUILDNAME})\n")

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

## #################################################################
## Collect headers and sources
## #################################################################

# Clean previously set variables in the cache
unset(GISMO_EXTRA_INSTANCE       CACHE)
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

# Add externals directory
add_subdirectory(external)

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

if(GISMO_WITH_MPFR OR GISMO_WITH_GMP)
  add_subdirectory(extensions/gsMultiPrecision)
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
endif(GISMO_WITH_UNUM)

if(GISMO_WITH_OCC)
  add_subdirectory(extensions/gsOpenCascade)
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

## #################################################################
## Unittests
## #################################################################

if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(unittests)
else()
  add_custom_target(unittests COMMAND "" COMMENT "Set CMake argument GISMO_BUILD_UNITTESTS=ON to enable unittests")
endif(GISMO_BUILD_UNITTESTS)

## #################################################################
## G+Smo sources
## #################################################################

# Add main sources directory
add_subdirectory(src)

## #################################################################
## Submodules
## #################################################################

option(GISMO_COMPFLOW "Enable compflow module" OFF)
if(GISMO_COMPFLOW)
  include(gsFetch)
  gismo_fetch_module(gsCompFlow)
endif()

option(GISMO_ELASTICITY "Enable elasticity module" OFF)
if(GISMO_ELASTICITY)
  include(gsFetch)
  gismo_fetch_module(gsElasticity)
endif()

option(GISMO_EXASTENCILS "Enable exastenciles module" OFF)
if(GISMO_EXASTENCILS)
  include(gsFetch)
  gismo_fetch_module(gsExaStencils)
endif()

option(GISMO_STRUCTURALANALYSIS "Enable Structural Analysis module" OFF)
if(GISMO_STRUCTURALANALYSIS)
  include(gsFetch)
  gismo_fetch_module(gsStructuralAnalysis)
endif()

option(GISMO_KLSHELL "Enable Kirchhoff-Love Shell" OFF)
if(GISMO_KLSHELL)
  include(gsFetch)
  gismo_fetch_module(gsKLShell)
endif()

option(GISMO_UNSUPPORTED "Enable unsupported module" OFF)
option(GISMO_MOTOR "Enable motor module" OFF)
if(GISMO_UNSUPPORTED OR GISMO_MOTOR)
  include(gsFetch)
  gismo_fetch_module(unsupported)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GISMO_DEV_INCLUDE_DIRS}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
endif()

if(GISMO_MOTOR)
  include(gsFetch)
  gismo_fetch_module(motor)
endif()

## #################################################################
## Shared and static library
## #################################################################

include(gsLibrary)

## #################################################################
## Documentation
## #################################################################

if(GISMO_BUILD_UNITTESTS OR NOSNIPPETS)
  add_subdirectory(doc/snippets EXCLUDE_FROM_ALL)
else()
  add_subdirectory(doc/snippets)
endif()
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

# Configure
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsConfig.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsConfig.h" )
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsExport.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsExport.h" )

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
