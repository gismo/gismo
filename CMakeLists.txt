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

project(gismo)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. Please make a separate folder for building, otherwise type \nmake\nthat will create a ./build folder and will compile in that folder. Cmake has created garbage files/dirs (you may manually remove them):\n CMakeCache.txt CMakeFiles")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

## #################################################################
## Cmake Inclusions and options
## #################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( gsVersion )
include( gismoUse  )

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
## Configuration and Compilation
## #################################################################

include( gsConfig )

## #################################################################
## Project options
## #################################################################

include( gsOptions )

## #################################################################
## Preprocessor definitions
## #################################################################

# Shared pointer
find_package (TR1 QUIET)

# Data directory
set (GISMO_DATA_DIR   "${PROJECT_SOURCE_DIR}/filedata/")

# Configure
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsConfig.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsConfig.h" )
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsExport.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsExport.h" )
configure_file ("${PROJECT_SOURCE_DIR}/cmake/CTestConfig.cmake.in"
                "${CMAKE_BINARY_DIR}/CTestConfig.cmake" )

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

if (GISMO_WITH_VTK)
  add_subdirectory(extensions/gsVtk)
  add_definitions(-DGISMO_WITH_VTK="${GISMO_WITH_VTK}")
endif(GISMO_WITH_VTK)

if(GISMO_WITH_IPOPT)
  find_library(IpOpt_LIB_Blas NAMES blas libblas cblas atlcblas libatlcblas )
  find_library(IpOpt_LIB_Lapack NAMES lapack liblapack atlclapack libatlclapack )
  find_library(IpOpt_LIB NAMES ipopt libipopt )
  #message(WARNING, "Found IpOpt_LIB_IpOpt='${IpOpt_LIB_IpOpt}'")
  #message(WARNING, "Found IpOpt_LIB_IpOpt-NOTFOUND='${IpOpt_LIB_IpOpt-NOTFOUND}'")

  if(IpOpt_LIB)
      #message(WARNING, "Using ipopt system library")
      find_library(IpOpt_LIB_Mumps NAMES mumps libmumps cmumps libcmumps )
      find_library(IpOpt_LIB_Metis NAMES metis libmetis )
      set(IpOpt_INCLUDE_DIRS /usr/include/coin CACHE INTERNAL "IpOpt include path")
      set(IpOpt_LIB_IpOpt ${IpOpt_LIB})
  else(IpOpt_LIB)
      #message(WARNING, "Downloading and compiling ipopt")
      include(ExternalProject)
      ExternalProject_Add(IpOpt
          SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
          CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/configure --prefix=<INSTALL_DIR> #--libdir=${CMAKE_BINARY_DIR}/lib
          URL http://www.coin-or.org/download/source/Ipopt/Ipopt-3.11.9.tgz
          PATCH_COMMAND      patch -t -s -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/external/ipopt.patch -d ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/get.Mumps
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Metis/get.Metis
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Blas/get.Blas
          COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib
          COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib64
          BUILD_COMMAND ${MAKE} 
          UPDATE_COMMAND     ""
      )

      #add_dependencies(IpOpt mumps metis)
      #add_dependencies(mumps metis)
      #get_property(IpOpt_location TARGET IpOpt PROPERTY LOCATION)
      #add_dependencies(${PROJECT_NAME} IpOpt)

      set(IpOpt_LIB_IpOpt "${CMAKE_BINARY_DIR}/lib/libipopt.so" CACHE INTERNAL "IpOpt shared library")
      set(IpOpt_LIB_Mumps "${CMAKE_BINARY_DIR}/lib/libcoinmumps.so" CACHE INTERNAL "IpOpt shared mumps library")
      set(IpOpt_LIB_Metis "${CMAKE_BINARY_DIR}/lib/libcoinmetis.so" CACHE INTERNAL "IpOpt shared metis library")

      set(IpOpt_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin
          CACHE INTERNAL "IpOpt include path")
    
  endif(IpOpt_LIB)
  include_directories("${IpOpt_INCLUDE_DIRS}")
endif(GISMO_WITH_IPOPT)

if(GISMO_WITH_SUPERLU)
  find_package(SuperLU REQUIRED)
endif(GISMO_WITH_SUPERLU)

if(GISMO_WITH_PARDISO)
   message(WARNING "PARDISO support is not cross-platform yet, more twwaking required.")
   #find_package(PARDISO REQUIRED)
   find_package(MKL REQUIRED)
endif(GISMO_WITH_PARDISO)

## #################################################################
## Build gismo library instantized for GISMO_COEFF_TYPE number type
## #################################################################

if(GISMO_BUILD_LIB)

  if (GISMO_EXTRA_DEBUG)
    set(gismo_SOURCES ${gismo_SOURCES} ${PROJECT_SOURCE_DIR}/src/misc/gsStackWalker.cpp)
  endif()

endif()

include(gsLibrary)

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
## Testing
## #################################################################

include(CTest)
enable_testing()

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


## #################################################################
## Documentation
## #################################################################

add_subdirectory(doc)
