### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012-2014 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0048 OLD)
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
message( "                  version ${${PROJECT_NAME}_VERSION} (revision ${GISMO_WC_REVISION})" )

## #################################################################
## Project options
## #################################################################

include( gsOptions )

## #################################################################
## Configuration and Compilation
## #################################################################

include( gsConfig )

## #################################################################
## Preprocessor definitions
## #################################################################

if(GISMO_BUILD_SHARED_LIB OR GISMO_BUILD_STATIC_LIB )
set (GISMO_BUILD_LIB ON)
endif(GISMO_BUILD_SHARED_LIB  OR GISMO_BUILD_STATIC_LIB)

#index type
set (GISMO_INDEX_TYPE "int")

# Shared pointe
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

# Note: if both shared and static libs are created, targets are linked against the shared lib.
#       if only one of the two is build, targets are linked using that one
#       if none of the two is build, targets are not linked against anything, only templated code is used.

## #################################################################
## Collect headers and sources
## #################################################################

# Clean previously set variables in the cache
unset(gismo_HEADERS    CACHE)
unset(gismo_SOURCES    CACHE)
unset(gismo_EXTENSIONS CACHE)
unset(gismo_MODULES    CACHE)
unset(gismo_examples   CACHE)
unset(gismo_optionals  CACHE)
#
unset(gismo_TMPL_HEADERS CACHE)

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

include(ExternalProject)

#ExternalProject_Add(metis
#    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/metis"
#    URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
#    PATCH_COMMAND      patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/cmake/metis-4.0-cmake.patch
#    UPDATE_COMMAND     ""
#    INSTALL_COMMAND    ""
#)
#ExternalProject_Get_Property(metis install_dir)
#add_library(Metis UNKNOWN IMPORTED)
#set_property(TARGET Metis PROPERTY IMPORTED_LOCATION ${install_dir}/src/metis-build/Lib/libMetis.a)
#set(METIS_LIBRARIES Metis)
#set(METIS_DIR ${install_dir} CACHE INTERNAL "Metis ROOT dir")
#set(METIS_INCLUDE_DIRS ${install_dir}/include CACHE INTERNAL "Metis include dirs")
#message("Metis install dir ${install_dir}")
#ExternalProject_Add(mumps
#    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/Mumps"
#    URL "http://mumps.enseeiht.fr/MUMPS_4.10.0.tar.gz"
#    PATCH_COMMAND      patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/external/mumps.patch COMMAND patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MUMPS_4.10.0-cmake.patch
#    UPDATE_COMMAND     ""
#    INSTALL_COMMAND    ""
#)

ExternalProject_Add(IpOpt
    SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/configure --prefix=<INSTALL_DIR> #--libdir=${CMAKE_BINARY_DIR}/lib
    URL http://www.coin-or.org/download/source/Ipopt/Ipopt-3.11.9.tgz
    PATCH_COMMAND      patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/external/ipopt.patch -d ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/get.Mumps
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Metis/get.Metis
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib64
    BUILD_COMMAND ${MAKE} 
    UPDATE_COMMAND     ""
)

#add_dependencies(IpOpt mumps metis)
#add_dependencies(mumps metis)
#get_property(IpOpt_location TARGET IpOpt PROPERTY LOCATION)

set(IpOpt_LIBS
    "${CMAKE_BINARY_DIR}/lib/libipopt.so;${CMAKE_BINARY_DIR}/lib/libcoinmumps.so;${CMAKE_BINARY_DIR}/lib/libcoinmetis.so"
    CACHE INTERNAL "IpOpt shared library")

set(IpOpt_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin
    CACHE INTERNAL "IpOpt include path")

include_directories("${IpOpt_INCLUDE_DIRS}")
 
endif(GISMO_WITH_IPOPT)


## #################################################################
## Build gismo library instantized for GISMO_COEFF_TYPE number type
## #################################################################

if(GISMO_BUILD_LIB)

  if (GISMO_EXTRA_DEBUG)
    set(gismo_SOURCES ${gismo_SOURCES} ${PROJECT_SOURCE_DIR}/src/misc/gsStackWalker.cpp)
  endif()

  include(gsLibrary)

endif()

## #################################################################
## Plugins
## #################################################################

if(GISMO_BUILD_AXL)
  add_subdirectory(plugins/gsAxel)
endif(GISMO_BUILD_AXL)

if(GISMO_BUILD_PVIEW)
  add_subdirectory(plugins/gsParaview)
endif(GISMO_BUILD_PVIEW)

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

enable_testing()
include(CTest)


## #################################################################
## Executables
## #################################################################

if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(examples)
  list(LENGTH gismo_examples len)
  message(STATUS "Number of examples to compile:      ${len}")
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
