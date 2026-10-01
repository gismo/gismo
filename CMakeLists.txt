
#####################################################################
## CMakeLists.txt --- root
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

if(CMAKE_VERSION VERSION_LESS "3.19")
  cmake_minimum_required(VERSION 2.8.12)
else()
  cmake_minimum_required(VERSION 3.1...3.10)
endif()

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
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()

set(CMAKE_BUILD_TYPE_INIT Release)
if(NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;MinSizeRel" CACHE STRING "" FORCE)
endif()
if(NOT CMAKE_BUILD_TYPE)
   # Set default build type to Release
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "RelWithDebInfo" "MinSizeRel")
endif()

set(gismo_VERSION_MAJOR 26) #year
set(gismo_VERSION_MINOR 09) #month
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

include(submodules.txt)
set(GISMO_OPTIONAL ${SUBMODULES_TXT} ${GISMO_OPTIONAL} CACHE STRING "List of optional G+Smo modules" FORCE)
if (CMAKE_VERSION VERSION_GREATER 3.11.5)
  list(TRANSFORM GISMO_OPTIONAL STRIP)
endif()
set(GISMO_OPTIONAL ${GISMO_OPTIONAL}
  CACHE STRING "List of optional G+Smo modules" FORCE)
#message(STATUS "Submodules ${GISMO_OPTIONAL}")

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
unset(GISMO_EXTERNALS            CACHE)
unset(${PROJECT_NAME}_HEADERS    CACHE)
unset(${PROJECT_NAME}_SOURCES    CACHE)
unset(${PROJECT_NAME}_EXTENSIONS CACHE)
unset(${PROJECT_NAME}_MODULES    CACHE)
unset(${PROJECT_NAME}_optionals  CACHE)
unset(${PROJECT_NAME}_LINKER     CACHE)
unset(gismo_TMPL_HEADERS         CACHE)
unset(gismo_UNITTESTS            CACHE)
set(${PROJECT_NAME}_optionals ""
    CACHE INTERNAL "${PROJECT_NAME} list of optionals")

set (GISMO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/external"
     "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/optional"
     CACHE INTERNAL "${PROJECT_NAME} include directories" )

# Offer the user the choice of overriding the installation directories
set(LIB_INSTALL_DIR     lib     CACHE STRING "Installation directory for libraries")
set(BIN_INSTALL_DIR     bin     CACHE STRING "Installation directory for executables")
set(INCLUDE_INSTALL_DIR include CACHE STRING "Installation directory for header files")

## #################################################################
## Extensions
## #################################################################

#functions to fetch(download) head revision from git or Subversion
include(gsFetch)

# Eigen 5 requires C++14. "98" needs its own test: it compares numerically greater than 14.
if("x${CMAKE_CXX_STANDARD}" STREQUAL "x98" OR CMAKE_CXX_STANDARD LESS 14)
  message(FATAL_ERROR "Eigen 5 requires C++14 or newer (CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD})")
endif()

if(GISMO_Eigen_MIRROR)
  # Without a usable clone, reset the fetch: its ExternalProject stamps in the build tree survive a deleted
  # clone (the re-download is then skipped), and that skipped fetch leaves an empty external/Eigen behind,
  # which gismo_fetch_directory() treats as present. A non-empty folder is never deleted here.
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/Eigen/signature_of_eigen3_matrix_library")
    file(REMOVE_RECURSE "${gismo_BINARY_DIR}/external/Eigen_fetch")
    file(GLOB _gismo_eigen_clone "${PROJECT_SOURCE_DIR}/external/Eigen/*")
    if(NOT _gismo_eigen_clone)
      file(REMOVE_RECURSE "${PROJECT_SOURCE_DIR}/external/Eigen")
    endif()
  endif()
  gismo_fetch_directory(Eigen
    GIT_REPOSITORY ${GISMO_Eigen_MIRROR_URL}
    GIT_TAG        5.0.0
    GIT_SHALLOW    TRUE
    DESTINATION    external)
  set(Eigen_DIR "${PROJECT_SOURCE_DIR}/external/Eigen")   # normal variable: shadows a user cache entry, leaves it intact
endif()

find_package(Eigen 5.0)                                     # no REQUIRED: our own FATAL_ERROR below names the remedies

if(NOT Eigen_FOUND OR NOT Eigen_VERSION VERSION_LESS 6)
  if(GISMO_Eigen_MIRROR)
    message(FATAL_ERROR "G+Smo requires Eigen 5.x, but ${PROJECT_SOURCE_DIR}/external/Eigen (cloned from ${GISMO_Eigen_MIRROR_URL}, tag 5.0.0) does not contain a usable Eigen 5.x. gismo_fetch_directory() never re-fetches an existing directory, so if that clone is missing, partial, or the wrong version, delete ${PROJECT_SOURCE_DIR}/external/Eigen by hand and reconfigure.")
  else()
    # cmake/FindEigen.cmake sets Eigen_REJECTED_DIR (and unsets EIGEN_INCLUDE_DIR)
    # only when the find failed outright. A too-new Eigen (found, but rejected by
    # the VERSION_LESS 6 test above) leaves EIGEN_INCLUDE_DIR cached and
    # Eigen_REJECTED_DIR empty, so fall back to it for the message. With Eigen_DIR
    # searched exclusively (NO_DEFAULT_PATH), an empty/NOTFOUND rejected dir means
    # nothing was found in Eigen_DIR itself, which the message should name.
    set(_gismo_eigen_rejected_dir "${Eigen_REJECTED_DIR}")
    if(NOT _gismo_eigen_rejected_dir)
      set(_gismo_eigen_rejected_dir "${EIGEN_INCLUDE_DIR}")
    endif()
    set(_gismo_eigen_rejected_version "")
    if(Eigen_VERSION)
      set(_gismo_eigen_rejected_version " (version '${Eigen_VERSION}')")
    endif()
    if((NOT _gismo_eigen_rejected_dir OR _gismo_eigen_rejected_dir MATCHES "-NOTFOUND$") AND Eigen_DIR)
      set(_gismo_eigen_found_desc "no Eigen found in Eigen_DIR='${Eigen_DIR}'")
    else()
      set(_gismo_eigen_found_desc "found: '${_gismo_eigen_rejected_dir}'${_gismo_eigen_rejected_version}")
    endif()
    message(FATAL_ERROR "G+Smo requires Eigen 5.x (${_gismo_eigen_found_desc}). Either pass -DEigen_DIR=<Eigen 5 source tree or install prefix>, or configure with -DGISMO_Eigen_MIRROR=ON to clone Eigen 5.0.0 into external/Eigen.")
  endif()
endif()
message(STATUS "Eigen ${Eigen_VERSION}: ${EIGEN_INCLUDE_DIR}")

set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories" FORCE)

if (EIGEN_USE_MKL_ALL)
  # See http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html
  find_package(MKL REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MKL_INCLUDE_DIR} ${MKL_INCLUDE}
    CACHE INTERNAL "${PROJECT_NAME} include directories")
endif()

if (";${GISMO_OPTIONAL};" MATCHES ";gsUniversal;")
  gismo_fetch_module(gsUniversal)
endif()

if (";${GISMO_OPTIONAL};" MATCHES ";gsCoDiPack;")
  gismo_fetch_module(gsCoDiPack)
endif()

# Ensure AutoDiff includes are available early if the module exists
if(EXISTS "${CMAKE_SOURCE_DIR}/optional/gsAutoDiff" AND autodiff_FOUND)
  get_target_property(AUTODIFF_INCLUDE_DIRS autodiff::autodiff INTERFACE_INCLUDE_DIRECTORIES)
  if(AUTODIFF_INCLUDE_DIRS)
    set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${AUTODIFF_INCLUDE_DIRS}
      CACHE INTERNAL "Gismo include directories" FORCE)
  endif()
endif()

if (";${GISMO_OPTIONAL};" MATCHES ";gsAutoDiff;")
  string(STRIP "gsAutoDiff" SUBMODULE)
  if(EXISTS "${PROJECT_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
    message(STATUS "Processing dependencies for optional module: ${SUBMODULE}")
    include("${PROJECT_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
  endif()
  gismo_fetch_module(gsAutoDiff)
endif()

#first time
include_directories(${GISMO_INCLUDE_DIRS})

# Add externals directory
add_subdirectory(external)

if(GISMO_WITH_MPI)
  find_package(MPI REQUIRED)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPI_INCLUDE_PATH}
  CACHE INTERNAL "${PROJECT_NAME} include directories)" FORCE)
  set(gismo_LINKER ${gismo_LINKER} ${MPI_CXX_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects" FORCE)
endif(GISMO_WITH_MPI)

if(${GISMO_COEFF_TYPE} STREQUAL "mpq_class")
  include(external/gsGmp.cmake)
endif()

if(${GISMO_COEFF_TYPE} STREQUAL "mpfr::mpreal")
  include(external/gsGmp.cmake)
  include(external/gsMpfr.cmake)
endif()

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

# Process external dependencies for optional modules BEFORE building sources
list(REMOVE_DUPLICATES GISMO_OPTIONAL)
list(SORT GISMO_OPTIONAL)
foreach(gssm ${GISMO_OPTIONAL})
  string(STRIP ${gssm} SUBMODULE)
  if(EXISTS "${PROJECT_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
    message(STATUS "Processing dependencies for optional module: ${SUBMODULE}")
    include("${PROJECT_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
  endif()
endforeach()

#second time
include_directories(${GISMO_INCLUDE_DIRS})

# external inclusion paths
#get_property(EXTERNAL_INCLUDE_DIRS GLOBAL PROPERTY EXTERNAL_INCLUDE_DIRS)
#include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external)

## #################################################################
## G+Smo sources
## #################################################################

# Add main sources directory
add_subdirectory(src)

## #################################################################
## Submodules
## #################################################################

option(GISMO_UNSUPPORTED "Enable unsupported module" OFF)
if(GISMO_UNSUPPORTED OR GISMO_MOTOR)
  gismo_fetch_module(unsupported)
  set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GISMO_DEV_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME} include directories")
endif()

option(GISMO_MOTOR "Enable motor module" OFF)
if(GISMO_MOTOR)
  #gismo_fetch_module(motor)
  set(GISMO_OPTIONAL ${GISMO_OPTIONAL} "motor"
    CACHE INTERNAL "List of optional G+Smo modules" FORCE)
endif()

if (CMAKE_VERSION VERSION_GREATER 3.11.5)
  list(TRANSFORM GISMO_OPTIONAL STRIP)
endif()
list(REMOVE_DUPLICATES GISMO_OPTIONAL)
list(SORT GISMO_OPTIONAL)
set(GISMO_OPTIONAL ${GISMO_OPTIONAL}
  CACHE STRING "List of optional G+Smo modules" FORCE)
foreach(gssm ${GISMO_OPTIONAL})
  string(STRIP ${gssm} SUBMODULE)
  #message("Fetch ${SUBMODULE}")
  gismo_fetch_module(${SUBMODULE})
  # Prevent special characters to be used in the variable name
  string(REPLACE "-" "_" SUBMODULE_NAME ${SUBMODULE})
  string(REPLACE "+" "_" SUBMODULE_NAME ${SUBMODULE_NAME})
  string(REPLACE "=" "_" SUBMODULE_NAME ${SUBMODULE_NAME})
  string(REPLACE "." "_" SUBMODULE_NAME ${SUBMODULE_NAME})
  string(REPLACE "!" "_" SUBMODULE_NAME ${SUBMODULE_NAME})
  string(REPLACE "*" "_" SUBMODULE_NAME ${SUBMODULE_NAME})
  set(CONFIGEXT_H "${CONFIGEXT_H}#define ${SUBMODULE_NAME}_ENABLED\n")
endforeach()

foreach(gssm ${GISMO_EXTERNALS})
  string(STRIP ${gssm} EXTMOD)
  #message("Fetch ext: ${EXTMOD}")
  set(CONFIGEXT_H "${CONFIGEXT_H}#define ${EXTMOD}_ENABLED\n")
endforeach()
file(WRITE ${gismo_BINARY_DIR}/gsCore/gsConfigExt.txt "${CONFIGEXT_H}")
execute_process( COMMAND ${CMAKE_COMMAND} -E compare_files gsConfigExt.txt gsConfigExt.h
  RESULT_VARIABLE comp_res
  WORKING_DIRECTORY ${gismo_BINARY_DIR}/gsCore)
if(NOT comp_res EQUAL 0)
  file(WRITE ${gismo_BINARY_DIR}/gsCore/gsConfigExt.h "${CONFIGEXT_H}")
endif()

message("  GISMO_OPTIONAL        ${GISMO_OPTIONAL}")
set(GISMO_OPTIONAL ${GISMO_OPTIONAL}
  CACHE STRING "List of optional G+Smo modules" FORCE)

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
## Python bindings with cppyy
## #################################################################

if(GISMO_WITH_CPPYY)
  add_subdirectory(cppyy)
endif()


## #################################################################
## Plugins
## #################################################################

add_custom_target(plugins)

if(GISMO_PLUGIN_AXL)
  add_subdirectory(plugins/gsAxel)
endif(GISMO_PLUGIN_AXL)

if(GISMO_BUILD_MEX)
  add_subdirectory(plugins/gsMex)
endif(GISMO_BUILD_MEX)

#if(GISMO_PLUGIN_PVIEW)
#  add_subdirectory(plugins/gsParaview)
#endif(GISMO_PLUGIN_PVIEW)

if(GISMO_PLUGIN_RHINO)
  add_subdirectory(plugins/RhinoTHBSplines)
endif(GISMO_PLUGIN_RHINO)

## #################################################################
## Unittests
## #################################################################

if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(unittests)
else()
  add_custom_target(unittests COMMAND "" COMMENT "Set CMake argument GISMO_BUILD_UNITTESTS=ON to enable unittests")
endif(GISMO_BUILD_UNITTESTS)

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

unset(GISMO_UNAME CACHE)
unset(GISMO_PASS CACHE)
unset(RESTORE_SUBMODULES CACHE)

#message("GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS}")
