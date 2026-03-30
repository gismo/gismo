######################################################################
## gsModule.cmake
## This file is part of the G+Smo library.
##
## Helper macro for building gismo MODULE (.so/.dll) shared libraries.
##
## Usage in an optional module's CMakeLists.txt:
##
##   if(GISMO_BUILD_LIB AND GISMO_BUILD_MODULE_LIB)
##     include(gsModule)
##     gismo_add_module(gsMyModule "1.0.0"
##       ${CMAKE_CURRENT_SOURCE_DIR}/gsMyModuleEntry.cpp)
##   endif()
##
## The macro:
##   1. Configures gsModuleVersion.h.in → <ModuleName>ModuleVersion.h
##      in the current binary directory.
##   2. Builds a CMake MODULE (.so) library named <ModuleName>_module.
##   3. Links the MODULE .so against libgismo (shared).
##   4. Installs the MODULE .so into <LIB_INSTALL_DIR>/gismo/modules/.
##
## Author: H.M. Verhelst
######################################################################

macro(gismo_add_module MODULE_NAME MODULE_VERSION)
  set(_gsmodule_sources ${ARGN})

  # --- Generate per-module version header ---
  set(GISMO_MODULE_NAME    "${MODULE_NAME}")
  set(GISMO_MODULE_VERSION "${MODULE_VERSION}")
  configure_file(
    "${gismo_SOURCE_DIR}/src/gsModules/gsModuleVersion.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}ModuleVersion.h"
    @ONLY)

  # --- Build the MODULE shared library ---
  add_library(${MODULE_NAME}_module MODULE
    ${_gsmodule_sources})

  # The MODULE .so links against the main gismo shared library so that
  # symbols like gsOptimizerRegistry::getMutableInstance() resolve to the
  # single authoritative copy in libgismo.so (not a local duplicate).
  target_link_libraries(${MODULE_NAME}_module PRIVATE gismo)

  target_include_directories(${MODULE_NAME}_module PRIVATE
    ${GISMO_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}   # for generated <Module>ModuleVersion.h
  )

  set_target_properties(${MODULE_NAME}_module PROPERTIES
    # Export only the two C-linkage entry points; everything else stays local.
    CXX_VISIBILITY_PRESET   default
    C_VISIBILITY_PRESET     default
    VISIBILITY_INLINES_HIDDEN OFF
    COMPILE_DEFINITIONS     gismo_EXPORTS
    OUTPUT_NAME             "${MODULE_NAME}_module"
    PREFIX                  ""    # avoid "lib" prefix on Linux
    FOLDER                  "G+Smo modules"
  )

  # Install into the directory that gsModuleLoader::loadAll() searches.
  install(TARGETS ${MODULE_NAME}_module
    LIBRARY DESTINATION "${LIB_INSTALL_DIR}/gismo/modules"
    RUNTIME DESTINATION "${LIB_INSTALL_DIR}/gismo/modules"  # DLL on Windows
    COMPONENT modules
  )

  message(STATUS "  gsModule: added MODULE target ${MODULE_NAME}_module v${MODULE_VERSION}")

endmacro()
