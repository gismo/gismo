######################################################################
## gsModuleLib.cmake
## This file is part of the G+Smo library.
##
## Builds a runtime-loadable gismo module library (.so/.dll).
##
## Usage in a module's CMakeLists.txt:
##
##   if(GISMO_BUILD_LIB AND GISMO_BUILD_MODULE_LIB)
##     include(gsModuleLib)
##     gismo_add_module_lib(gsMyModule "1.0.0"
##       ${CMAKE_CURRENT_SOURCE_DIR}/gsMyModuleEntry.cpp)
##   endif()
##
## The function configures <name>ModuleVersion.h from
## src/gsModules/gsModuleVersion.h.in, builds a MODULE library named
## <name>_module linked against libgismo, and installs it into
## ${GISMO_MODULE_INSTALL_DIR} - the directory that
## gsModuleManager::loadAll() searches at runtime.
######################################################################

function(gismo_add_module_lib MODULE_NAME MODULE_VERSION)
  set(_sources ${ARGN})

  # Per-module version header for the ABI handshake
  set(GISMO_MODULE_NAME    "${MODULE_NAME}")
  set(GISMO_MODULE_VERSION "${MODULE_VERSION}")
  configure_file(
    "${gismo_SOURCE_DIR}/src/gsModules/gsModuleVersion.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}ModuleVersion.h" @ONLY)

  add_library(${MODULE_NAME}_module MODULE ${_sources})

  # Links against the shared gismo library so that GISMO_EXPORT'ed
  # registry singletons resolve to the single copy in libgismo
  target_link_libraries(${MODULE_NAME}_module PRIVATE gismo)

  target_include_directories(${MODULE_NAME}_module PRIVATE
    ${GISMO_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}) # generated <name>ModuleVersion.h

  set_target_properties(${MODULE_NAME}_module PROPERTIES
    # only the two C-linkage entry points need to be visible, but the
    # default-hidden preset would also hide them; keep default and rely
    # on GISMO_EXPORT at the entry points
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    OUTPUT_NAME "${MODULE_NAME}_module"
    PREFIX ""  # no "lib" prefix: the loader skips libgismo*
    FOLDER "G+Smo modules")

  install(TARGETS ${MODULE_NAME}_module
    LIBRARY DESTINATION "${LIB_INSTALL_DIR}/gismo/modules"
    RUNTIME DESTINATION "${LIB_INSTALL_DIR}/gismo/modules" # DLL on Windows
    COMPONENT modules)

  message(STATUS "  gsModuleLib: added runtime module ${MODULE_NAME}_module v${MODULE_VERSION}")
endfunction()
