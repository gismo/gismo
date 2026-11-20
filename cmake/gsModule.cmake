######################################################################
## gsModule.cmake
## This file is part of the G+Smo library.
##
## Functions to build optional modules either as dynamically loaded
## libraries or into the G+Smo library (built-in).
##
## Author: Angelos Mantzaflaris
######################################################################

include(CMakeParseArguments)

# Directory of the dynamic modules in the build tree and after installation
set(GISMO_MODULE_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/gismo/modules")      # = <build>/lib/gismo/modules
if(IS_ABSOLUTE "${LIB_INSTALL_DIR}")
  set(GISMO_MODULE_INSTALL_DIR "${LIB_INSTALL_DIR}/gismo/modules")
else()
  set(GISMO_MODULE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}/gismo/modules")
endif()

if(GISMO_DYNAMIC_MODULES AND NOT GISMO_BUILD_LIB)
  message(WARNING "GISMO_DYNAMIC_MODULES requires GISMO_BUILD_LIB=ON; building ${GISMO_DYNAMIC_MODULES} into the library")
  set(GISMO_DYNAMIC_MODULES "")       # plain variable shadows the cache entry
endif()

# The module names end up in C symbols, hence they must be C identifiers
function(gismo_check_module_name name)
  if(NOT "${name}" MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "Module name '${name}' is not a valid C identifier")
  endif()
endfunction()

# gismo_add_module_library(<name> SOURCES <files> [LINK <libs>])
#
# Builds the optional module <name> as a dynamically loaded library,
# placed in GISMO_MODULE_DIR and installed to GISMO_MODULE_INSTALL_DIR.
# Requires GISMO_BUILD_LIB.
function(gismo_add_module_library name)
  gismo_check_module_name(${name})
  cmake_parse_arguments(GML "" "" "SOURCES;LINK" ${ARGN})
  if(NOT GISMO_BUILD_LIB)
    message(FATAL_ERROR "gismo_add_module_library(${name}) requires GISMO_BUILD_LIB=ON")
  endif()
  if(NOT GML_SOURCES)
    message(FATAL_ERROR "gismo_add_module_library(${name}) needs SOURCES")
  endif()

  add_library(${name} MODULE ${GML_SOURCES})
  # The literal target name is used, since PROJECT_NAME is the module's here
  target_link_libraries(${name} gismo ${GML_LINK})
  target_compile_definitions(${name} PRIVATE GISMO_MODULE_DYNAMIC)

  set_target_properties(${name} PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${GISMO_MODULE_DIR}"
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET "hidden"
    FOLDER "G+Smo modules")
  foreach(config ${CMAKE_CONFIGURATION_TYPES}) # For Visual studio
    string(TOUPPER ${config} CONFIG)
    set_target_properties(${name} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${GISMO_MODULE_DIR}")
  endforeach()
  # Relative path from lib/gismo/modules to lib
  if(UNIX AND NOT APPLE)
    set_target_properties(${name} PROPERTIES INSTALL_RPATH "\$ORIGIN/../..")
  elseif(APPLE)
    set_target_properties(${name} PROPERTIES INSTALL_RPATH "@loader_path/../..")
  endif()

  install(TARGETS ${name}
    LIBRARY DESTINATION "${GISMO_MODULE_INSTALL_DIR}" COMPONENT modules)

  set_property(GLOBAL APPEND PROPERTY GISMO_DYNAMIC_MODULE_LIST ${name})
endfunction()

# gismo_add_builtin_module(<name>)
#
# Registers the optional module <name>, whose object code is part of
# the G+Smo library, in the list of built-in modules.
function(gismo_add_builtin_module name)
  gismo_check_module_name(${name})
  set_property(GLOBAL APPEND PROPERTY GISMO_BUILTIN_MODULE_LIST ${name})
endfunction()

# gismo_configure_modules()
#
# Generates gsModuleConfig.h and gsBuiltinModules.cpp in the build
# tree. To be called after all optional modules have been added.
function(gismo_configure_modules)
  get_property(_builtin GLOBAL PROPERTY GISMO_BUILTIN_MODULE_LIST)
  if(_builtin)
    list(REMOVE_DUPLICATES _builtin)
  endif()
  get_property(_dynamic GLOBAL PROPERTY GISMO_DYNAMIC_MODULE_LIST)
  if(_dynamic)
    list(REMOVE_DUPLICATES _dynamic)
  endif()

  foreach(m ${GISMO_DYNAMIC_MODULES})
    list(FIND _dynamic ${m} _idx)
    if(_idx EQUAL -1)
      message(WARNING "${m} is listed in GISMO_DYNAMIC_MODULES but was not built as a dynamic module")
    endif()
  endforeach()

  set(GISMO_BUILTIN_MODULE_DECLS "")
  set(GISMO_BUILTIN_MODULE_NAMES "")
  set(GISMO_BUILTIN_MODULE_CALLS "")
  foreach(m ${_builtin})
    set(GISMO_BUILTIN_MODULE_DECLS "${GISMO_BUILTIN_MODULE_DECLS}extern \"C\" void gismo_module_register_${m}();\n")
    set(GISMO_BUILTIN_MODULE_NAMES "${GISMO_BUILTIN_MODULE_NAMES}\"${m}\", ")
    set(GISMO_BUILTIN_MODULE_CALLS "${GISMO_BUILTIN_MODULE_CALLS}    if (0 == std::strcmp(name, \"${m}\")) { gismo_module_register_${m}(); return true; }\n")
  endforeach()

  configure_file("${gismo_SOURCE_DIR}/src/gsSystem/gsModuleConfig.h.in"
    "${gismo_BINARY_DIR}/gsSystem/gsModuleConfig.h" @ONLY)
  configure_file("${gismo_SOURCE_DIR}/src/gsSystem/gsBuiltinModules.cpp.in"
    "${gismo_BINARY_DIR}/gsSystem/gsBuiltinModules.cpp" @ONLY)

  # Installed modules include gsModule.h, which includes this header
  install(FILES "${gismo_BINARY_DIR}/gsSystem/gsModuleConfig.h"
    DESTINATION include/gismo/gsSystem/)

  message("  GISMO_DYNAMIC_MODULES ${_dynamic}")
  message("  Built-in modules      ${_builtin}")
endfunction()
