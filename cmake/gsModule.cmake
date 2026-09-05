######################################################################
## gsModule.cmake
## This file is part of the G+Smo library.
##
## In-tree module declaration.
##
## Every module directory src/<module>/ has a CMakeLists.txt calling
##
##   gismo_add_module(<module>
##                    [DEPENDS <module>...]   # public dependencies
##                    [SUBDIRS <dir>...]      # nested source directories
##                    [NO_COMPONENT])
##
## which
##
##  1. compiles the module's sources (.h/.hpp/.cpp/_.cpp and the
##     GISMO_EXTRA_INSTANCE instantiations) into the OBJECT library
##     <module> and appends the objects to gismo_MODULES. libgismo is
##     still linked from all module objects in gsLibrary.cmake; building
##     one library per module is a later step that keeps this
##     declaration unchanged;
##
##  2. records the declared dependencies. gismo_check_module_graph()
##     (called from src/CMakeLists.txt) rejects unknown modules and
##     cycles and writes the graph to <build>/gismo_module_graph.txt;
##     .github/scripts/check_includes.py checks that the public headers
##     (.h) of a module only include modules it declares, transitively.
##     Known violations are listed in cmake/dag_undeclared.txt, a list
##     that may only shrink;
##
##  3. creates the component target gismo::<module>: an interface
##     target that links libgismo and the components of the declared
##     dependencies and carries the include directories, so that
##
##       find_package(gismo COMPONENTS gsNurbs)
##       target_link_libraries(app gismo::gsNurbs)
##
##     works now and keeps working when the modules become separate
##     libraries. gismo::gismo is the umbrella target for everything.
##
## In header-only builds (GISMO_BUILD_LIB=OFF) the components link
## gismo_static; executables built that way also compile the
## non-template sources, see add_gismo_pure_executable in gismoUse.cmake.
##
## C interface: a module keeps its C functions in src/<module>/cinterface/
## (gsC<Class>.h/.cpp and the umbrella Cgismo<Module>.h). With
## GISMO_BUILD_CINTERFACE=ON they are compiled into the module and the
## umbrellas are collected into the generated <gsCInterface/Cgismo.h>
## (gismo_generate_cinterface_header, called from the root CMakeLists
## after the optional modules).
######################################################################

function(gismo_add_module NAME)
  cmake_parse_arguments(GM "NO_COMPONENT" "" "DEPENDS;SUBDIRS" ${ARGN})
  if(GM_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "gismo_add_module(${NAME}): unexpected arguments: ${GM_UNPARSED_ARGUMENTS}")
  endif()
  get_property(_modules GLOBAL PROPERTY GISMO_MODULES)
  if(NAME IN_LIST _modules)
    message(FATAL_ERROR "gismo_add_module(${NAME}): module declared twice")
  endif()
  set_property(GLOBAL APPEND PROPERTY GISMO_MODULES ${NAME})
  set_property(GLOBAL PROPERTY GISMO_MODULE_DEPENDS_${NAME} "${GM_DEPENDS}")
  set_property(GLOBAL PROPERTY GISMO_MODULE_DIR_${NAME} "${CMAKE_CURRENT_SOURCE_DIR}")

  ## Collect files
  set(_dirs ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(_sd ${GM_SUBDIRS})
    list(APPEND _dirs ${CMAKE_CURRENT_SOURCE_DIR}/${_sd})
  endforeach()
  set(${NAME}_H)
  set(${NAME}_HPP)
  set(${NAME}_CPP)
  set(${NAME}_INS)
  foreach(_d ${_dirs})
    aux_header_directory     (${_d} _h)
    aux_tmpl_header_directory(${_d} _hpp)
    aux_cpp_noins_directory  (${_d} _cpp)
    list(APPEND ${NAME}_H   ${_h})
    list(APPEND ${NAME}_HPP ${_hpp})
    list(APPEND ${NAME}_CPP ${_cpp})
    if(GISMO_BUILD_LIB)
      aux_instance_directory (${_d} _ins)
      list(APPEND ${NAME}_INS ${_ins})
    endif()
  endforeach()

  ## C interface of the module (primary instance only)
  if(GISMO_BUILD_CINTERFACE AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cinterface)
    aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/cinterface _ch)
    aux_cpp_directory   (${CMAKE_CURRENT_SOURCE_DIR}/cinterface _ccpp)
    list(APPEND ${NAME}_H   ${_ch})
    list(APPEND ${NAME}_CPP ${_ccpp})
    file(GLOB _umbrellas RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/..
      ${CMAKE_CURRENT_SOURCE_DIR}/cinterface/Cgismo*.h)
    set_property(GLOBAL APPEND PROPERTY GISMO_CINTERFACE_UMBRELLAS ${_umbrellas})
  endif()

  ## Compile the module (linked into libgismo by gsLibrary.cmake)
  add_library(${NAME} OBJECT
    ${${NAME}_H}
    ${${NAME}_HPP}
    ${${NAME}_CPP} # static/non templated part
    ${${NAME}_INS}
    )

  set_target_properties(${NAME} PROPERTIES
    # GISMO_PRIMARY_INSTANCE marks the real_t=GISMO_COEFF_TYPE build of a
    # module (extra scalar instances below do not get it); used to emit
    # process-unique symbols such as the XML-registration anchors
    COMPILE_DEFINITIONS "gismo_EXPORTS;GISMO_PRIMARY_INSTANCE"
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    FOLDER "G+Smo modules"
    )

  if(TARGET gmp)
    add_dependencies(${NAME} gmp)
  endif()
  if(TARGET mpfr)
    add_dependencies(${NAME} mpfr)
  endif()

  ## Add extra instances (template instantiations)
  # NOTE: The static/non templated part of the library is
  # compiled above only for GISMO_COEFF_TYPE
  math(EXPR ii "1")
  foreach(ins ${GISMO_EXTRA_INSTANCE})
    add_library(${NAME}_${ii} OBJECT
      ${${NAME}_H}
      ${${NAME}_HPP}
      ${${NAME}_INS}
      )
    set_target_properties(${NAME}_${ii} PROPERTIES
      COMPILE_DEFINITIONS "gismo_EXPORTS;real_t=${ins}"
      POSITION_INDEPENDENT_CODE ON
      LINKER_LANGUAGE CXX
      FOLDER "G+Smo modules"
      )
    set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${NAME}_${ii}>
      CACHE INTERNAL "G+Smo modules" )
    math(EXPR ii "${ii}+1")
  endforeach()

  if (GISMO_BUILD_PCH)
    target_precompiled_header(${NAME} gsPrecompiledHeader)
  endif()

  set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${NAME}>
    CACHE INTERNAL "G+Smo modules" )

  ## Install the public headers
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DESTINATION include/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.h" )

  ## Component target gismo::<module>
  if(NOT GM_NO_COMPONENT)
    add_library(gismo_${NAME} INTERFACE)
    add_library(gismo::${NAME} ALIAS gismo_${NAME})
    set_target_properties(gismo_${NAME} PROPERTIES EXPORT_NAME ${NAME})
    # the module's code is in libgismo (one library for all modules)
    if(GISMO_BUILD_LIB)
      target_link_libraries(gismo_${NAME} INTERFACE gismo)
    else()
      target_link_libraries(gismo_${NAME} INTERFACE gismo_static)
    endif()
    foreach(_dep ${GM_DEPENDS})
      target_link_libraries(gismo_${NAME} INTERFACE gismo_${_dep})
    endforeach()
    install(TARGETS gismo_${NAME} EXPORT gismoTargets)
  endif()
endfunction()

## Checks the declared module graph: every dependency must be a declared
## module, the graph must be acyclic. Sets GISMO_MODULE_ORDER (modules
## in dependency order) and writes <build>/gismo_module_graph.txt.
function(gismo_check_module_graph)
  get_property(_modules GLOBAL PROPERTY GISMO_MODULES)
  foreach(_m ${_modules})
    get_property(_deps GLOBAL PROPERTY GISMO_MODULE_DEPENDS_${_m})
    foreach(_d ${_deps})
      if(_d STREQUAL _m)
        message(FATAL_ERROR "module ${_m} depends on itself")
      endif()
      if(NOT _d IN_LIST _modules)
        message(FATAL_ERROR "module ${_m} depends on ${_d}, which is not a declared module. Declared: ${_modules}")
      endif()
    endforeach()
  endforeach()

  # topological order; if no module can be removed, there is a cycle
  set(_todo ${_modules})
  set(_order)
  while(_todo)
    set(_progress FALSE)
    foreach(_m ${_todo})
      get_property(_deps GLOBAL PROPERTY GISMO_MODULE_DEPENDS_${_m})
      set(_ready TRUE)
      foreach(_d ${_deps})
        if(_d IN_LIST _todo)
          set(_ready FALSE)
        endif()
      endforeach()
      if(_ready)
        list(APPEND _order ${_m})
        list(REMOVE_ITEM _todo ${_m})
        set(_progress TRUE)
      endif()
    endforeach()
    if(NOT _progress)
      message(FATAL_ERROR "cycle in the declared module graph among: ${_todo}")
    endif()
  endwhile()
  set(GISMO_MODULE_ORDER ${_order} CACHE INTERNAL "G+Smo modules in dependency order")

  set(_text "# Declared G+Smo module graph (written by gismo_check_module_graph)\n")
  foreach(_m ${_order})
    get_property(_deps GLOBAL PROPERTY GISMO_MODULE_DEPENDS_${_m})
    foreach(_d ${_deps})
      string(APPEND _text "${_m} -> ${_d}\n")
    endforeach()
  endforeach()
  file(WRITE "${CMAKE_BINARY_DIR}/gismo_module_graph.txt" "${_text}")
  message(STATUS "G+Smo modules (dependency order): ${_order}")
endfunction()

## Generates <gsCInterface/Cgismo.h> in the build tree from the umbrella
## headers of all modules (and of optional modules with a src/cinterface
## directory). Called from the root CMakeLists after the optional modules
## are known.
function(gismo_generate_cinterface_header)
  get_property(_umbrellas GLOBAL PROPERTY GISMO_CINTERFACE_UMBRELLAS)
  set(gsModule_includes "")
  foreach(_u ${_umbrellas})
    string(APPEND gsModule_includes "#include <${_u}>\n")
  endforeach()
  set(gsOptional_includes "")
  foreach(gssm ${GISMO_OPTIONAL})
    string(STRIP ${gssm} SUBMODULE)
    if (EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/src/cinterface")
      string(APPEND gsOptional_includes "#ifdef ${SUBMODULE}_ENABLED\n")
      file(GLOB SUBMODULE_HEADERS ${gismo_SOURCE_DIR}/optional/${SUBMODULE}/src/cinterface/*.h)
      foreach(HEADER ${SUBMODULE_HEADERS})
        get_filename_component(HEADER_NAME ${HEADER} NAME)
        string(APPEND gsOptional_includes "#include <${SUBMODULE}/cinterface/${HEADER_NAME}>\n")
      endforeach()
      string(APPEND gsOptional_includes "#endif\n")
    endif()
  endforeach()
  # generated into the build tree (never the source tree); resolves as
  # <gsCInterface/Cgismo.h> because ${gismo_BINARY_DIR} is an include dir
  configure_file("${gismo_SOURCE_DIR}/src/gsCInterface/Cgismo.h.in"
                 "${gismo_BINARY_DIR}/gsCInterface/Cgismo.h" @ONLY)
  install(FILES "${gismo_BINARY_DIR}/gsCInterface/Cgismo.h"
          DESTINATION include/${PROJECT_NAME}/gsCInterface)
  list(LENGTH _umbrellas _n)
  message(STATUS "C interface: Cgismo.h generated from ${_n} module umbrellas")
endfunction()
