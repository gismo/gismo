######################################################################
## CMakeLists.txt --- src
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

include_directories(${GISMO_INCLUDE_DIRS})

if (GISMO_BUILD_PCH)
  add_precompiled_header(gsPrecompiledHeader misc/gsPrecompiledHeader.h)
endif()

include(gsModule)

## Modules. Each src/<module>/CMakeLists.txt declares the module and its
## dependencies with gismo_add_module(); the order below is the
## dependency order of the declared graph (checked at the end).
add_subdirectory(gsCInterface) # type-blind core of the C interface
add_subdirectory(gsParallel)
add_subdirectory(gsMatrix)
add_subdirectory(gsUtils)
add_subdirectory(gsCore)
add_subdirectory(gsIO)
add_subdirectory(gsDomain)
add_subdirectory(gsTensor)
add_subdirectory(gsFunctionExpr)
add_subdirectory(gsContainers)
add_subdirectory(gsNurbs)
add_subdirectory(gsMesh2)
add_subdirectory(gsMSplines)
add_subdirectory(gsHSplines)
add_subdirectory(gsPde)
add_subdirectory(gsExpressions)
add_subdirectory(gsAssembler)
add_subdirectory(gsOptimizer)
add_subdirectory(gsSolver)
add_subdirectory(gsModeling)
add_subdirectory(gsMultiGrid)
add_subdirectory(gsIeti)
add_subdirectory(gsModules)

## Every gs* directory of src/ must be a declared module
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
get_property(GISMO_DECLARED_MODULES GLOBAL PROPERTY GISMO_MODULES)
foreach(subdir ${SUBDIRS})
  if(NOT subdir IN_LIST GISMO_DECLARED_MODULES)
    message(FATAL_ERROR "src/${subdir} is not a declared module: add src/${subdir}/CMakeLists.txt calling gismo_add_module(${subdir} ...) and an add_subdirectory(${subdir}) line in src/CMakeLists.txt")
  endif()
endforeach()

gismo_check_module_graph()
