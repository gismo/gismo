### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## XBraid extension
project(gsXBraidExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

# Look for pre-installed XBraid libraries
find_package(XBRAID QUIET)

if (NOT XBRAID_FOUND)

  set(XBRAID_VER "v3.0.0")  
  gismo_fetch_directory(XBraid
    URL https://github.com/XBraid/xbraid/archive/${XBRAID_VER}.zip
    DESTINATION    external
    )  

  if( (NOT GISMO_BUILD_LIB) )
    aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
    if(${PROJECT_NAME}_INS)
      LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
    endif()
  endif()

  set(${PROJECT_NAME}_HEADERS "${${PROJECT_NAME}_HEADERS}"
    ${gismo_externals}/XBraid/braid/_braid.h
    ${gismo_externals}/XBraid/braid/base.h
    ${gismo_externals}/XBraid/braid/status.h
    ${gismo_externals}/XBraid/braid/tape.h
    ${gismo_externals}/XBraid/braid/util.h
    ${gismo_externals}/XBraid/braid/braid.h
    ${gismo_externals}/XBraid/braid/braid_status.h
    ${gismo_externals}/XBraid/braid/braid_test.h)

  if(NOT GISMO_WITH_MPI )
    set(${PROJECT_NAME}_HEADERS "${${PROJECT_NAME}_HEADERS}"
      ${gismo_externals}/XBraid/braid/mpistubs.h)
    set(${PROJECT_NAME}_SOURCES "${${PROJECT_NAME}_SOURCES}"
      ${gismo_externals}/XBraid/braid/mpistubs.c)
  endif()

  set(${PROJECT_NAME}_SOURCES "${${PROJECT_NAME}_SOURCES}"
    ${gismo_externals}/XBraid/braid/access.c
    ${gismo_externals}/XBraid/braid/adjoint.c
    ${gismo_externals}/XBraid/braid/base.c
    ${gismo_externals}/XBraid/braid/braid.c
    ${gismo_externals}/XBraid/braid/braid_status.c
    ${gismo_externals}/XBraid/braid/braid_test.c
    ${gismo_externals}/XBraid/braid/communication.c
    ${gismo_externals}/XBraid/braid/distribution.c
    ${gismo_externals}/XBraid/braid/drive.c
    ${gismo_externals}/XBraid/braid/grid.c
    ${gismo_externals}/XBraid/braid/hierarchy.c
    ${gismo_externals}/XBraid/braid/interp.c
    ${gismo_externals}/XBraid/braid/norm.c
    ${gismo_externals}/XBraid/braid/refine.c
    ${gismo_externals}/XBraid/braid/relax.c
    ${gismo_externals}/XBraid/braid/residual.c
    ${gismo_externals}/XBraid/braid/restrict.c
    ${gismo_externals}/XBraid/braid/space.c
    ${gismo_externals}/XBraid/braid/step.c
    ${gismo_externals}/XBraid/braid/tape.c
    ${gismo_externals}/XBraid/braid/util.c
    ${gismo_externals}/XBraid/braid/uvector.c)
  
  include_directories(${gismo_externals}/XBraid/braid)

endif (NOT XBRAID_FOUND)

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_HEADERS}
  ${${PROJECT_NAME}_HPPFILES}
  ${${PROJECT_NAME}_SOURCES}
 )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo extensions" )

if( GISMO_WITH_MPI )
  target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})
else()
  add_definitions("-Dbraid_SEQUENTIAL")
endif()

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "Gismo extensions to be included")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
  DESTINATION include/gismo/gsXBraid
  FILES_MATCHING PATTERN "*.h")
