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

  aux_header_directory(${gismo_externals}/XBraid/braid ${PROJECT_NAME}_HEADERS_XBRAID)
  aux_source_directory(${gismo_externals}/XBraid/braid ${PROJECT_NAME}_SOURCES_XBRAID)
  aux_tmpl_header_directory(${gismo_externals} ${PROJECT_NAME}_HPPFILES_XBRAID)
  
  include_directories(${gismo_externals}/XBraid/braid)

endif (NOT XBRAID_FOUND)

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_HEADERS}
  ${${PROJECT_NAME}_HPPFILES}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS_XBRAID}
  ${${PROJECT_NAME}_HPPFILES_XBRAID}
  ${${PROJECT_NAME}_SOURCES_XBRAID}
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
