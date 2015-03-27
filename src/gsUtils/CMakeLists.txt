######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

get_filename_component(GM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#project(${GM_NAME})

## #################################################################
## Collect files
## #################################################################

aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_H  )
aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_HPP)
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/gsMesh ${GM_NAME}_MESH_H  )
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/gsMesh ${GM_NAME}_MESH_HPP)
aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR}/gsMesh ${GM_NAME}_MESH_CPP)

## #################################################################
## Add module
## #################################################################

if(NOT GISMO_BUILD_LIB)
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_INS)
  LIST( REMOVE_ITEM ${GM_NAME}_CPP ${${GM_NAME}_INS})
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR}/gsMesh ${GM_NAME}_MINS)
  LIST( REMOVE_ITEM ${GM_NAME}_MESH_CPP ${${GM_NAME}_MINS})
endif()

  include_directories(${GISMO_INCLUDE_DIRS})

  add_library(${GM_NAME} OBJECT
  ${${GM_NAME}_H}
  ${${GM_NAME}_HPP}
  ${${GM_NAME}_CPP}
  ${${GM_NAME}_MESH_H}
  ${${GM_NAME}_MESH_HPP}
  ${${GM_NAME}_MESH_CPP}
  )

  set_target_properties(${GM_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS 
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules"
  #EXCLUDE_FROM_ALL 1
  #EXCLUDE_FROM_DEFAULT_BUILD 1
  )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}>
    CACHE INTERNAL "G+Smo modules" )
