######################################################################
## CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

get_filename_component(GM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${GM_NAME})

## #################################################################
## Collect files
## #################################################################

aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_H  )
aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPP)

## #################################################################
## Add module
## #################################################################

include_directories(${GISMO_INCLUDE_DIRS})

if(GISMO_BUILD_LIB)

  aux_generate_instance( ${CMAKE_CURRENT_SOURCE_DIR})
  aux_instance_directory   (${PROJECT_BINARY_DIR}/instance ${PROJECT_NAME}_INS)

  add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H} 
  ${${PROJECT_NAME}_HPP} 
  ${${PROJECT_NAME}_CPP}
  ${${PROJECT_NAME}_INS}
  )

  set_target_properties(${PROJECT_NAME}
  PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "G+Smo modules" )

endif()
