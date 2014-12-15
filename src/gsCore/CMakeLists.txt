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

if (NOT GISMO_EXTRA_DEBUG)
  list(REMOVE_ITEM ${GM_NAME}_CPP ${CMAKE_CURRENT_SOURCE_DIR}/gsStackWalker.cpp)
endif()

## #################################################################
## Add module
## #################################################################

if(GISMO_BUILD_LIB)

  #include_directories(${GISMO_INCLUDE_DIRS})

  add_library(${GM_NAME} OBJECT
  ${${GM_NAME}_H}
  ${${GM_NAME}_HPP}
  #${${GM_NAME}_CPP}
  )

  set_target_properties(${GM_NAME} PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  EXCLUDE_FROM_ALL 1
  EXCLUDE_FROM_DEFAULT_BUILD 1)

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}>
    CACHE INTERNAL "G+Smo modules" )
set (gismo_SOURCES ${gismo_SOURCES} ${${GM_NAME}_CPP}
     CACHE INTERNAL "gismo source files")

endif()
