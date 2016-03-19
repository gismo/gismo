######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

include_directories(${GISMO_INCLUDE_DIRS})

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
#LIST( REMOVE_ITEM SUBDIRS misc)
LIST( APPEND SUBDIRS gsUtils/gsMesh)

FOREACH(subdir ${SUBDIRS})

  get_filename_component(GM_NAME ${subdir} NAME)
  #message("module: ${GM_NAME}")

  ## Collect files
  aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_H  )
  aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_CPP)
  aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_HPP)

  if( (NOT GISMO_BUILD_LIB) )
    aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_INS)
    if(${subdir}_INS)
      LIST( REMOVE_ITEM ${subdir}_CPP ${${subdir}_INS})
    endif()
  endif()

  ## Add module
  add_library(${GM_NAME} OBJECT
    ${${subdir}_H}
    ${${subdir}_HPP}
    ${${subdir}_CPP}
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

  ## install module
  if (NOT ${GM_NAME} MATCHES gsMesh) # avoid installing twice gsMesh
    install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}"
            DESTINATION include/${PROJECT_NAME}
            FILES_MATCHING PATTERN "*.h" )
  endif()

ENDFOREACH()

if (GISMO_BUILD_PCH)

  include(cotire)

  #set(CMAKE_C_FLAGS    "${CMAKE_CXX_FLAGS} -H") # for GCC: test if PCH is used

  set_property(DIRECTORY PROPERTY COTIRE_ADD_UNITY_BUILD FALSE)
  LIST( REMOVE_ITEM SUBDIRS gsCore)
  set_target_properties(gsCore PROPERTIES
    COTIRE_PREFIX_HEADER_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/gsCore/gsLinearAlgebra.h"
  )
  cotire(gsCore)
  get_target_property(gismo_prefixHeader gsCore COTIRE_CXX_PREFIX_HEADER)


  FOREACH(subdir ${SUBDIRS})
    get_filename_component(GM_NAME ${subdir} NAME)
    add_dependencies(${GM_NAME} gsCore_pch)
    set_target_properties(${GM_NAME} PROPERTIES
      COTIRE_CXX_PREFIX_HEADER_INIT "${gismo_prefixHeader}")
    cotire(${GM_NAME})
  ENDFOREACH()

endif(GISMO_BUILD_PCH)

