######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

include_directories(${GISMO_INCLUDE_DIRS})

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
LIST( REMOVE_ITEM SUBDIRS misc)

FOREACH(subdir ${SUBDIRS})

  #message("module: ${subdir}")

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
  add_library(${subdir} OBJECT
    ${${subdir}_H}
    ${${subdir}_HPP}
    ${${subdir}_CPP}
    )

  set_target_properties(${subdir} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    FOLDER "G+Smo modules"
    #EXCLUDE_FROM_ALL 1
    #EXCLUDE_FROM_DEFAULT_BUILD 1
    )

  set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${subdir}>
      CACHE INTERNAL "G+Smo modules" )

  ## install module
  install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}"
          DESTINATION include/${PROJECT_NAME}/
          FILES_MATCHING PATTERN "*.h"
          )

ENDFOREACH()

if (GISMO_BUILD_PCH)

  #set(CMAKE_C_FLAGS    "${CMAKE_CXX_FLAGS} -H") # for GCC: test if PCH is used

  include(cotire)
  set_property(DIRECTORY PROPERTY COTIRE_ADD_UNITY_BUILD FALSE)
  LIST( REMOVE_ITEM SUBDIRS gsCore)
  set_target_properties(gsCore PROPERTIES
    COTIRE_PREFIX_HEADER_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/gsCore/gsLinearAlgebra.h"
  )
  cotire(gsCore)
  get_target_property(gismo_prefixHeader gsCore COTIRE_CXX_PREFIX_HEADER)

  FOREACH(subdir ${SUBDIRS})

    add_dependencies(${subdir} gsCore_pch)
    set_target_properties(${subdir} PROPERTIES
      COTIRE_CXX_PREFIX_HEADER_INIT "${gismo_prefixHeader}")
    #cotire(${subdir})

  ENDFOREACH()

cotire(${SUBDIRS})

endif(GISMO_BUILD_PCH)
