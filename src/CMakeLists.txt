######################################################################
## CMakeLists.txt --- src
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

include_directories(${GISMO_INCLUDE_DIRS})

if(GISMO_WITH_CODIPACK)
  include_directories(${CODIPACK_INCLUDE_DIR})
endif()

if(GISMO_WITH_UNUM)
  include_directories(${UNUM_INCLUDE_DIR})
endif()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
#LIST( REMOVE_ITEM SUBDIRS misc)
LIST( APPEND SUBDIRS gsUtils/gsMesh)

if (GISMO_BUILD_PCH)
  add_precompiled_header(gsPrecompiledHeader misc/gsPrecompiledHeader.h)
endif()

FOREACH(subdir ${SUBDIRS})

  get_filename_component(GM_NAME ${subdir} NAME)
  #message("module: ${GM_NAME}")

  ## Collect files
  aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_H  )
  aux_cpp_noins_directory  (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_CPP)
  aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_HPP)

  if(GISMO_BUILD_LIB)
    aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_INS)
  endif()

  ## Add module
  add_library(${GM_NAME} OBJECT
    ${${subdir}_H}
    ${${subdir}_HPP}
    ${${subdir}_CPP} # static/non templated part
    ${${subdir}_INS}
    )

  set_target_properties(${GM_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    FOLDER "G+Smo modules"
    )

  ## Add extra instances (template instantiations)
  # NOTE: The static/non templated part of the library is
  # compiled above only for GISMO_COEFF_TYPE
  math(EXPR ii "1")
  foreach(ins ${GISMO_EXTRA_INSTANCE})
    add_library(${GM_NAME}${ii} OBJECT
      ${${subdir}_H}
      ${${subdir}_HPP}
      ${${subdir}_INS}
      )
    set_target_properties(${GM_NAME}${ii} PROPERTIES
      COMPILE_DEFINITIONS "gismo_EXPORTS;real_t=${ins}"
      POSITION_INDEPENDENT_CODE ON
      LINKER_LANGUAGE CXX
      FOLDER "G+Smo modules"
      )
    set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}${ii}>
      CACHE INTERNAL "G+Smo modules" )
    math(EXPR ii "${ii}+1")
  endforeach()

  if(GISMO_BUILD_MULTIPRECISION)
    add_dependencies(${GM_NAME} gsMultiPrecisionExtension)
  endif()

  if (GISMO_BUILD_PCH)
    target_precompiled_header(${GM_NAME} gsPrecompiledHeader)
  endif()

  set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}>
    CACHE INTERNAL "G+Smo modules" )

  ## install module
  if (NOT ${GM_NAME} MATCHES gsMesh) # avoid installing twice gsMesh
    install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}"
      DESTINATION include/${PROJECT_NAME}
      FILES_MATCHING PATTERN "*.h" )
  endif()

ENDFOREACH()
