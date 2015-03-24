######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

#include (GenerateExportHeader)

## #################################################################
## Add library targets
## #################################################################

#message("Using ${${PROJECT_NAME}_EXTENSIONS}")
#message("Using ${${PROJECT_NAME}_MODULES}")
#message("Using ${${PROJECT_NAME}_SOURCES}")

if(GISMO_BUILD_LIB)

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-implicit-templates")
#endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
 set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES} 
     "${gismo_SOURCE_DIR}/src/misc/gsDllMain.cpp")
endif()

#add_library(${PROJECT_NAME} ${LIB_TYPE}
#           ${${PROJECT_NAME}_MODULES}
#           ${${PROJECT_NAME}_EXTENSIONS}
#           GENERATE_EXPORT_HEADER( ${PROJECT_NAME}
#           BASE_NAME ${PROJECT_NAME}
#           EXPORT_MACRO_NAME ${PROJECT_NAME}_EXPORTS
#           EXPORT_FILE_NAME gsExport.h
#           STATIC_DEFINE GISMO_BUILT_STATIC_LIB
#	   )
#
#  add_library(${PROJECT_NAME}_obj OBJECT
#    #${${PROJECT_NAME}_MODULES}
#    ${${PROJECT_NAME}_SOURCES}
#    )
#
#  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    set_target_properties(${PROJECT_NAME}_obj PROPERTIES 
#    DEFINE_SYMBOL  "${PROJECT_NAME}_EXPORTS"
#    DEFINE_SYMBOL  "GISMO_BUILD_SHARED_LIB"
#    )
#  endif()
#
#    GENERATE_EXPORT_HEADER(${PROJECT_NAME}_obj
#    #BASE_NAME ${PROJECT_NAME}_obj
#    #EXPORT_MACRO_NAME ${PROJECT_NAME}_obj_EXPORTS
#    EXPORT_FILE_NAME gsExport.h
#    #STATIC_DEFINE GISMO_BUILT_STATIC_LIB
#    )
#
#  set_target_properties(${PROJECT_NAME}_obj PROPERTIES 
#  POSITION_INDEPENDENT_CODE ON
#)
#
#  add_library(${PROJECT_NAME} SHARED
#    #${${PROJECT_NAME}_MODULES}
#    $<TARGET_OBJECTS:${PROJECT_NAME}_obj>
#    ${${PROJECT_NAME}_EXTENSIONS}
#    )
#
#  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    set_target_properties(${PROJECT_NAME}_obj PROPERTIES 
#    UNDEFINE_SYMBOL  "${PROJECT_NAME}_EXPORTS"
#    DEFINE_SYMBOL  "GISMO_BUILD_SHARED_LIB"
#    )
#  endif()

  add_library(${PROJECT_NAME} SHARED
    ${${PROJECT_NAME}_MODULES}
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_EXTENSIONS}
    )

  set_target_properties(${PROJECT_NAME} PROPERTIES 
  PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}.h" 
  POSITION_INDEPENDENT_CODE ON
  #COMPILE_DEFINITIONS ${PROJECT_NAME}_EXPORTS # Used for DLL exporting (defined by default by CMake)
  )

  if (GISMO_WITH_ONURBS AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
    target_link_libraries(${PROJECT_NAME} Rpcrt4)
  endif()

  if (GISMO_WITH_PSOLID)
    target_link_libraries(${PROJECT_NAME} ${PARASOLID_LIBRARY})
  endif()

  if (GISMO_WITH_SUPERLU)
    target_link_libraries(${PROJECT_NAME} ${SUPERLU_LIBRARIES})
  endif()

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${DBGHELP_LIBRARY}) 	
  ENDIF() 

 set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

if( WIN32 ) # Copy the dll to the bin folder to allow executables to find it
    add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:${PROJECT_NAME}> is copied in the bin folder for convenience.' )
endif( WIN32 )

endif(GISMO_BUILD_LIB)

  add_library(${PROJECT_NAME}_static STATIC
  ${${PROJECT_NAME}_MODULES}
  #$<TARGET_OBJECTS:${PROJECT_NAME}_obj>
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_EXTENSIONS}
  )

  set_target_properties(${PROJECT_NAME}_static PROPERTIES 
  POSITION_INDEPENDENT_CODE ON
  EXCLUDE_FROM_ALL 1
  EXCLUDE_FROM_DEFAULT_BUILD 1
  COMPILE_DEFINITIONS ${PROJECT_NAME}_STATIC
  )
 
  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
     TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static ${DBGHELP_LIBRARY}) 	
  ENDIF() 

  set_target_properties(${PROJECT_NAME}_static 
  PROPERTIES LINKER_LANGUAGE CXX
  OUTPUT_NAME ${PROJECT_NAME}_static )


set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/)


## #################################################################
## Installation
## #################################################################

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")


if(GISMO_BUILD_LIB)

  install(TARGETS ${PROJECT_NAME}
  # IMPORTANT: Add the ${PROJECT_NAME} library to the "export-set"
  EXPORT gismoTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}"
  #ARCHIVE DESTINATION lib
  )

endif(GISMO_BUILD_LIB)