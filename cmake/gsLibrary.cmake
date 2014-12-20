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

#message("Using ${gismo_EXTENSIONS}")
#message("Using ${gismo_MODULES}")
#message("Using ${gismo_SOURCES}")

if(GISMO_BUILD_LIB)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
 set(gismo_SOURCES ${gismo_SOURCES} "${PROJECT_SOURCE_DIR}/src/misc/gsLibInit.cpp")
endif()

#add_library(gismo ${LIB_TYPE}
#           ${gismo_MODULES}
#           ${gismo_EXTENSIONS}
#           GENERATE_EXPORT_HEADER( gismo
#           BASE_NAME gismo
#           EXPORT_MACRO_NAME gismo_EXPORT
#           EXPORT_FILE_NAME gsExport.h
#           STATIC_DEFINE GISMO_BUILT_AS_STATIC
#	   )

  add_library(gismo SHARED
    #${gismo_MODULES}
    ${gismo_SOURCES}
    ${gismo_EXTENSIONS}
    )

  set_target_properties(gismo PROPERTIES 
  PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/gismo.h" 
  POSITION_INDEPENDENT_CODE ON
  COMPILE_DEFINITIONS GISMO_BUILD_SHARED_LIB # Used for DLL exporting on windows
)

  if (GISMO_WITH_PSOLID)
    target_link_libraries(gismo ${PARASOLID_LIBRARY})
  endif()

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
    TARGET_LINK_LIBRARIES(gismo ${DBGHELP_LIBRARY}) 	
  ENDIF() 

 set_target_properties(gismo 
 PROPERTIES LINKER_LANGUAGE CXX)

if( WIN32 ) # Copy the dll to the bin folder to allow executables to find it
    add_custom_command(
    TARGET gismo
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:gismo> ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:gismo> is copied in the bin folder for convenience.' )
endif( WIN32 )

  add_library(gismo_static STATIC
  #${gismo_MODULES}
  ${gismo_SOURCES}
  ${gismo_EXTENSIONS}
  )

  set_target_properties(gismo_static PROPERTIES 
  POSITION_INDEPENDENT_CODE ON
  EXCLUDE_FROM_ALL 1
  EXCLUDE_FROM_DEFAULT_BUILD 1
  )
 
  IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    SET_TARGET_PROPERTIES(gismo_static PROPERTIES 
                         POSITION_INDEPENDENT_CODE ON)
  ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
     TARGET_LINK_LIBRARIES(gismo_static ${DBGHELP_LIBRARY}) 	
  ENDIF() 

  set_target_properties(gismo_static 
  PROPERTIES LINKER_LANGUAGE CXX
  OUTPUT_NAME gismo_static )

endif(GISMO_BUILD_LIB)


set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/)


## #################################################################
## Installation
## #################################################################

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")


install(TARGETS ${PROJECT_NAME}
  # IMPORTANT: Add the gismo library to the "export-set"
  EXPORT gismoTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/gismo"
  #ARCHIVE DESTINATION lib
  )
