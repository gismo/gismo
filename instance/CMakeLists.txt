## 
### Commentary: 
##
## These .cpp files are used to construct a shared library eg. with
## double arithmetic
## 
######################################################################

## #################################################################
## Collect sources
## #################################################################

aux_generate_instance( ${CMAKE_CURRENT_SOURCE_DIR})

# Collect source file names (done in aux_generate_instance)
# aux_source_directory(${PROJECT_BINARY_DIR}/instance INSTANCE_SOURCES)

## #################################################################
## Add library targets
## #################################################################

#TO DO: The signature
#
#    add_library(<name> OBJECT <src>...)
#
#  creates a special "object library" target.  An object library compiles
#  source files but does not archive or link their object files into a
#  library.  Instead other targets created by add_library or
#  add_executable may reference the objects using an expression of the
#  form $<TARGET_OBJECTS:objlib> as a source, where "objlib" is the
#  object library name.  For example:
#
#    add_library(... $<TARGET_OBJECTS:objlib> ...)
#    add_executable(... $<TARGET_OBJECTS:objlib> ...)
#
#  will include objlib's object files in a library and an executable
#  along with those compiled from their own sources.  Object libraries
#  may contain only sources (and headers) that compile to object files.
#  They may contain custom commands generating such sources, but not
#  PRE_BUILD, PRE_LINK, or POST_BUILD commands.  Object libraries cannot
#  be imported, exported, installed, or linked.


if(GISMO_BUILD_SHARED_LIB)

 add_library(gismo SHARED
  ${gismo_SOURCES}
  ${INSTANCE_SOURCES}
  ${gismo_HEADERS}
  ${gismo_TMPL_HEADERS}
 )

  set_target_properties(gismo PROPERTIES 
  PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/gismo.h" 
  COMPILE_DEFINITIONS GISMO_BUILD_SHARED_LIB # Used for DLL exporting on windows
)

  target_link_libraries(gismo ${gismo_EXTENSIONS})

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
    TARGET_LINK_LIBRARIES(gismo ${DBGHELP_LIBRARY}) 	
  ENDIF() 

 set_target_properties(gismo 
 PROPERTIES LINKER_LANGUAGE CXX)

  if( WIN32 ) # Copy the dll to the bin folder to allow executables to find it

    #Find the dll
    FIND_FILE(DLLFILE gismo.dll ${CMAKE_BINARY_DIR} 
    	      PATH_SUFFIXES lib lib/${CMAKE_BUILD_TYPE}
              NO_DEFAULT_PATH )

    if(NOT ${DLLFILE} MATCHES "DLLFILE-NOTFOUND") 
    add_custom_command(
    TARGET gismo
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${DLLFILE} ${CMAKE_BINARY_DIR}/bin
    )
    endif()

  endif( WIN32 )

endif(GISMO_BUILD_SHARED_LIB)

if(GISMO_BUILD_STATIC_LIB)

 add_library(gismo_static STATIC
  ${gismo_SOURCES}
  ${INSTANCE_SOURCES}
  ${gismo_HEADERS}
  ${gismo_TMPL_HEADERS}
  )

  set_target_properties(gismo_static PROPERTIES 
  COMPILE_DEFINITIONS GISMO_BUILD_STATIC_LIB)
 
  target_link_libraries(gismo_static ${gismo_EXTENSIONS})

  IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
    SET_TARGET_PROPERTIES(gismo_static PROPERTIES COMPILE_FLAGS -fPIC)
  ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
     TARGET_LINK_LIBRARIES(gismo_static ${DBGHELP_LIBRARY}) 	
  ENDIF() 

  set_target_properties(gismo_static 
  PROPERTIES LINKER_LANGUAGE CXX
  OUTPUT_NAME gismo_static )

endif(GISMO_BUILD_STATIC_LIB)


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
