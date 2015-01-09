######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

######################################################################
# Macros for adding executables linked to G+Smo
######################################################################

# add an executable
macro(add_gismo_executable FILE)
    if( GISMO_BUILD_LIB )
        add_gismo_shared_executable(${FILE})
    #elseif ( GISMO_BUILD_STATIC_LIB )
    #    add_gismo_static_executable(${FILE})
    else ( GISMO_BUILD_LIB )
    	add_gismo_pure_executable(${FILE})
    endif( GISMO_BUILD_LIB )
endmacro(add_gismo_executable)

# add an executable compiled with pure template headers
macro(add_gismo_pure_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (pure template): ${FNAME}")
    add_executable(${FNAME} ${FNAME} ${gismo_SOURCES} ${gismo_EXTENSIONS})
    # Allow CMake to follow dependencies on hpp files
    set_property( TARGET ${FNAME} PROPERTY 
    IMPLICIT_DEPENDS_INCLUDE_TRANSFORM "GISMO_HPP_HEADER(%)=\"%\"")
    SET_TARGET_PROPERTIES(${FNAME} PROPERTIES COMPILE_FLAGS -UGISMO_BUILD_LIB)
endmacro(add_gismo_pure_executable)

# add an executable compiled with the shared library
macro(add_gismo_shared_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (dynamically linked): ${FNAME}")
    add_executable(${FNAME} ${FNAME})
    target_link_libraries(${FNAME} gismo)
    if (GISMO_BUILD_COVERAGE)
      target_link_libraries(${FNAME} gcov)
    endif(GISMO_BUILD_COVERAGE)
endmacro(add_gismo_shared_executable)

# add an executable compiled statically with the library
macro(add_gismo_static_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (statically linked): ${FNAME}")
    add_executable(${FNAME} ${FNAME})
    target_link_libraries(${FNAME} gismo_static)
endmacro(add_gismo_static_executable)

macro(mark_gismo_optional FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    set(gismo_optionals ${gismo_optionals} ${FNAME}
     CACHE INTERNAL "${PROJECT_NAME} list of optionals")
endmacro(mark_gismo_optional)

######################################################################
# File collection utilities
######################################################################

# list all subdirectories of the current directory
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        SET(dirlist ${dirlist} ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO(SUBDIRLIST)

# collect .cpp files
#aux_source_directory

# collect .h files
macro(aux_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.h)
endmacro(aux_header_directory)

# collect .hpp files
macro(aux_tmpl_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.hpp)
endmacro(aux_tmpl_header_directory)

# collect .cxx files
macro(aux_cxx_source_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.cxx)
endmacro(aux_cxx_source_directory)

# collect _.cpp (instance) files
macro(aux_instance_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*_.cpp)
endmacro(aux_instance_directory)
