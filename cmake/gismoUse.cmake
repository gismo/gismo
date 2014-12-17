######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

#macro to add an executable
macro(add_gismo_executable FILE)

    if( GISMO_BUILD_SHARED_LIB )
        add_gismo_shared_executable(${FILE})
    elseif ( GISMO_BUILD_STATIC_LIB )
        add_gismo_static_executable(${FILE})
    else ( GISMO_BUILD_STATIC_LIB )
    	add_gismo_pure_executable(${FILE})
    endif( GISMO_BUILD_SHARED_LIB )

endmacro(add_gismo_executable)

#macro to add an executable compiled with pure template headers
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

#macro to add an executable compiled with the shared library
macro(add_gismo_shared_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (dynamically linked): ${FNAME}")
    add_executable(${FNAME} ${FNAME})
    target_link_libraries(${FNAME} gismo)
    
    # For devel
    #target_link_libraries(${FNAME} gismo-dev)

    if (GISMO_BUILD_COVERAGE)
      target_link_libraries(${FNAME} gcov)
    endif(GISMO_BUILD_COVERAGE)

endmacro(add_gismo_shared_executable)

#macro to add an executable compiled statically with the library
macro(add_gismo_static_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (statically linked): ${FNAME}")
    add_executable(${FNAME} ${FNAME})
    target_link_libraries(${FNAME} gismo_static)
endmacro(add_gismo_static_executable)

#macro to add source files when needed
macro(add_gismo_sources VAR)
  if( (NOT GISMO_BUILD_STATIC_LIB) AND  (NOT GISMO_BUILD_SHARED_LIB) )
    set(${VAR} ${VAR} ${gismo_SOURCES} )
  endif( (NOT GISMO_BUILD_STATIC_LIB) AND  (NOT GISMO_BUILD_SHARED_LIB) )
endmacro(add_gismo_sources)

#macro for h
macro(aux_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.h)
endmacro(aux_header_directory)

#macro for hpp
macro(aux_tmpl_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.hpp)
endmacro(aux_tmpl_header_directory)

#macro for cxx
macro(aux_cxx_source_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*.cxx)
endmacro(aux_cxx_source_directory)

#add to the list
macro(add_gismo_test FILE)
    add_gismo_executable(${file})
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    set(gismo_tests ${gismo_tests} ${FNAME}
     CACHE INTERNAL "${PROJECT_NAME} list of tests")
endmacro(add_gismo_test)

macro(add_gismo_example FILE)
    add_gismo_executable(${file})
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    set(gismo_examples ${gismo_examples} ${FNAME}
     CACHE INTERNAL "${PROJECT_NAME} list of examples")
endmacro(add_gismo_example)

macro(mark_gismo_optional FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    set(gismo_optionals ${gismo_optionals} ${FNAME}
     CACHE INTERNAL "${PROJECT_NAME} list of optionals")
endmacro(mark_gismo_optional)

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


#macro for generating instance files
macro(aux_generate_instance DIR)
	FILE(GLOB IFILES ${DIR}/*_.cpp.in)

  foreach(file ${IFILES})

    get_filename_component(FNAME ${file} NAME_WE) # name without extension

    configure_file ("${file}"
                    "${PROJECT_BINARY_DIR}/instance/${FNAME}.cpp" )

    set(INSTANCE_SOURCES ${INSTANCE_SOURCES} "${PROJECT_BINARY_DIR}/instance/${FNAME}.cpp"
     CACHE INTERNAL "${PROJECT_NAME} instance files to be compiled")

  endforeach(file ${IFILES})

endmacro(aux_generate_instance)

#macro for collecting instance files
macro(aux_instance_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/*_.cpp)
endmacro(aux_instance_directory)
