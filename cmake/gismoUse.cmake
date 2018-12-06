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
    add_executable(${FNAME} ${FILE} ${gismo_SOURCES} ${gismo_EXTENSIONS} ${gismo_dev_EXTENSIONS})
    target_link_libraries(${FNAME} gismo_static)
    if(UNIX AND NOT APPLE)
    	target_link_libraries(${FNAME} dl)
    endif(UNIX AND NOT APPLE)
    # Allow CMake to follow dependencies on hpp files
    set_property(TARGET ${FNAME} PROPERTY
    IMPLICIT_DEPENDS_INCLUDE_TRANSFORM "GISMO_HPP_HEADER(%)=\"%\"")
    SET_TARGET_PROPERTIES(${FNAME} PROPERTIES COMPILE_FLAGS -UGISMO_BUILD_LIB)
endmacro(add_gismo_pure_executable)

# add an executable compiled with the shared library
macro(add_gismo_shared_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (dynamically linked): ${FNAME}")
    add_executable(${FNAME} ${FILE})
    target_link_libraries(${FNAME} gismo)
    if (GISMO_BUILD_COVERAGE)
      target_link_libraries(${FNAME} gcov)
    endif(GISMO_BUILD_COVERAGE)
    if(UNIX AND NOT APPLE)
    	target_link_libraries(${FNAME} dl)
    endif(UNIX AND NOT APPLE)
	#set_property(TARGET ${FNAME} PROPERTY FOLDER "tests-gismo")
endmacro(add_gismo_shared_executable)

# add an executable compiled statically with the library
macro(add_gismo_static_executable FILE)
    get_filename_component(FNAME ${FILE} NAME_WE) # name without extension
    add_test(${FNAME} ${CMAKE_BINARY_DIR}/bin/${FNAME} )
    #message(STATUS "exec (statically linked): ${FNAME}")
    add_executable(${FNAME} ${FILE})
    target_link_libraries(${FNAME} gismo_static)
    if(UNIX AND NOT APPLE)
    	target_link_libraries(${FNAME} dl)
    endif(UNIX AND NOT APPLE)
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
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/gs*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO(SUBDIRLIST)

# collect .cpp files
macro(aux_cpp_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/[^.]*.cpp)
endmacro(aux_cpp_directory)

# collect gs*.cpp files
macro(aux_gs_cpp_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/gs*.cpp)
endmacro(aux_gs_cpp_directory)

# collect .h files
macro(aux_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/[^.]*.h)
endmacro(aux_header_directory)

# collect .hpp files
macro(aux_tmpl_header_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/[^.]*.hpp)
endmacro(aux_tmpl_header_directory)

# collect .cxx files
macro(aux_cxx_source_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/[^.]*.cxx)
endmacro(aux_cxx_source_directory)

# collect _.cpp (instance) files
macro(aux_instance_directory DIR VAR)
	FILE(GLOB ${ARGV1} ${DIR}/[^.]*_.cpp)
endmacro(aux_instance_directory)

function(get_repo_info repository revision) #REPO_REVISION
  if (NOT DEFINED ${repository})
    if (EXISTS "${gismo_SOURCE_DIR}/.svn")
      set(${repository} "svn" PARENT_SCOPE)
      #find_program(Subversion_SVN_EXECUTABLE NAMES svn svn.bat)
      find_package(Subversion)
      set(${revision} ${Project_WC_REVISION} PARENT_SCOPE)
      #elseif (EXISTS "${gismo_SOURCE_DIR}/.git/svn/refs")
    elseif (EXISTS "${gismo_SOURCE_DIR}/.git")
      set(${repository} "git" PARENT_SCOPE)
      find_program(git_executable NAMES git git.exe git.cmd)
      execute_process(COMMAND ${git_executable} log -1 --pretty=format:%H
	WORKING_DIRECTORY ${gismo_SOURCE_DIR} TIMEOUT 5
	RESULT_VARIABLE git_res OUTPUT_VARIABLE git_rev)
      set(${revision} ${git_rev} PARENT_SCOPE)
    #else()
    #  message("GISMO_REPO undefined.")
    endif()
  endif()
    #set( ${repo_exe}
  endfunction()
