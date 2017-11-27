######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 OLD)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

project(doc-snippets)

include(gsConfig)

include(CTest)
enable_testing()


# Collect source file names
file(GLOB snippets_SRCS "*.cpp")

add_custom_target(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "unittests")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

foreach(snippet_src ${snippets_SRCS})

  get_filename_component(snippet_target ${snippet_src} NAME_WE)
  #set(snippet_target snippet_${snippet_target})

  file(READ ${snippet_src} snippet_source_code)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compile_snippet.cpp.in
                 ${CMAKE_CURRENT_BINARY_DIR}/${snippet_target}.cpp)

  add_executable(${snippet_target}  ${CMAKE_CURRENT_BINARY_DIR}/${snippet_target}.cpp)
  add_test(${snippet_target} ${EXECUTABLE_OUTPUT_PATH}/${snippet_target})
  set_target_properties(${snippet_target} PROPERTIES LABELS "${PROJECT_NAME};unittests" FOLDER "unittests")

  if(GISMO_BUILD_LIB)
    target_link_libraries(${snippet_target} gismo)
  endif()

if(DOXYGEN_FOUND)
  add_custom_command(
    TARGET ${snippet_target}
    POST_BUILD
    COMMAND ${snippet_target}
    ARGS >${CMAKE_CURRENT_BINARY_DIR}/${snippet_target}.out
  )
endif(DOXYGEN_FOUND)

  add_dependencies(${PROJECT_NAME} ${snippet_target})

  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${snippet_src}
                              PROPERTIES OBJECT_DEPENDS ${snippet_src})
endforeach(snippet_src)

#list(LENGTH snippets_SRCS len)
#message(STATUS "Number of snippets to compile: ${len}")
