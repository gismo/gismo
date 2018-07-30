######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(doc-snippets)

include(gsConfig)

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

find_package(Doxygen)
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
