######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(unittests)

include(gsConfig)

#update all:
#git submodule update --init --recursive
#git pull --recurse-submodules

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)

#if((NOT GSUNITTEST_INCLUDE_DIR) OR (NOT EXISTS ${GSUNITTEST_INCLUDE_DIR}))
    # clone submodule if its not there
execute_process(COMMAND git submodule update --init -- extensions/gsUnitTest
  WORKING_DIRECTORY ${PARENT_DIR} )

if( NOT EXISTS "${PARENT_DIR}/extensions/gsUnitTest/CMakeLists.txt")
message(SEND_ERROR "The git submodule is not available.")
endif()

# set include directory
set(GSUNITTEST_INCLUDE_DIR ${PARENT_DIR}/extensions/gsUnitTest
  CACHE INTERNAL "gsUnitTest include directory")

install(DIRECTORY ${GSUNITTEST_INCLUDE_DIR} DESTINATION
  DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.h")

# header
#  add_library(ho_lib INTERFACE)
#  target_include_directories(foo INTERFACE
#    $<BUILD_INTERFACE:${GSUNITTEST_INCLUDE_DIR}>
#    $<INSTALL_INTERFACE:include/${PROJECT_NAME}/ho_lib>)
#  # need to export target as well
#  #install(TARGETS ho_lib EXPORT the_export_set DESTINATION ${some_dest})

    # build the project (if not header-only)
    add_subdirectory(${PARENT_DIR}/extensions/gsUnitTest ${CMAKE_BINARY_DIR}/extensions/gsUnitTest)
#endif()

include_directories(${GSUNITTEST_INCLUDE_DIR})

# Collect source file names
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR} unittests_SRCS)

add_executable(${PROJECT_NAME} gismo_unittest.h ${unittests_SRCS} main.cpp)
target_link_libraries(${PROJECT_NAME} gismo UnitTestPP)
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}")

if (NOT NOSNIPPETS)
  add_dependencies(${PROJECT_NAME} doc-snippets) # doc-snippets are tests as well
endif()

foreach(file ${unittests_SRCS})
    get_filename_component(testname ${file} NAME_WE) # name without extension
    add_test(NAME unit_${testname} COMMAND $<TARGET_FILE:unittests> ${file})
    set_property(TEST unit_${testname} PROPERTY LABELS unittests)
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
