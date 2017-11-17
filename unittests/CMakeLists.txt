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

project(unittests)

# TO DO: Use file
include( ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/gsConfig.cmake)

find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )
include( ${GISMO_USE_FILE})

find_package (gismo_dev REQUIRED)
include_directories( ${GISMO_DEV_INCLUDE_DIRS} )
link_directories   ( ${GISMO_DEV_LIBRARY_DIR} )

add_subdirectory(UnitTest++)

# to do UnitTestPPConfig.cmake
#FIND_PACKAGE (UnitTestPP REQUIRED)
#INCLUDE_DIRECTORIES(${UnitTestPP_INCLUDE_DIRS})
#LINK_DIRECTORIES(${UnitTestPP_LIBRARY_DIRS})
set(UnitTestPP_LIBRARIES UnitTestPP)

include(CTest)
enable_testing()

# Collect source file names
FILE (GLOB unittests_SRCS *.cpp)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)
add_definitions(-DGISMO_DDATA_DIR="${PARENT_DIR}/devel/filedata/")

add_executable(unittests ${unittests_SRCS} common/main.cpp common/test.h
  common/CheckFiles.h common/CheckFiles.cpp 
  common/CheckMacros.h common/CheckMatrix.h 
  common/SuppressOutput.h common/SuppressOutput.cpp 
  common/random_refinement.h)
target_link_libraries(unittests ${UnitTestPP_LIBRARIES} ${GISMO_LIBRARIES} ${GISMO_DEV_LIBRARIES})
set_property(TARGET unittests PROPERTY LABELS unittests)

foreach(file ${unittests_SRCS})
    get_filename_component(testname ${file} NAME_WE) # name without extension
    add_test(NAME unit_${testname} COMMAND $<TARGET_FILE:unittests> ${file})
    set_property(TEST unit_${testname} PROPERTY LABELS unittests)
endforeach(file ${FILES})

#add_test(unittests ${CMAKE_BINARY_DIR}/bin/unittests)
#add_custom_target(check COMMAND $<TARGET_FILE:unittests> {CMAKE_CTEST_COMMAND} -V)
#add_dependencies(check unittests)

# run unit tests as post build step
#add_custom_command(TARGET unittests
# POST_BUILD COMMAND $<TARGET_FILE:unittests>
# COMMENT "Running unit tests")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
