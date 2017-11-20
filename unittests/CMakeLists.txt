######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(unittests)

include(gsConfig)

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

add_executable(unittests ${unittests_SRCS} common/main.cpp common/test.h
  common/CheckFiles.h common/CheckFiles.cpp 
  common/CheckMacros.h common/CheckMatrix.h 
  common/SuppressOutput.h common/SuppressOutput.cpp 
  common/random_refinement.h)
target_link_libraries(unittests ${UnitTestPP_LIBRARIES} gismo)
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
