######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(unittests)

include(gsConfig)

#add_subdirectory(UnitTest++)
#if (NOT IPOPT_FOUND)
#  set(IPOPT_VER "3.12.6")
#  message("IpOpt not found (define IpOpt_DIR if it exists on the system).")
#  message("Will attempt to download and compile IpOpt ${IPOPT_VER} (do: make IpOpt; cmake .")
include(ExternalProject)
ExternalProject_Add(gsUnitTest
          SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest
          URL https://github.com/gismo/gsUnitTest/archive/master.zip
          URL_MD5 "a20cd363fa98760a097ab5f56d54774c"
          #DOWNLOAD_NO_PROGRESS 0
          CMAKE_ARGS   -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
#-DCMAKE_BUILD_TYPE=Release
)
      set(GSUNITTEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest CACHE INTERNAL "")
      set(GSUNITTEST_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest-prefix/lib/libUnitTestPP.a CACHE INTERNAL "")

#else()
#    #add_custom_target (IpOpt)
#endif(NOT IPOPT_FOUND)

include_directories(${GSUNITTEST_INCLUDE_DIR})

# to do UnitTestPPConfig.cmake
#FIND_PACKAGE (UnitTestPP REQUIRED)
#INCLUDE_DIRECTORIES(${UnitTestPP_INCLUDE_DIRS})
#LINK_DIRECTORIES(${UnitTestPP_LIBRARY_DIRS})
#set(UnitTestPP_LIBRARIES UnitTestPP)

include(CTest)
enable_testing()

# Collect source file names
FILE (GLOB unittests_SRCS *.cpp)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)

add_executable(unittests ${unittests_SRCS} main.cpp)
#  gs/main.cpp gs/test.h
#  gs/CheckFiles.h gs/CheckFiles.cpp 
#  gs/CheckMacros.h gs/CheckMatrix.h 
#  gs/SuppressOutput.h gs/SuppressOutput.cpp 
#  gs/random_refinement.h)
target_link_libraries(unittests ${GSUNITTEST_LIBRARY} gismo)
set_property(TARGET unittests PROPERTY LABELS unittests)

add_dependencies(${PROJECT_NAME} gsUnitTest)


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
