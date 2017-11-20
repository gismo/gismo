######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(unittests)

include(gsConfig)

include(ExternalProject)
ExternalProject_Add(gsUnitTest
          SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest
          URL https://github.com/gismo/gsUnitTest/archive/master.zip
          URL_MD5 "913dc7cc931c89ff3e561f7f43644eee"
          #DOWNLOAD_NO_PROGRESS 0
          CMAKE_ARGS   -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DCMAKE_BUILD_TYPE=Release
)
set(GSUNITTEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest CACHE INTERNAL "")
set(GSUNITTEST_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/gsUnitTest-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}UnitTestPP${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")

include_directories(${GSUNITTEST_INCLUDE_DIR})

include(CTest)
enable_testing()

# Collect source file names
FILE (GLOB unittests_SRCS *.cpp)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)

add_executable(unittests ${unittests_SRCS} main.cpp)
target_link_libraries(unittests ${GSUNITTEST_LIBRARY} gismo)
set_property(TARGET unittests PROPERTY LABELS unittests)

add_dependencies(${PROJECT_NAME} gsUnitTest)


foreach(file ${unittests_SRCS})
    get_filename_component(testname ${file} NAME_WE) # name without extension
    add_test(NAME unit_${testname} COMMAND $<TARGET_FILE:unittests> ${file})
    set_property(TEST unit_${testname} PROPERTY LABELS unittests)
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

if (NOT NOSNIPPETS)
  add_dependencies(unittests doc-snippets)
endif()
