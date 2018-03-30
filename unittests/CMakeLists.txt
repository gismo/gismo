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
          URL_MD5 "73ca72375178fca77fa22dbbf3624526"
          #DOWNLOAD_NO_PROGRESS 0
          CMAKE_ARGS   -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
          BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}UnitTestPP${CMAKE_STATIC_LIBRARY_SUFFIX}
)
ExternalProject_Get_Property(gsUnitTest INSTALL_DIR)
add_library(gsUnitTest_lib STATIC IMPORTED GLOBAL)
set_property(TARGET gsUnitTest_lib PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}UnitTestPP${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GSUNITTEST_INCLUDE_DIR ${INSTALL_DIR}/include/UnitTest++ CACHE INTERNAL "")
include_directories(${GSUNITTEST_INCLUDE_DIR})

include(CTest)
enable_testing()

# Collect source file names
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR} unittests_SRCS)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)


add_executable(${PROJECT_NAME} gismo_unittest.h ${unittests_SRCS} main.cpp)
set_property(TARGET ${PROJECT_NAME} PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES gsUnitTest_lib)
target_link_libraries(${PROJECT_NAME} gismo gsUnitTest_lib)
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}")

#CTEST_BUILD_TARGET=ALL
add_dependencies(${PROJECT_NAME}  gsUnitTest) # gsUnitTest_lib
if (NOT NOSNIPPETS)
  add_dependencies(${PROJECT_NAME} doc-snippets) # doc-snippets are tests as well
endif()

foreach(file ${unittests_SRCS})
    get_filename_component(testname ${file} NAME_WE) # name without extension
    add_test(NAME unit_${testname} COMMAND $<TARGET_FILE:unittests> ${file})
    set_property(TEST unit_${testname} PROPERTY LABELS unittests)
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
