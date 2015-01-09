### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

project(gismo-examples)

#find_package(gismo REQUIRED)
#include_directories( ${GISMO_INCLUDE_DIRS})
#link_directories   ( ${GISMO_LIBRARY_DIR} )
#include( ${GISMO_USE_FILE})

# Collect source file names
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
#message("found files: ${FILES}")

foreach(file ${FILES})
    add_gismo_executable(${file})
endforeach(file ${FILES})

list(LENGTH FILES len)
message(STATUS "Number of examples to compile: ${len}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
