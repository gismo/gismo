### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Author: Harald Weiner
## Copyright (C) 2012-2015 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

project(examples)

# Add a grouping target that builds all tests this is equivalent to
# running "make" inside the tests build folder
add_custom_target(${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY LABELS "${PROJECT_NAME};gismo")

#find_package(gismo REQUIRED)
#include_directories( ${GISMO_INCLUDE_DIRS})
#link_directories   ( ${GISMO_LIBRARY_DIR} )
#include( ${GISMO_USE_FILE})

# Collect source file names
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
#message("found files: ${FILES}")

foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME};gismo")
    add_dependencies(${PROJECT_NAME} ${tarname})
endforeach(file ${FILES})

list(LENGTH FILES len)
message(STATUS "Number of examples to compile: ${len}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

include( CMakeEmscripten )

# install the example executables to <CMAKE_INSTALL_PREFIX>/bin/
message("install_bin_dir='${CMAKE_INSTALL_PREFIX}/${BIN_INSTALL_DIR}/'")
install(DIRECTORY "${CMAKE_BINARY_DIR}/bin/"
		DESTINATION "${CMAKE_INSTALL_PREFIX}/${BIN_INSTALL_DIR}/"
		FILE_PERMISSIONS
					OWNER_READ OWNER_WRITE OWNER_EXECUTE
					GROUP_READ GROUP_WRITE GROUP_EXECUTE
					WORLD_READ WORLD_EXECUTE
		PATTERN "*.*"
		PERMISSIONS
					OWNER_READ OWNER_WRITE OWNER_EXECUTE
					GROUP_READ GROUP_WRITE GROUP_EXECUTE
					WORLD_READ WORLD_EXECUTE
)
