######################################################################
##
## CMakeLists.txt ---
## Sample CMake project using G+Smo
##
## To ensure that G+Smo is located correctly, set the variable gismo_DIR
## to a path containing the binary and include files, eg:
##
## cmake . -Dgismo_DIR=/path/to/gismo/build
##
######################################################################

cmake_minimum_required(VERSION 2.8.8)

project(main)

# G+Smo
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )
include( ${GISMO_USE_FILE})
message(STATUS "G+Smo shared library in ${GISMO_LIBRARY_DIR}")

# Source file
add_executable(main main.cpp)
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(main gismo)
message(STATUS "Type \"make\" to build main")
