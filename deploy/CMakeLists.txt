######################################################################
## CMakeLists.txt ---
## Sample project using G+Smo
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
