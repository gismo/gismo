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

cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()
if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 NEW)
endif()
if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

project(main)

include(CTest)
#set(BUILD_TESTING "Enable testing" ${GISMO_BUILD_UNITTESTS})

set(CMAKE_CXX_STANDARD 11)
message(STATUS "CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}")

# G+Smo
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )

include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})
include(gsConfig) #compiler configuration

# Application files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
#message("found files: ${FILES}")

foreach(file ${FILES})

  add_gismo_shared_executable(${file})
  get_filename_component(tarname ${file} NAME_WE) # name without extension
  set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
  set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")

  # install the example executables (optionally)
  install(TARGETS ${tarname} DESTINATION bin COMPONENT exe OPTIONAL)

endforeach(file ${FILES})

message(STATUS "G+Smo shared library found in ${GISMO_LIBRARY_DIR}")
message(STATUS "Type \"make\" to build")
