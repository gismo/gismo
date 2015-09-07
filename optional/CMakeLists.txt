### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

project(optional)
#add_custom_target(${PROJECT_NAME})
#set_property(TARGET ${PROJECT_NAME} PROPERTY LABELS ${PROJECT_NAME})

## Parasolid option
if(GISMO_WITH_PSOLID)

  FILE(GLOB parasolid_exe ${CMAKE_CURRENT_SOURCE_DIR}/gsParasolid*.cpp)
  foreach(file ${parasolid_exe})
    add_gismo_executable(${file})
    mark_gismo_optional (${file})
    get_filename_component(namewe ${file} NAME_WE)
    add_dependencies(examples ${namewe})
    set_property(TEST ${namewe} PROPERTY LABELS examples)
  endforeach()
  
  #set_property(TARGET gsParasolid_test APPEND PROPERTY 
  #INCLUDE_DIRECTORIES ${PARASOLID_INCLUDE_DIR} )

endif(GISMO_WITH_PSOLID)


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
