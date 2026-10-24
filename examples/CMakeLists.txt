######################################################################
## CMakeLists.txt --- examples
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Harald Weiner
######################################################################

project(examples)

set(CMAKE_DIRECTORY_LABELS "${PROJECT_NAME}") #CMake 3.10
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

# Add a grouping target that builds all tests this is equivalent to
# running "make" inside the tests build folder
add_custom_target(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}")

#find_package(gismo REQUIRED)
#include_directories( ${GISMO_INCLUDE_DIRS})
#link_directories   ( ${GISMO_LIBRARY_DIR} )
#include( ${GISMO_USE_FILE})

# Collect source file names
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
#message("found files: ${FILES}")

if ("x${CMAKE_CXX_STANDARD}" STREQUAL "x98")
  set(IGNORE_FILES "kirchhoff-Love_example;")
  foreach(ign ${IGNORE_FILES})
    foreach(file ${FILES})
      string(REGEX MATCH "${ign}.cpp$" fmatch ${file})
      if(fmatch)
	list(REMOVE_ITEM FILES ${file})
      endif()
    endforeach()
  endforeach()
endif()

foreach(file ${FILES})
  get_filename_component(tarname ${file} NAME_WE) # name without extension
  add_gismo_executable(${file})
  set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
  set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
  add_dependencies(${PROJECT_NAME} ${tarname})

  # install the example executables (optionally)
  install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

if(GISMO_BUILD_UNITTESTS AND TARGET rh_refinement_example)
  # Keep generated files separate from other examples and avoid launching a viewer.
  set(rh_surface_test_dir "${CMAKE_CURRENT_BINARY_DIR}/rh_refinement_sphere_surface/build")
  file(MAKE_DIRECTORY "${rh_surface_test_dir}")
  add_test(NAME rh_refinement_sphere_surface
    COMMAND $<TARGET_FILE:rh_refinement_example>
      -d pde/surfacepoisson_sphere_bvp.xml -a 0 -b 1958
      -m 1 -i 2 --plot --exactGeo --no-open)
  set_tests_properties(rh_refinement_sphere_surface PROPERTIES
    WORKING_DIRECTORY "${rh_surface_test_dir}"
    LABELS "examples;regression"
    TIMEOUT 60)
endif()

list(LENGTH FILES len)
message(STATUS "Number of examples to compile: ${len}")

