### CMakeLists.txt --- AutoDiff library
## Author: H.M. Verhelst
#
# The library is header-only, the extra code is downloaded in place of
# the source folder of G+Smo
#
######################################################################

# ## AutoDiff extension
# project(gsAutoDiffExtension)

# cmake_minimum_required(VERSION 3.16)

# project(app)

# find_package(autodiff)

# # add_executable(app main.cpp)

# target_link_libraries(app autodiff::autodiff)


######################################################################
## CMakeLists.txt --- gsAutoDiff
## This file is part of the G+Smo library.
##
## Author: H.M. Verhelst
######################################################################

## gsAutoDiff module
project(gsAutoDiff)

# Apply G+Smo config
include(gsConfig)

# Check if AutoDiff was found by dependencies.cmake
if(NOT autodiff_FOUND)
  message(WARNING "AutoDiff not found. Skipping gsAutoDiff module.")
  return()
endif()

# Get autodiff include directories from the target
get_target_property(AUTODIFF_INCLUDE_DIRS autodiff::autodiff INTERFACE_INCLUDE_DIRECTORIES)
if(AUTODIFF_INCLUDE_DIRS)
  # Add autodiff include directory to G+Smo standard include directories
  set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${AUTODIFF_INCLUDE_DIRS}
    CACHE INTERNAL "Gismo include directories" FORCE)
  message(STATUS "Added autodiff include directories: ${AUTODIFF_INCLUDE_DIRS}")
else()
  message(WARNING "Could not get autodiff include directories")
endif()

# Add needed instantiations to the precompiled library
foreach(instance ${GISMO_COEFF_TYPE} ${GISMO_EXTRA_COEFF_TYPE})
  set(GISMO_EXTRA_INSTANCE ${GISMO_EXTRA_INSTANCE}
    "autodiff::detail::Dual<${instance},${instance}>"
    CACHE INTERNAL "Additional instantiations")
endforeach()

## Collect files
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"  )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

#Symlink include dir
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src"
  DESTINATION include/gismo/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.h" )

# add example files
if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(examples)
else()
  add_subdirectory(examples EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

# add unittests
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")

message(STATUS "gsAutoDiff module configured successfully")