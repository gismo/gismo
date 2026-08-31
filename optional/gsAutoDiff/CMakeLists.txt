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

# Check if AutoDiff was found by dependencies.cmake; if not, fetch it
if(NOT autodiff_FOUND)
  message(STATUS "AutoDiff not found on system. Fetching from https://github.com/autodiff/autodiff/ ...")
  include(gsFetch)
  gismo_fetch_directory(autodiff
    GIT_REPOSITORY https://github.com/autodiff/autodiff.git
    GIT_TAG        v1.1.2
    DESTINATION    external
  )

  set(AUTODIFF_FETCH_DIR "${gismo_SOURCE_DIR}/external/autodiff")
  if(EXISTS "${AUTODIFF_FETCH_DIR}/autodiff")
    if(NOT TARGET autodiff::autodiff)
      add_library(autodiff::autodiff INTERFACE IMPORTED GLOBAL)
      target_include_directories(autodiff::autodiff INTERFACE "${AUTODIFF_FETCH_DIR}")
    endif()
    set(autodiff_FOUND TRUE)
    set(autodiff_INCLUDE_DIR "${AUTODIFF_FETCH_DIR}")
    message(STATUS "AutoDiff fetched successfully: ${AUTODIFF_FETCH_DIR}")
  else()
    message(WARNING "Failed to fetch AutoDiff. Skipping gsAutoDiff module.")
    return()
  endif()
endif()

# AutoDiff mode options (only available when gsAutoDiff is enabled)
option(GISMO_AUTODIFF_FORWARD "Enable forward mode automatic differentiation (dual_t)" ON)
option(GISMO_AUTODIFF_BACKWARD "Enable backward/reverse mode automatic differentiation (var_t)" OFF)

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

# Add needed instantiations to the precompiled library based on enabled modes
set(AUTODIFF_INSTANCES "")
if(GISMO_AUTODIFF_FORWARD)
  set(AUTODIFF_INSTANCES "gismo::dual_t")
  message(STATUS "gsAutoDiff: Forward mode (dual_t) enabled")
endif()
if(GISMO_AUTODIFF_BACKWARD)
  if(AUTODIFF_INSTANCES)
    set(AUTODIFF_INSTANCES "${AUTODIFF_INSTANCES};gismo::var_t")
  else()
    set(AUTODIFF_INSTANCES "gismo::var_t")
  endif()
  message(STATUS "gsAutoDiff: Backward mode (var_t) enabled")
endif()

if(GISMO_AUTODIFF_FORWARD)
  add_definitions(-DGISMO_AUTODIFF_FORWARD)
endif()
if(GISMO_AUTODIFF_BACKWARD)
  add_definitions(-DGISMO_AUTODIFF_BACKWARD)
endif()

if(AUTODIFF_INSTANCES)
  set(GISMO_EXTRA_INSTANCE ${GISMO_EXTRA_INSTANCE} 
    "${AUTODIFF_INSTANCES}"
    CACHE INTERNAL "Additional instantiations")
endif()

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