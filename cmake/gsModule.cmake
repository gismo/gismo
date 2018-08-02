######################################################################
## gsModue.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2018
######################################################################

include(CMakeParseArguments)

#latest CMake has FetchContent
function(gismo_fetch_module)

  set(GF_NAME "${ARGV0}")
  set(oneValueArgs
    # Protect the following options
    SOURCE_DIR
    BINARY_DIR
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
    TEST_COMMAND )
  cmake_parse_arguments(GF "${GF_NAME}" "${oneValueArgs}" "" ${ARGN})
  #message( GF_UNPARSED_ARGUMENTS "= ${GF_UNPARSED_ARGUMENTS}")

  message(STATUS "Enabling remote module ${GF_NAME}")
  set(GF_SOURCE_DIR   "${CMAKE_SOURCE_DIR}/extensions/${GF_NAME}")
  set(GF_BINARY_DIR   "${CMAKE_BINARY_DIR}/extensions/${GF_NAME}_fetch")
  set(GF_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/extensions/${GF_NAME}_fetch")
  set(${GF_NAME}_SOURCE_DIR "${GF_SOURCE_DIR}" PARENT_SCOPE)
  set(${GF_NAME}_BINARY_DIR "${GF_BINARY_DIR}" PARENT_SCOPE)
  file(REMOVE "${GF_DOWNLOAD_DIR}/CMakeCache.txt")

#  if(NOT EXISTS ${GF_DOWNLOAD_DIR}/CMakeLists.txt)
    file(WRITE ${GF_DOWNLOAD_DIR}/CMakeLists.txt "cmake_minimum_required(VERSION 2.8.2)\n project(${GF_NAME}_fetch NONE)\n include(ExternalProject)\n ExternalProject_Add(${GF_NAME}_fetch\n ${GF_UNPARSED_ARGUMENTS}\n SOURCE_DIR          \"${GF_SOURCE_DIR}\"\n BINARY_DIR          \"${GF_BINARY_DIR}\"\n CONFIGURE_COMMAND   \"\"\n BUILD_COMMAND       \"\"\n INSTALL_COMMAND     \"\"\n TEST_COMMAND        \"\")\n")
#  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}"
    -D "CMAKE_MAKE_PROGRAM:FILE=${CMAKE_MAKE_PROGRAM}" .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(FATAL_ERROR "Configure step for ${GF_NAME} failed: ${result}")
  endif()

  #! Update step requires the git sources to be available
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(FATAL_ERROR "Build step for ${GF_NAME} failed: ${result}")
  endif()

  message(STATUS "Enabling remote module ${GF_NAME} - done")
endfunction()

function(gismo_fetch_git_module)

get_repo_info(GISMO_REPO GISMO_REPO_REV)
#message("The repository is ${GISMO_REPO}.")
#message("The revision is ${GISMO_REPO_REV}.")

if( EXISTS "${CMAKE_SOURCE_DIR}/extensions/${ARGV0}/CMakeLists.txt")
  message(STATUS "Found local module ${ARGV0}")
else()

  set(git_pr https) #ssh

  if("x${GISMO_REPO}" STREQUAL "xgit")
    gismo_fetch_module(${ARGN}
      GIT_REPOSITORY  ${git_pr}://git@github.com/gismo/${ARGV0}.git)
  elseif("x${GISMO_REPO}" STREQUAL "xsvn")
    gismo_fetch_module(${ARGN}
      SVN_REPOSITORY  https://github.com/gismo/${ARGV0}/trunk)
  else()
    gismo_fetch_module(${ARGN}
      URL https://github.com/gismo/${ARGV0}/archive/master.zip)
  endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/extensions/${ARGV0}/CMakeLists.txt")
  message(FATAL_ERROR "Module ${ARGV0} was not fetched correctly.")
endif()

endfunction()
