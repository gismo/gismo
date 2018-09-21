######################################################################
## gsFetch.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2018
######################################################################

include(CMakeParseArguments)

#latest CMake has FetchContent
function(gismo_fetch_directory)

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

  if( EXISTS "${gismo_SOURCE_DIR}/extensions/${GF_NAME}/CMakeLists.txt")
    message(STATUS "Enabling remote module ${GF_NAME} - found")
    return()
  endif()

  message(STATUS "Enabling remote module ${GF_NAME}")
  set(GF_SOURCE_DIR   "${gismo_SOURCE_DIR}/extensions/${GF_NAME}")
  set(GF_BINARY_DIR   "${gismo_BINARY_DIR}/extensions/${GF_NAME}_fetch")
  set(GF_DOWNLOAD_DIR "${gismo_BINARY_DIR}/extensions/${GF_NAME}_fetch")
  set(${GF_NAME}_SOURCE_DIR "${GF_SOURCE_DIR}" PARENT_SCOPE)
  set(${GF_NAME}_BINARY_DIR "${GF_BINARY_DIR}" PARENT_SCOPE)
  file(REMOVE "${GF_DOWNLOAD_DIR}/CMakeCache.txt")

  #  if(NOT EXISTS ${GF_DOWNLOAD_DIR}/CMakeLists.txt)
  file(WRITE ${GF_DOWNLOAD_DIR}/CMakeLists.txt "if(POLICY CMP0048)\ncmake_policy(SET CMP0048 NEW)\nendif()\nif(POLICY CMP0054)\ncmake_policy(SET CMP0054 NEW)\nendif()\ncmake_minimum_required(VERSION 2.8.8)\nproject(${GF_NAME}_fetch NONE)\ninclude(ExternalProject)\nExternalProject_Add(${GF_NAME}_fetch\n ${GF_UNPARSED_ARGUMENTS}\n SOURCE_DIR          \"${GF_SOURCE_DIR}\"\n BINARY_DIR          \"${GF_BINARY_DIR}\"\n CONFIGURE_COMMAND   \"\"\n BUILD_COMMAND       \"\"\n INSTALL_COMMAND     \"\"\n TEST_COMMAND        \"\")\n")
  #  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}"
    -D "CMAKE_MAKE_PROGRAM:FILE=${CMAKE_MAKE_PROGRAM}" .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Configure step for ${GF_NAME} failed: ${result}")
  endif()

  #! Update step requires the git sources to be available
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Build step for ${GF_NAME} failed: ${result}")
  endif()

  if(NOT EXISTS "${gismo_SOURCE_DIR}/extensions/${GF_NAME}/CMakeLists.txt")
    message(SEND_ERROR "Enabling remote module ${GF_NAME} - not found")
  else()
    message(STATUS "Enabling remote module ${GF_NAME} - downloaded")
  endif()

endfunction()


function(gismo_fetch_module)

  get_repo_info(GISMO_REPO GISMO_REPO_REV) # or set manually
  if (NOT DEFINED GISMO_FETCH_PROT)
    set(GISMO_FETCH_PROT https) #ssh
  endif()

  #message("Fetch ${ARGV0} (repository: ${GISMO_REPO}, revision: ${GISMO_REPO_REV}, protocol: ${GISMO_FETCH_PROT}, username: ${GISMO_UNAME}, password: ${GISMO_PASS})")

  if("x${GISMO_REPO}" STREQUAL "xgit")
    if("x${GISMO_FETCH_PROT}" STREQUAL "xssh")
      set(git_repo git@github.com:gismo/${ARGV0}.git)
    elseif("x${GISMO_FETCH_PROT}" STREQUAL "xhttps")
      set(git_repo https://github.com/gismo/${ARGV0}.git)
    #else()
    endif()
    gismo_fetch_directory(${ARGN}
      GIT_REPOSITORY  ${git_repo})
  elseif("x${GISMO_REPO}" STREQUAL "xsvn")
    #if("x${GISMO_FETCH_PROT}" STREQUAL "xssh") message(ERROR "GitHub does not support svn+ssh") endif()
    gismo_fetch_directory(${ARGN}
      SVN_REPOSITORY https://github.com/gismo/${ARGV0}/trunk
      SVN_USERNAME ${GISMO_UNAME} # Username for Subversion checkout and update
      SVN_PASSWORD ${GISMO_PASS}  # Password for Subversion checkout and update
      SVN_TRUST_CERT 1            # Trust the Subversion server site certificate
      )
  else()
    gismo_fetch_directory(${ARGN}
      URL https://github.com/gismo/${ARGV0}/archive/master.zip)
  endif()

  if(EXISTS "${gismo_SOURCE_DIR}/extensions/${ARGN}/CMakeLists.txt")
    add_subdirectory(${gismo_SOURCE_DIR}/extensions/${ARGN} ${gismo_BINARY_DIR}/extensions/${ARGN})
  endif()

endfunction()
