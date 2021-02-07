######################################################################
## gsFetch.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

include(CMakeParseArguments)

#latest CMake has FetchContent
function(gismo_fetch_directory)
  set(GF_NAME "${ARGV0}")
  set(oneValueArgs
    # Protect the following options
    DESTINATION
    SOURCE_DIR
    BINARY_DIR
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
    TEST_COMMAND )
  cmake_parse_arguments(GF "${GF_NAME}" "${oneValueArgs}" "" ${ARGN})
  #message( GF_UNPARSED_ARGUMENTS "= ${GF_UNPARSED_ARGUMENTS}")

  file(GLOB RESULT ${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME})
  list(LENGTH RESULT RESULT_LENGTH)
  if(NOT RESULT_LENGTH EQUAL 0)
    message(STATUS "Enabling remote module ${GF_NAME} (${GF_DESTINATION}) - found")
    return()
  endif()

  message(STATUS "Enabling remote module ${GF_NAME} (${GF_DESTINATION})")
  set(GF_SOURCE_DIR   "${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME}")
  set(GF_BINARY_DIR   "${gismo_BINARY_DIR}/${GF_DESTINATION}/${GF_NAME}_fetch")
  set(GF_DOWNLOAD_DIR "${gismo_BINARY_DIR}/${GF_DESTINATION}/${GF_NAME}_fetch")
  set(${GF_NAME}_SOURCE_DIR "${GF_SOURCE_DIR}" PARENT_SCOPE)
  set(${GF_NAME}_BINARY_DIR "${GF_BINARY_DIR}" PARENT_SCOPE)
  file(REMOVE "${GF_DOWNLOAD_DIR}/CMakeCache.txt")

  #  if(NOT EXISTS ${GF_DOWNLOAD_DIR}/CMakeLists.txt)
  file(WRITE ${GF_DOWNLOAD_DIR}/CMakeLists.txt "if(POLICY CMP0048)\ncmake_policy(SET CMP0048 NEW)\nendif()\nif(POLICY CMP0054)\ncmake_policy(SET CMP0054 NEW)\nendif()\ncmake_minimum_required(VERSION 2.8.12)\nproject(${GF_NAME}_fetch NONE)\ninclude(ExternalProject)\nExternalProject_Add(${GF_NAME}_fetch\n ${GF_UNPARSED_ARGUMENTS}\n SOURCE_DIR          \"${GF_SOURCE_DIR}\"\n BINARY_DIR          \"${GF_BINARY_DIR}\"\n CONFIGURE_COMMAND   \"\"\n BUILD_COMMAND       \"\"\n INSTALL_COMMAND     \"\"\n TEST_COMMAND        \"\"\n UPDATE_DISCONNECTED TRUE)\n")
  #  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}"
    -D "CMAKE_MAKE_PROGRAM:FILE=${CMAKE_MAKE_PROGRAM}" .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Configure step for ${GF_NAME} failed: ${result}")
  endif()

  # make sure that directory exists
  cmake_parse_arguments(GF "${GF_NAME}" "SVN_REPOSITORY" "" ${ARGN})
  if(DEFINED GF_SVN_REPOSITORY AND NOT EXISTS "${GF_SOURCE_DIR}/.svn")
    execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} clean
      OUTPUT_QUIET
      WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Build step for ${GF_NAME} failed: ${result}")
  endif()

  file(GLOB RESULT ${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME})
  list(LENGTH RESULT RESULT_LENGTH)
  if(RESULT_LENGTH EQUAL 0)
    message(SEND_ERROR "Enabling remote module ${GF_NAME} (${GF_DESTINATION}) - not found")
  else()
    message(STATUS "Enabling remote module ${GF_NAME} (${GF_DESTINATION}) - downloaded")
  endif()

endfunction()

# called to fetch/download a submodule form git (working) and svn (in progress)
# (ARGV0) SUBMODULE:  name of submodule
# (ARGVN)             Used to give Update Command, unittests/CMakeLists.txt-15-#07.08.18
function(gismo_fetch_module SUBMODULE)

# TODO: online/offline mode

  get_repo_info(GISMO_REPO GISMO_REPO_REV) # or set manually

  #if (NOT DEFINED GISMO_FETCH_PROT)
  #  set(GISMO_FETCH_PROT https) #ssh
  #endif()

  #message("Fetch ${SUBMODULE} (repository: ${GISMO_REPO}, revision: ${GISMO_REPO_REV}, protocol: ${GISMO_FETCH_PROT}, username: ${GISMO_UNAME}, password: ${GISMO_PASS})")

  if("x${GISMO_REPO}" STREQUAL "xgit")
    #if("x${GISMO_FETCH_PROT}" STREQUAL "xssh")
    #  set(git_repo git@github.com:gismo/${SUBMODULE}.git)
    #elseif("x${GISMO_FETCH_PROT}" STREQUAL "xhttps")
    #  set(git_repo https://github.com/gismo/${SUBMODULE}.git)
    #endif()
    # gismo_fetch_directory(${ARGN} GIT_REPOSITORY  ${git_repo} DESTINATION  extensions)
    
    if(NOT EXISTS "${gismo_SOURCE_DIR}/extensions/${SUBMODULE}/CMakeLists.txt")
      message(STATUS "Initializing remote submodule ${SUBMODULE}")
      find_package(Git REQUIRED)

      # init SUBMODULE
      execute_process(COMMAND "${GIT_EXECUTABLE}" "submodule" "update" "--init" "extensions/${SUBMODULE}"
        WORKING_DIRECTORY ${gismo_SOURCE_DIR}
        #RESULT_VARIABLE gresult
        #OUTPUT_QUIET
        )
    endif()

  elseif("x${GISMO_REPO}" STREQUAL "xsvn")
    #if("x${GISMO_FETCH_PROT}" STREQUAL "xssh") message(ERROR "GitHub does not support svn+ssh") endif()
    gismo_fetch_directory(${SUBMODULE}
      SVN_REPOSITORY https://github.com/gismo/${SUBMODULE}/trunk
      SVN_USERNAME ${GISMO_UNAME} # Username for Subversion checkout and update
      SVN_PASSWORD ${GISMO_PASS}  # Password for Subversion checkout and update
      SVN_TRUST_CERT 1            # Trust the Subversion server site certificate
      DESTINATION  extensions
      )
  else()
    gismo_fetch_directory(${SUBMODULE}
      URL https://github.com/gismo/${SUBMODULE}/archive/master.zip
      DESTINATION  extensions
      )
  endif()

  # get list of programs to compile
  if(EXISTS "${gismo_SOURCE_DIR}/extensions/${SUBMODULE}/CMakeLists.txt")
    add_subdirectory(${gismo_SOURCE_DIR}/extensions/${SUBMODULE} ${gismo_BINARY_DIR}/extensions/${SUBMODULE})
    if(EXISTS "${gismo_SOURCE_DIR}/extensions/${SUBMODULE}/filedata")
      string(REGEX MATCH "extensions/${SUBMODULE}/filedata" fmatch ${GISMO_SEARCH_PATHS})
      if(NOT fmatch)
        set(GISMO_SEARCH_PATHS "${GISMO_SEARCH_PATHS};${gismo_SOURCE_DIR}/extensions/${SUBMODULE}/filedata/" CACHE INTERNAL "File search paths")
      endif()
    endif()
  else()
    #WARNING
  endif()

endfunction()
