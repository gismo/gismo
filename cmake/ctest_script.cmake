######################################################################
## ctest_script.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2016 RICAM-Linz.
##
## To execute:
##
##   ctest -S /path/to/ctest_script.cmake
##
## It is recommended to make a copy of the file (especially using git).
## For extra information
##
##   ctest -S /path/to/ctest_script.cmake -V
##
## or even  -VV
##
## Set execution options in the Configuration part.
## For multiple tests (eg. different compilers) make multiple copies
## of this file and adjust options.
##
######################################################################

## #################################################################
## Configuration
## #################################################################

# ID for this computer that shows up on the dashboard.
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
set(CTEST_SITE "${HOSTNAME}")

# The Generator for CMake.
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# Set compiler
#set(CXX g++)
#set(CC  gcc)

# Build type
set(CTEST_BUILD_CONFIGURATION RelWithDebInfo)

# Name of this build
find_program(UNAME NAMES uname)
exec_program("${UNAME}" ARGS "-s" OUTPUT_VARIABLE osname)
exec_program("${UNAME}" ARGS "-m" OUTPUT_VARIABLE "cpu")
set(CTEST_BUILD_NAME "${osname}-${cpu} ${CTEST_CMAKE_GENERATOR} / ${CTEST_BUILD_CONFIGURATION} ${CXX}")

# Test type (Nightly, Continuous, Experimental)
set(dashboard_model "Experimental")

# For continuous builds, number of seconds to stay alive
set(dashboard_runtime 40000)

# Build flags
set( CTEST_BUILD_FLAGS "-j2")

# Source folder
set(CTEST_SOURCE_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/..)

#Build folder
set(CTEST_BINARY_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/build_ctest)

# Update type (eg. svn or git)
set( UPDATE_TYPE git)
set( CTEST_UPDATE_COMMAND "git") #"${CTEST_SOURCE_DIRECTORY}/cmake/svn_github.sh"

# Timeouts
set(CTEST_TEST_TIMEOUT 200 CACHE STRING 
    "Maximum time allowed before CTest will kill the test.") 
set(DART_TESTING_TIMEOUT 200 CACHE STRING 
    "Maximum time allowed before CTest will kill the test." FORCE)

# Build options
set(gismo_build_options 
    -DGISMO_BUILD_LIB=ON
    -DGISMO_WITH_ONURBS=ON
    #-DGISMO_WITH_IPOPT=ON
    #-DIpOpt_DIR=/path/to/ipopt
    #-DGISMO_WITH_PSOLID=ON
    #-DParasolid_DIR=/path/to/parasolid
    #-DGISMO_BUILD_AXL=ON -DAxel_DIR=/path/to/axel
    #-DGISMO_PLAINDOX=ON
    #-DGISMO_BUILD_COVERAGE=ON
)

# Coverage analysis
#set(dashboard_do_coverage TRUE)
#set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")
#set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")

# Memory check with valgrind
#set(dashboard_do_memcheck true)
#set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")



## #################################################################
## Test routines
## #################################################################

macro(run_ctests)
  ctest_configure(OPTIONS "${gismo_build_options}")
  ctest_submit(PARTS Update Notes Configure)
  ctest_build()
  ctest_build(TARGET doc-snippets APPEND)
  ctest_submit(PARTS Build)
  ctest_test()
  ctest_submit(PARTS Test)

  if(dashboard_do_coverage)
     message("Running coverage..")
     ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
     ctest_submit(PARTS Coverage)
  endif()

  if(dashboard_do_memcheck)
    message("Running memcheck..")
    ctest_memcheck()
    ctest_submit(PARTS MemCheck)
  endif()
endmacro(run_ctests)

#if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
#  message("Starting fresh configuration...")
#  write_cache()
#endif()

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${dashboard_model})

if(NOT "${dashboard_model}" STREQUAL "Continuous")

ctest_update()
run_ctests()

else() #continuous model

while(${CTEST_ELAPSED_TIME} LESS ${dashboard_runtime})
  set(START_TIME ${CTEST_ELAPSED_TIME})
  ctest_update(RETURN_VALUE count)
  #message(STATUS "Found ${count} changed files.")
  if( ${count} GREATER 0 )
    run_ctests()
  endif()
  ctest_sleep(${START_TIME} 300    ${CTEST_ELAPSED_TIME})
endwhile()

endif(NOT "${dashboard_model}" STREQUAL "Continuous")
