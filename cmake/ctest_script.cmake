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
## For more verbosity add the flag
##
##   ctest -S /path/to/ctest_script.cmake -V
##
## or even -VV.
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
#set(CTEST_CMAKE_GENERATOR "Ninja")
#set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
#set(CTEST_CMAKE_GENERATOR "NMake Makefiles JOM")
#set(CTEST_CMAKE_GENERATOR "MinGW Makefiles")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 12 2013")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015 Win64")
#set(CTEST_CMAKE_GENERATOR "Xcode")
#set(CTEST_CMAKE_GENERATOR "CodeBlocks")
#set(CTEST_CMAKE_GENERATOR "Sublime Text 2")
#set(CTEST_CMAKE_GENERATOR "Eclipse CDT4")

# Set environment/compiler
#set(ENV{LD_LIBRARY_PATH} /path/to/vendor/lib)
#set(ENV{CC}  "gcc")
#set(ENV{CXX} "g++")
#exec_program(source ARGS "/path/to/iccvars.sh intel64")
#set(ENV{CC}  "icc")
#set(ENV{CXX} "icpc")
#set(ENV{CC}  "clang")
#set(ENV{CXX} "clang++")

# Build type
set(CTEST_BUILD_CONFIGURATION RelWithDebInfo)

# Name of this build
find_program(UNAME NAMES uname)
exec_program("${UNAME}" ARGS "-s" OUTPUT_VARIABLE osname)
exec_program("${UNAME}" ARGS "-m" OUTPUT_VARIABLE "cpu")
set(CTEST_BUILD_NAME "${osname}-${cpu} ${CTEST_CMAKE_GENERATOR}/${CTEST_BUILD_CONFIGURATION} $ENV{CXX}")

# Test type (Nightly, Continuous, Experimental)
set(test_model "Experimental")

# For continuous builds, number of seconds to stay alive
set(test_runtime 40000)

# Build flags
set(CTEST_BUILD_FLAGS "-j2")

# Source folder
set(CTEST_SOURCE_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/..)

#Build folder
set(CTEST_BINARY_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/build_ctest)

# Update type (eg. svn or git)
set(UPDATE_TYPE git)
set(CTEST_UPDATE_COMMAND "git")

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
#set(test_coverage TRUE)
#set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")
#set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")

# Memory check with valgrind
#set(test_memcheck true)
#set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "/path_to/suppression_file.supp")
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --show-reachable=yes --track-origins=yes")


## #################################################################
## Test routines
## #################################################################

macro(run_ctests)
  ctest_configure(OPTIONS "${gismo_build_options}")
  ctest_submit(PARTS Update Notes Configure)
  ctest_build(TARGET unittests APPEND)
  ctest_submit(PARTS Build)
  ctest_test()
  ctest_submit(PARTS Test)

  if(test_coverage)
     message("Running coverage..")
     ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
     ctest_submit(PARTS Coverage)
  endif()

  if(test_memcheck)
    message("Running memcheck..")
    ctest_memcheck()
    ctest_submit(PARTS MemCheck)
  endif()
endmacro(run_ctests)

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${test_model})

if(NOT "${test_model}" STREQUAL "Continuous")

ctest_update()
run_ctests()

else() #continuous model

while(${CTEST_ELAPSED_TIME} LESS ${test_runtime})
  set(START_TIME ${CTEST_ELAPSED_TIME})
  ctest_update(RETURN_VALUE count)
  #message(STATUS "Found ${count} changed files.")
  if( ${count} GREATER 0 )
    run_ctests()
  endif()
  ctest_sleep(${START_TIME} 300 ${CTEST_ELAPSED_TIME})
endwhile()

endif(NOT "${test_model}" STREQUAL "Continuous")
